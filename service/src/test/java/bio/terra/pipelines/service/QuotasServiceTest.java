package bio.terra.pipelines.service;

import static org.junit.jupiter.api.Assertions.*;
import static org.mockito.ArgumentMatchers.any;
import static org.mockito.Mockito.*;

import bio.terra.common.exception.BadRequestException;
import bio.terra.common.exception.InternalServerErrorException;
import bio.terra.pipelines.app.configuration.internal.PipelineConfigurations;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.common.utils.QuotaAllocationSourceEnum;
import bio.terra.pipelines.common.utils.QuotaAllocationStatusEnum;
import bio.terra.pipelines.common.utils.QuotaUnitsEnum;
import bio.terra.pipelines.db.entities.QuotaAllocation;
import bio.terra.pipelines.db.entities.UserQuota;
import bio.terra.pipelines.db.repositories.QuotaAllocationsRepository;
import bio.terra.pipelines.db.repositories.UserQuotasRepository;
import bio.terra.pipelines.model.UserQuotaTotals;
import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import bio.terra.pipelines.testutils.TestUtils;
import java.util.ArrayList;
import java.util.List;
import java.util.Optional;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;
import org.junit.jupiter.api.Test;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.dao.DataIntegrityViolationException;

class QuotasServiceTest extends BaseEmbeddedDbTest {
  @Autowired QuotasService quotasService;
  @Autowired UserQuotasRepository userQuotasRepository;
  @Autowired QuotaAllocationsRepository quotaAllocationsRepository;
  @Autowired PipelineConfigurations pipelineConfigurations;

  @Test
  void getQuotaUnitsForPipeline() {
    QuotaUnitsEnum quotaUnits =
        quotasService.getQuotaUnitsForPipeline(PipelinesEnum.ARRAY_IMPUTATION);

    assertEquals(QuotaUnitsEnum.SAMPLES, quotaUnits);
  }

  @Test
  void validateQuotaConfigurationForAllPipelines() {
    // Test quotas defined in pipelines-config.yml (test config)
    PipelineConfigurations.PipelineQuotaConfiguration arrayImputationQuota =
        pipelineConfigurations.getQuotaForPipeline(PipelinesEnum.ARRAY_IMPUTATION);
    assertEquals(2500, arrayImputationQuota.getDefaultQuota());
    assertEquals(500, arrayImputationQuota.getMinQuotaConsumed());
    assertEquals(QuotaUnitsEnum.SAMPLES, arrayImputationQuota.getQuotaUnits());

    PipelineConfigurations.PipelineQuotaConfiguration lowPassImputationQuota =
        pipelineConfigurations.getQuotaForPipeline(PipelinesEnum.LOW_PASS_IMPUTATION);
    assertEquals(100, lowPassImputationQuota.getDefaultQuota());
    assertEquals(10, lowPassImputationQuota.getMinQuotaConsumed());
    assertEquals(QuotaUnitsEnum.SAMPLES, lowPassImputationQuota.getQuotaUnits());
  }

  @Test
  void getQuotaForUserAndPipeline() {
    // add row to user_quotas table
    createAndSaveTestUserQuota(TestUtils.TEST_USER_1_ID, PipelinesEnum.ARRAY_IMPUTATION, 30, 100);

    // call service and assert correct UserQuota is returned
    assertTrue(
        quotasService
            .getQuotaForUserAndPipeline(TestUtils.TEST_USER_1_ID, PipelinesEnum.ARRAY_IMPUTATION)
            .isPresent());
    // assert no UserQuota is returned for a non-existing user
    assertTrue(
        quotasService
            .getQuotaForUserAndPipeline("randomUser", PipelinesEnum.ARRAY_IMPUTATION)
            .isEmpty());
  }

  @Test
  void getOrCreateQuotaForUserAndPipeline() {
    // add row to user_quotas table
    UserQuota userQuota =
        createAndSaveTestUserQuota(
            TestUtils.TEST_USER_1_ID, PipelinesEnum.ARRAY_IMPUTATION, 30, 100);

    // call service and assert correct UserQuota is returned
    UserQuota returnedUserQuota =
        quotasService.getOrCreateQuotaForUserAndPipeline(
            TestUtils.TEST_USER_1_ID, PipelinesEnum.ARRAY_IMPUTATION);

    assertEquals(userQuota.getUserId(), returnedUserQuota.getUserId());
    assertEquals(userQuota.getPipelineName(), returnedUserQuota.getPipelineName());
    assertEquals(userQuota.getQuotaConsumed(), returnedUserQuota.getQuotaConsumed());
    assertEquals(userQuota.getQuota(), returnedUserQuota.getQuota());
  }

  @Test
  void addRowForNonExistentUserQuotaAndSecondUser() {
    // assert nothing exists in the user_quotas table
    assertTrue(
        userQuotasRepository
            .findByUserIdAndPipelineName(TestUtils.TEST_USER_1_ID, PipelinesEnum.ARRAY_IMPUTATION)
            .isEmpty());

    // get default quota for pipeline
    int defaultQuotaForImputationBeaglePipeline =
        quotasService.getPipelineQuota(PipelinesEnum.ARRAY_IMPUTATION).getDefaultQuota();

    // call service with same inputs and a new row should exist in user_quotas table
    UserQuota userQuota =
        quotasService.getOrCreateQuotaForUserAndPipeline(
            TestUtils.TEST_USER_1_ID, PipelinesEnum.ARRAY_IMPUTATION);

    assertEquals(TestUtils.TEST_USER_1_ID, userQuota.getUserId());
    assertEquals(PipelinesEnum.ARRAY_IMPUTATION, userQuota.getPipelineName());
    assertEquals(0, userQuota.getQuotaConsumed());
    assertEquals(defaultQuotaForImputationBeaglePipeline, userQuota.getQuota());

    // assert user + pipeline exists in the user_quotas table
    assertTrue(
        userQuotasRepository
            .findByUserIdAndPipelineName(TestUtils.TEST_USER_1_ID, PipelinesEnum.ARRAY_IMPUTATION)
            .isPresent());

    // call service with second user to test for unique constraints
    quotasService.getOrCreateQuotaForUserAndPipeline(
        TestUtils.TEST_USER_2_ID, PipelinesEnum.ARRAY_IMPUTATION);

    // assert user + pipeline exists in the user_quotas table for the second user
    assertTrue(
        userQuotasRepository
            .findByUserIdAndPipelineName(TestUtils.TEST_USER_2_ID, PipelinesEnum.ARRAY_IMPUTATION)
            .isPresent());
  }

  private UserQuota createAndSaveTestUserQuota(
      String userId, PipelinesEnum pipelineName, int quotaConsumed, int quota) {
    UserQuota userQuota = new UserQuota();
    userQuota.setUserId(userId);
    userQuota.setPipelineName(pipelineName);
    userQuota.setQuotaConsumed(quotaConsumed);
    userQuota.setQuota(quota);
    return userQuotasRepository.save(userQuota);
  }

  @Test
  void updateQuotaConsumed() {
    // add row to user_quotas table
    UserQuota userQuota =
        createAndSaveTestUserQuota(
            TestUtils.TEST_USER_1_ID, PipelinesEnum.ARRAY_IMPUTATION, 30, 100);

    // call service to update quota consumed
    int newQuotaConsumed = 50;
    UserQuota updatedUserQuota = quotasService.updateQuotaConsumed(userQuota, newQuotaConsumed);

    assertEquals(userQuota.getUserId(), updatedUserQuota.getUserId());
    assertEquals(userQuota.getPipelineName(), updatedUserQuota.getPipelineName());
    assertEquals(newQuotaConsumed, updatedUserQuota.getQuotaConsumed());
    assertEquals(userQuota.getQuota(), updatedUserQuota.getQuota());
  }

  @Test
  void updateQuotaConsumedZero() {
    // add row to user_quotas table
    UserQuota userQuota =
        createAndSaveTestUserQuota(
            TestUtils.TEST_USER_1_ID, PipelinesEnum.ARRAY_IMPUTATION, 30, 100);

    // call service to update quota consumed
    int newQuotaConsumed = 0;
    UserQuota updatedUserQuota = quotasService.updateQuotaConsumed(userQuota, newQuotaConsumed);

    assertEquals(userQuota.getUserId(), updatedUserQuota.getUserId());
    assertEquals(userQuota.getPipelineName(), updatedUserQuota.getPipelineName());
    assertEquals(newQuotaConsumed, updatedUserQuota.getQuotaConsumed());
    assertEquals(userQuota.getQuota(), updatedUserQuota.getQuota());
  }

  @Test
  void updateQuotaConsumedLessThanZero() {
    // add row to user_quotas table
    UserQuota userQuota =
        createAndSaveTestUserQuota(
            TestUtils.TEST_USER_1_ID, PipelinesEnum.ARRAY_IMPUTATION, 30, 100);

    // call service to update quota consumed
    int newQuotaConsumed = -1;
    assertThrows(
        InternalServerErrorException.class,
        () -> quotasService.updateQuotaConsumed(userQuota, newQuotaConsumed));
  }

  @Test
  void updateQuotaConsumedGreaterThanQuota() {
    // add row to user_quotas table
    UserQuota userQuota =
        createAndSaveTestUserQuota(
            TestUtils.TEST_USER_1_ID, PipelinesEnum.ARRAY_IMPUTATION, 30, 100);

    // call service to update quota consumed
    int newQuotaConsumed = 130;
    assertThrows(
        InternalServerErrorException.class,
        () -> quotasService.updateQuotaConsumed(userQuota, newQuotaConsumed));
  }

  @Test
  void updateQuotaLimit() {
    // add row to user_quotas table
    UserQuota userQuota =
        createAndSaveTestUserQuota(
            TestUtils.TEST_USER_1_ID, PipelinesEnum.ARRAY_IMPUTATION, 30, 100);

    // call service to update quota limit
    int newQuotaLimit = 150;
    UserQuota updatedUserQuota = quotasService.adminUpdateQuotaLimit(userQuota, newQuotaLimit);

    assertEquals(userQuota.getUserId(), updatedUserQuota.getUserId());
    assertEquals(userQuota.getPipelineName(), updatedUserQuota.getPipelineName());
    assertEquals(userQuota.getQuotaConsumed(), updatedUserQuota.getQuotaConsumed());
    assertEquals(newQuotaLimit, updatedUserQuota.getQuota());
  }

  @Test
  void updateQuotaLimitLessThanQuotaConsumed() {
    // add row to user_quotas table
    UserQuota userQuota =
        createAndSaveTestUserQuota(
            TestUtils.TEST_USER_1_ID, PipelinesEnum.ARRAY_IMPUTATION, 30, 100);

    // call service to update quota limit
    int newQuotaLimit = 20;
    assertThrows(
        InternalServerErrorException.class,
        () -> quotasService.adminUpdateQuotaLimit(userQuota, newQuotaLimit));
  }

  @Test
  void verifySufficientUserQuotaForPipelineRun() {
    // user has 1000 quota, consumed 0, pipeline needs 500
    createAndSaveTestUserQuota(TestUtils.TEST_USER_1_ID, PipelinesEnum.ARRAY_IMPUTATION, 500, 1000);

    // should not throw
    quotasService.validateUserHasEnoughQuota(
        TestUtils.TEST_USER_1_ID, PipelinesEnum.ARRAY_IMPUTATION);
  }

  @Test
  void verifyInsufficientUserQuotaForPipelineRun() {
    // user has 1000 quota, consumed 700, pipeline needs 500
    createAndSaveTestUserQuota(TestUtils.TEST_USER_1_ID, PipelinesEnum.ARRAY_IMPUTATION, 700, 1000);

    // should throw exception
    BadRequestException exception =
        assertThrows(
            BadRequestException.class,
            () ->
                quotasService.validateUserHasEnoughQuota(
                    TestUtils.TEST_USER_1_ID, PipelinesEnum.ARRAY_IMPUTATION));

    // Validate the exception message
    assertEquals(
        "Insufficient quota to run the pipeline. Quota available: 300, Minimum quota required: 500. "
            + "Please visit https://services.terra.bio/pipelines/quotas to purchase additional quota.",
        exception.getMessage());
  }

  /**
   * This test simulates a race condition where multiple concurrent threads try to create the same
   * user quota simultaneously. This verifies that the method handles the race condition correctly
   * and doesn't throw duplicate key violations.
   */
  @Test
  void testConcurrentUserQuotaCreation() throws Exception {
    // Use a unique user ID for this test to avoid conflicts with other tests
    String concurrentUserId = "concurrent-test-user-groot";
    PipelinesEnum pipelineName = PipelinesEnum.ARRAY_IMPUTATION;

    // Verify no quota exists for this user
    assertTrue(
        userQuotasRepository.findByUserIdAndPipelineName(concurrentUserId, pipelineName).isEmpty());

    // Number of concurrent threads to simulate the race condition
    int numberOfThreads = 10;
    ExecutorService executor = Executors.newFixedThreadPool(numberOfThreads);
    CountDownLatch startLatch = new CountDownLatch(1);
    CountDownLatch endLatch = new CountDownLatch(numberOfThreads);

    List<Future<UserQuota>> futures = new ArrayList<>();

    // Create multiple threads that will all try to create the same quota simultaneously
    for (int i = 0; i < numberOfThreads; i++) {
      Future<UserQuota> future =
          executor.submit(
              () -> {
                try {
                  // Wait for all threads to be ready before starting
                  startLatch.await();

                  // All threads call this method at the same time
                  return quotasService.getOrCreateQuotaForUserAndPipeline(
                      concurrentUserId, pipelineName);
                } finally {
                  endLatch.countDown();
                }
              });
      futures.add(future);
    }

    // Release all threads at once to maximize chance of race condition
    startLatch.countDown();

    // Wait for all threads to complete (with timeout)
    assertTrue(endLatch.await(10, TimeUnit.SECONDS), "All threads should complete within timeout");

    // Verify all futures completed successfully and returned valid results
    int expectedDefaultQuota = quotasService.getPipelineQuota(pipelineName).getDefaultQuota();
    List<Long> quotaIds = new ArrayList<>();
    for (Future<UserQuota> future : futures) {
      UserQuota quota = future.get(1, TimeUnit.SECONDS);
      assertNotNull(quota, "UserQuota should not be null");
      assertNotNull(quota.getId(), "UserQuota ID should not be null");
      assertEquals(concurrentUserId, quota.getUserId());
      assertEquals(pipelineName, quota.getPipelineName());
      assertEquals(expectedDefaultQuota, quota.getQuota());
      assertEquals(0, quota.getQuotaConsumed());
      quotaIds.add(quota.getId());
    }

    // Verify that only ONE row was created in the database despite concurrent attempts
    List<UserQuota> dbQuotas = userQuotasRepository.findByUserId(concurrentUserId);
    assertEquals(
        1,
        dbQuotas.size(),
        "Only one quota should exist in database despite concurrent creation attempts");

    // Verify all threads received the SAME entity (same ID) from the database
    long uniqueIdCount = quotaIds.stream().distinct().count();
    assertEquals(
        1, uniqueIdCount, "All threads should have returned the same UserQuota entity (same ID)");

    executor.shutdown();
  }

  /**
   * Unit test to verify IllegalStateException is thrown in the edge case where
   * DataIntegrityViolationException occurs but the retry fetch returns empty. This test has been
   * added to satisfy SonarQube's requirement for coverage of the catch block in
   * getOrCreateQuotaForUserAndPipeline.
   */
  @Test
  void testRaceConditionRetryReturnsEmpty() {
    // Create a mock repository to simulate the edge case
    UserQuotasRepository mockRepo = mock(UserQuotasRepository.class);
    QuotaAllocationsRepository mockAllocationsRepo = mock(QuotaAllocationsRepository.class);
    QuotasService serviceWithMock =
        new QuotasService(mockRepo, mockAllocationsRepo, pipelineConfigurations);

    String userId = "test-user-edge-case-rocket";
    PipelinesEnum pipeline = PipelinesEnum.ARRAY_IMPUTATION;

    // First call: findByUserIdAndPipelineName returns empty (no quota exists)
    // Second call: after catching exception, findByUserIdAndPipelineName still returns empty to
    // simulate the edge case
    when(mockRepo.findByUserIdAndPipelineName(userId, pipeline))
        .thenReturn(Optional.empty())
        .thenReturn(Optional.empty());

    // save() throws DataIntegrityViolationException (simulating race condition)
    when(mockRepo.save(any(UserQuota.class)))
        .thenThrow(new DataIntegrityViolationException("duplicate key"));

    // Execute and verify IllegalStateException is thrown
    IllegalStateException exception =
        assertThrows(
            IllegalStateException.class,
            () -> serviceWithMock.getOrCreateQuotaForUserAndPipeline(userId, pipeline));

    // Verify exception message
    assertEquals(
        "Quota should exist but was not found for user test-user-edge-case-rocket and pipeline ARRAY_IMPUTATION.",
        exception.getMessage());

    // Verify the retry happened (findByUserIdAndPipelineName called twice)
    verify(mockRepo, times(2)).findByUserIdAndPipelineName(userId, pipeline);
    verify(mockRepo, times(1)).save(any(UserQuota.class));
  }

  @Test
  void allocateQuotaForUserAndPipelineOk() {
    String userId = "i-am-groot";
    PipelinesEnum pipeline = PipelinesEnum.LOW_PASS_IMPUTATION;

    QuotaAllocation allocation =
        quotasService.allocateQuotaForUserAndPipeline(
            userId,
            pipeline,
            QuotaAllocationSourceEnum.PRODUCTION_FREE,
            500,
            "Additional free quota for testing");

    assertNotNull(allocation.getId());
    assertEquals(userId, allocation.getUserId());
    assertEquals(pipeline, allocation.getPipelineName());
    assertEquals(QuotaAllocationSourceEnum.PRODUCTION_FREE, allocation.getQuotaSource());
    assertEquals(500, allocation.getQuotaAllocated());
    assertEquals(0, allocation.getQuotaConsumed());
    assertEquals(QuotaAllocationStatusEnum.ACTIVE, allocation.getQuotaStatus());
    assertEquals("Additional free quota for testing", allocation.getComments());
    assertNotNull(allocation.getCreated());
    assertNotNull(allocation.getUpdated());
  }

  @Test
  void allocateQuotaForUserAndPipelineAllowsNullComments() {
    QuotaAllocation allocation =
        quotasService.allocateQuotaForUserAndPipeline(
            "rocket-raccoon",
            PipelinesEnum.SV_IMPUTATION,
            QuotaAllocationSourceEnum.INTERNAL_PAID,
            500,
            null);

    assertNotNull(allocation.getId());
    assertNull(allocation.getComments());
  }

  @Test
  void allocateQuotaAllowsMultipleNonDefaultFreeAllocationsOfSameSource() {
    String userId = "i-am-groot";
    PipelinesEnum pipeline = PipelinesEnum.ARRAY_IMPUTATION;

    // the partial unique index only restricts DEFAULT_FREE, so a user can have multiple
    // EXTERNAL_PAID allocations for the same pipeline
    QuotaAllocation first =
        quotasService.allocateQuotaForUserAndPipeline(
            userId, pipeline, QuotaAllocationSourceEnum.FRIENDS_AND_FAMILY_FREE, 500, null);
    QuotaAllocation second =
        quotasService.allocateQuotaForUserAndPipeline(
            userId, pipeline, QuotaAllocationSourceEnum.EXTERNAL_PAID, 300, "first purchase");

    assertNotEquals(first.getId(), second.getId());

    UserQuotaTotals totals = quotasService.getQuotaTotalsForUserAndPipeline(userId, pipeline);
    assertEquals(800, totals.totalAllocated());
    assertEquals(0, totals.totalConsumed());
  }

  @Test
  void getQuotaTotalsReturnsZeroWhenNoAllocationsExist() {
    UserQuotaTotals totals =
        quotasService.getQuotaTotalsForUserAndPipeline(
            "groot-has-zero-quota", PipelinesEnum.ARRAY_IMPUTATION);

    assertEquals(0, totals.totalAllocated());
    assertEquals(0, totals.totalConsumed());
  }

  @Test
  void getQuotaTotalsIsScopedToUserAndPipeline() {
    String userId = "i-am-groot";
    String otherUserId = "rocket-raccoon";

    quotasService.allocateQuotaForUserAndPipeline(
        userId, PipelinesEnum.ARRAY_IMPUTATION, QuotaAllocationSourceEnum.DEFAULT_FREE, 2500, null);
    // same user, different pipeline - should not be included
    quotasService.allocateQuotaForUserAndPipeline(
        userId,
        PipelinesEnum.LOW_PASS_IMPUTATION,
        QuotaAllocationSourceEnum.DEFAULT_FREE,
        100,
        null);
    // different user, same pipeline - should not be included
    quotasService.allocateQuotaForUserAndPipeline(
        otherUserId,
        PipelinesEnum.ARRAY_IMPUTATION,
        QuotaAllocationSourceEnum.DEFAULT_FREE,
        10000,
        null);

    UserQuotaTotals totals =
        quotasService.getQuotaTotalsForUserAndPipeline(userId, PipelinesEnum.ARRAY_IMPUTATION);
    assertEquals(2500, totals.totalAllocated());
    assertEquals(0, totals.totalConsumed());
  }

  /**
   * This test verifies the partial unique index on (user_id, pipeline_name) WHERE quota_source =
   * 'DEFAULT_FREE' exists and is enforced, and that allocateQuotaForUserAndPipeline surfaces that
   * violation as a BadRequestException.
   */
  @Test
  void allocateQuotaRejectsDuplicateDefaultFree() {
    String userId = "i-am-groot";
    PipelinesEnum pipeline = PipelinesEnum.ARRAY_IMPUTATION;

    // first DEFAULT_FREE allocation should succeed
    QuotaAllocation firstAllocation =
        quotasService.allocateQuotaForUserAndPipeline(
            userId, pipeline, QuotaAllocationSourceEnum.DEFAULT_FREE, 2500, null);
    assertNotNull(firstAllocation.getId());

    // a second DEFAULT_FREE allocation for the same user/pipeline should violate the partial
    // unique index and be translated into a BadRequestException
    BadRequestException exception =
        assertThrows(
            BadRequestException.class,
            () ->
                quotasService.allocateQuotaForUserAndPipeline(
                    userId, pipeline, QuotaAllocationSourceEnum.DEFAULT_FREE, 2500, null));

    assertEquals(
        "User already has a DEFAULT_FREE quota allocation for pipeline ARRAY_IMPUTATION",
        exception.getMessage());

    // a non-DEFAULT_FREE allocation for the same user/pipeline is unaffected by the partial index
    // and should succeed
    QuotaAllocation secondAllocation =
        quotasService.allocateQuotaForUserAndPipeline(
            userId, pipeline, QuotaAllocationSourceEnum.EXTERNAL_PAID, 500, "extra quota");
    assertNotNull(secondAllocation.getId());

    // only the two successful allocations (DEFAULT_FREE + EXTERNAL_PAID) should have persisted;
    // the failed duplicate DEFAULT_FREE insert should not have left a row behind
    UserQuotaTotals totals = quotasService.getQuotaTotalsForUserAndPipeline(userId, pipeline);
    assertEquals(3000, totals.totalAllocated());
    assertEquals(0, totals.totalConsumed());
  }
}
