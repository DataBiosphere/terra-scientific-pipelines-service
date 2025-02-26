package bio.terra.pipelines.service;

import static org.junit.jupiter.api.Assertions.*;

import bio.terra.common.exception.InternalServerErrorException;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.common.utils.QuotaUnitsEnum;
import bio.terra.pipelines.db.entities.PipelineQuota;
import bio.terra.pipelines.db.entities.UserQuota;
import bio.terra.pipelines.db.repositories.PipelineQuotasRepository;
import bio.terra.pipelines.db.repositories.UserQuotasRepository;
import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import bio.terra.pipelines.testutils.TestUtils;
import org.junit.jupiter.api.Test;
import org.springframework.beans.factory.annotation.Autowired;

class QuotasServiceTest extends BaseEmbeddedDbTest {
  @Autowired QuotasService quotasService;
  @Autowired UserQuotasRepository userQuotasRepository;
  @Autowired PipelineQuotasRepository pipelineQuotasRepository;

  @Test
  void getPipelineQuota() {
    PipelineQuota pipelineQuota = quotasService.getPipelineQuota(PipelinesEnum.ARRAY_IMPUTATION);

    assertEquals(PipelinesEnum.ARRAY_IMPUTATION, pipelineQuota.getPipelineName());
    assertEquals(10000, pipelineQuota.getDefaultQuota());
    assertEquals(500, pipelineQuota.getMinQuotaConsumed());
    assertEquals(QuotaUnitsEnum.SAMPLES, pipelineQuota.getQuotaUnits());
  }

  @Test
  void getQuotaForUserAndPipeline() {
    // add row to user_quotas table
    createAndSaveUserQuota(TestUtils.TEST_USER_ID_1, PipelinesEnum.ARRAY_IMPUTATION, 30, 100);

    // call service and assert correct UserQuota is returned
    assertTrue(
        quotasService
            .getQuotaForUserAndPipeline(TestUtils.TEST_USER_ID_1, PipelinesEnum.ARRAY_IMPUTATION)
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
        createAndSaveUserQuota(TestUtils.TEST_USER_ID_1, PipelinesEnum.ARRAY_IMPUTATION, 30, 100);

    // call service and assert correct UserQuota is returned
    UserQuota returnedUserQuota =
        quotasService.getOrCreateQuotaForUserAndPipeline(
            TestUtils.TEST_USER_ID_1, PipelinesEnum.ARRAY_IMPUTATION);

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
            .findByUserIdAndPipelineName(TestUtils.TEST_USER_ID_1, PipelinesEnum.ARRAY_IMPUTATION)
            .isEmpty());

    // get default quota for pipeline
    int defaultQuotaForImputationBeaglePipeline =
        pipelineQuotasRepository
            .findByPipelineName(PipelinesEnum.ARRAY_IMPUTATION)
            .getDefaultQuota();

    // call service with same inputs and a new row should exist in user_quotas table
    UserQuota userQuota =
        quotasService.getOrCreateQuotaForUserAndPipeline(
            TestUtils.TEST_USER_ID_1, PipelinesEnum.ARRAY_IMPUTATION);

    assertEquals(TestUtils.TEST_USER_ID_1, userQuota.getUserId());
    assertEquals(PipelinesEnum.ARRAY_IMPUTATION, userQuota.getPipelineName());
    assertEquals(0, userQuota.getQuotaConsumed());
    assertEquals(defaultQuotaForImputationBeaglePipeline, userQuota.getQuota());

    // assert user + pipeline exists in the user_quotas table
    assertTrue(
        userQuotasRepository
            .findByUserIdAndPipelineName(TestUtils.TEST_USER_ID_1, PipelinesEnum.ARRAY_IMPUTATION)
            .isPresent());

    // call service with second user to test for unique constraints
    quotasService.getOrCreateQuotaForUserAndPipeline(
        TestUtils.TEST_USER_ID_2, PipelinesEnum.ARRAY_IMPUTATION);

    // assert user + pipeline exists in the user_quotas table for the second user
    assertTrue(
        userQuotasRepository
            .findByUserIdAndPipelineName(TestUtils.TEST_USER_ID_2, PipelinesEnum.ARRAY_IMPUTATION)
            .isPresent());
  }

  UserQuota createAndSaveUserQuota(
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
        createAndSaveUserQuota(TestUtils.TEST_USER_ID_1, PipelinesEnum.ARRAY_IMPUTATION, 30, 100);

    // call service to update quota consumed
    int newQuotaConsumed = 50;
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
        createAndSaveUserQuota(TestUtils.TEST_USER_ID_1, PipelinesEnum.ARRAY_IMPUTATION, 30, 100);

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
        createAndSaveUserQuota(TestUtils.TEST_USER_ID_1, PipelinesEnum.ARRAY_IMPUTATION, 30, 100);

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
        createAndSaveUserQuota(TestUtils.TEST_USER_ID_1, PipelinesEnum.ARRAY_IMPUTATION, 30, 100);

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
        createAndSaveUserQuota(TestUtils.TEST_USER_ID_1, PipelinesEnum.ARRAY_IMPUTATION, 30, 100);

    // call service to update quota limit
    int newQuotaLimit = 20;
    assertThrows(
        InternalServerErrorException.class,
        () -> quotasService.adminUpdateQuotaLimit(userQuota, newQuotaLimit));
  }
}
