package bio.terra.pipelines.service;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import bio.terra.pipelines.common.utils.PipelinesEnum;
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
  void getQuotaForUserAndPipeline() {
    // add row to user_quotas table
    UserQuota userQuota =
        createAndSaveUserQuota(TestUtils.TEST_USER_ID_1, PipelinesEnum.IMPUTATION_BEAGLE, 30, 100);

    // call service and assert correct UserQuota is returned
    UserQuota returnedUserQuota =
        quotasService.getQuotaForUserAndPipeline(
            TestUtils.TEST_USER_ID_1, PipelinesEnum.IMPUTATION_BEAGLE);

    assertEquals(userQuota.getUserId(), returnedUserQuota.getUserId());
    assertEquals(userQuota.getPipelineName(), returnedUserQuota.getPipelineName());
    assertEquals(userQuota.getQuotaConsumed(), returnedUserQuota.getQuotaConsumed());
    assertEquals(userQuota.getQuota(), returnedUserQuota.getQuota());
  }

  @Test
  void addRowForNonExistentUserQuota() {
    // assert nothing exists in the user_quotas table
    assertTrue(
        userQuotasRepository
            .findByUserIdAndPipelineName(TestUtils.TEST_USER_ID_1, PipelinesEnum.IMPUTATION_BEAGLE)
            .isEmpty());

    // get default quota for pipeline
    int defaultQuotaForImputationBeaglePipeline =
        pipelineQuotasRepository
            .findByPipelineName(PipelinesEnum.IMPUTATION_BEAGLE)
            .getDefaultQuota();

    // call service with same inputs and a new row should exist in user_quotas table
    UserQuota userQuota =
        quotasService.getQuotaForUserAndPipeline(
            TestUtils.TEST_USER_ID_1, PipelinesEnum.IMPUTATION_BEAGLE);

    assertEquals(TestUtils.TEST_USER_ID_1, userQuota.getUserId());
    assertEquals(PipelinesEnum.IMPUTATION_BEAGLE, userQuota.getPipelineName());
    assertEquals(0, userQuota.getQuotaConsumed());
    assertEquals(defaultQuotaForImputationBeaglePipeline, userQuota.getQuota());

    // assert user + pipeline exists in the user_quotas table
    assertTrue(
        userQuotasRepository
            .findByUserIdAndPipelineName(TestUtils.TEST_USER_ID_1, PipelinesEnum.IMPUTATION_BEAGLE)
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
}