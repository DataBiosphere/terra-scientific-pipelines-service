package bio.terra.pipelines.db.entities;

import static org.junit.jupiter.api.Assertions.assertEquals;

import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.testutils.BaseTest;
import bio.terra.pipelines.testutils.TestUtils;
import org.junit.jupiter.api.Test;

class UserQuotaTest extends BaseTest {

  @Test
  void testUserQuotaConstructor() {
    UserQuota userQuota =
        new UserQuota(PipelinesEnum.ARRAY_IMPUTATION, TestUtils.TEST_USER_1_ID, 1000, 100);
    assertEquals(PipelinesEnum.ARRAY_IMPUTATION, userQuota.getPipelineName());
    assertEquals(TestUtils.TEST_USER_1_ID, userQuota.getUserId());
    assertEquals(1000, userQuota.getQuota());
    assertEquals(100, userQuota.getQuotaConsumed());
  }
}
