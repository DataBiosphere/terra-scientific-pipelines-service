package bio.terra.pipelines.controller;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertInstanceOf;
import static org.mockito.ArgumentMatchers.any;
import static org.mockito.ArgumentMatchers.eq;
import static org.mockito.Mockito.when;
import static org.springframework.test.web.servlet.request.MockMvcRequestBuilders.get;
import static org.springframework.test.web.servlet.result.MockMvcResultMatchers.content;
import static org.springframework.test.web.servlet.result.MockMvcResultMatchers.status;

import bio.terra.common.iam.SamUser;
import bio.terra.common.iam.SamUserFactory;
import bio.terra.pipelines.app.configuration.external.SamConfiguration;
import bio.terra.pipelines.app.controller.GlobalExceptionHandler;
import bio.terra.pipelines.app.controller.QuotasController;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.db.entities.UserQuota;
import bio.terra.pipelines.db.exception.InvalidPipelineException;
import bio.terra.pipelines.generated.model.ApiQuotaWithDetails;
import bio.terra.pipelines.service.QuotasService;
import bio.terra.pipelines.testutils.MockMvcUtils;
import com.fasterxml.jackson.databind.ObjectMapper;
import jakarta.servlet.http.HttpServletRequest;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.boot.test.autoconfigure.web.servlet.WebMvcTest;
import org.springframework.http.MediaType;
import org.springframework.test.context.ContextConfiguration;
import org.springframework.test.context.bean.override.mockito.MockitoBean;
import org.springframework.test.web.servlet.MockMvc;
import org.springframework.test.web.servlet.MvcResult;

@ContextConfiguration(classes = {QuotasController.class, GlobalExceptionHandler.class})
@WebMvcTest
class QuotasControllerTest {
  @MockitoBean QuotasService quotasServiceMock;
  @MockitoBean SamUserFactory samUserFactoryMock;
  @MockitoBean SamConfiguration samConfiguration;

  @Autowired private MockMvc mockMvc;

  private final SamUser testUser = MockMvcUtils.TEST_SAM_USER;
  private final UserQuota testUserQuota = MockMvcUtils.TEST_USER_QUOTA_1;

  @BeforeEach
  void beforeEach() {
    when(samConfiguration.baseUri()).thenReturn("baseSamUri");
    when(samUserFactoryMock.from(any(HttpServletRequest.class), eq("baseSamUri")))
        .thenReturn(testUser);
    when(quotasServiceMock.getOrCreateQuotaForUserAndPipeline(
            testUser.getSubjectId(), PipelinesEnum.ARRAY_IMPUTATION))
        .thenReturn(testUserQuota);
  }

  @Test
  void getQuotaOk() throws Exception {

    MvcResult result =
        mockMvc
            .perform(get("/api/quotas/v1/" + PipelinesEnum.ARRAY_IMPUTATION.getValue()))
            .andExpect(status().isOk())
            .andExpect(content().contentType(MediaType.APPLICATION_JSON))
            .andReturn();

    ApiQuotaWithDetails response =
        new ObjectMapper()
            .readValue(result.getResponse().getContentAsString(), ApiQuotaWithDetails.class);

    assertEquals(testUserQuota.getQuota(), response.getQuotaLimit());
    assertEquals(testUserQuota.getQuotaConsumed(), response.getQuotaConsumed());
    assertEquals(testUserQuota.getPipelineName().getValue(), response.getPipelineName());
  }

  @Test
  void getQuotasBadPipeline() throws Exception {
    String pipelineName = "bad-pipeline-id";

    mockMvc
        .perform(get("/api/quotas/v1/" + pipelineName))
        .andExpect(status().isBadRequest())
        .andExpect(
            result ->
                assertInstanceOf(InvalidPipelineException.class, result.getResolvedException()));
  }
}
