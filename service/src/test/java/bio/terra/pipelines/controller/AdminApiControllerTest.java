package bio.terra.pipelines.controller;

import static bio.terra.pipelines.testutils.MockMvcUtils.TEST_WORKSPACE_UUID;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.mockito.ArgumentMatchers.any;
import static org.mockito.Mockito.*;
import static org.springframework.test.web.servlet.request.MockMvcRequestBuilders.get;
import static org.springframework.test.web.servlet.request.MockMvcRequestBuilders.patch;
import static org.springframework.test.web.servlet.result.MockMvcResultMatchers.content;
import static org.springframework.test.web.servlet.result.MockMvcResultMatchers.status;

import bio.terra.common.exception.ForbiddenException;
import bio.terra.common.iam.BearerTokenFactory;
import bio.terra.common.iam.SamUser;
import bio.terra.common.iam.SamUserFactory;
import bio.terra.pipelines.app.configuration.external.SamConfiguration;
import bio.terra.pipelines.app.controller.AdminApiController;
import bio.terra.pipelines.app.controller.GlobalExceptionHandler;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.dependencies.sam.SamService;
import bio.terra.pipelines.generated.model.ApiAdminPipeline;
import bio.terra.pipelines.generated.model.ApiUpdatePipelineRequestBody;
import bio.terra.pipelines.service.PipelinesService;
import bio.terra.pipelines.testutils.MockMvcUtils;
import com.fasterxml.jackson.core.JsonProcessingException;
import com.fasterxml.jackson.databind.ObjectMapper;
import jakarta.servlet.http.HttpServletRequest;
import java.util.UUID;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.boot.test.autoconfigure.web.servlet.WebMvcTest;
import org.springframework.boot.test.mock.mockito.MockBean;
import org.springframework.http.MediaType;
import org.springframework.test.context.ContextConfiguration;
import org.springframework.test.web.servlet.MockMvc;
import org.springframework.test.web.servlet.MvcResult;

@ContextConfiguration(classes = {AdminApiController.class, GlobalExceptionHandler.class})
@WebMvcTest()
class AdminApiControllerTest {
  @MockBean PipelinesService pipelinesServiceMock;
  @MockBean SamUserFactory samUserFactoryMock;
  @MockBean BearerTokenFactory bearerTokenFactory;
  @MockBean SamConfiguration samConfiguration;
  @MockBean SamService samServiceMock;

  @Autowired private MockMvc mockMvc;
  private final SamUser testUser = MockMvcUtils.TEST_SAM_USER;
  private final String testUserId = testUser.getSubjectId();

  @BeforeEach
  void beforeEach() {
    when(samUserFactoryMock.from(any(HttpServletRequest.class), any())).thenReturn(testUser);
    doNothing().when(samServiceMock).checkAdminAuthz(testUser);
  }

  @Test
  void updatePipelineWorkspaceIdOk() throws Exception {
    when(pipelinesServiceMock.updatePipelineWorkspaceId(
            PipelinesEnum.IMPUTATION_MINIMAC4, TEST_WORKSPACE_UUID))
        .thenReturn(MockMvcUtils.getTestPipeline());
    MvcResult result =
        mockMvc
            .perform(
                patch(
                        String.format(
                            "/api/admin/v1/pipeline/%s",
                            PipelinesEnum.IMPUTATION_MINIMAC4.getValue()))
                    .contentType(MediaType.APPLICATION_JSON)
                    .content(createTestJobPostBody(TEST_WORKSPACE_UUID)))
            .andExpect(status().isOk())
            .andExpect(content().contentType(MediaType.APPLICATION_JSON))
            .andReturn();

    ApiAdminPipeline response =
        new ObjectMapper()
            .readValue(result.getResponse().getContentAsString(), ApiAdminPipeline.class);

    // this is all mocked data so really not worth checking values, really just testing that it's a
    // 200 status with a properly formatted response
    assertEquals(TEST_WORKSPACE_UUID, response.getWorkspaceId());
  }

  @Test
  void updatePipelineWorkspaceIdRequireWorkspaceId() throws Exception {
    mockMvc
        .perform(
            patch(
                    String.format(
                        "/api/admin/v1/pipeline/%s", PipelinesEnum.IMPUTATION_MINIMAC4.getValue()))
                .contentType(MediaType.APPLICATION_JSON))
        .andExpect(status().isBadRequest());
  }

  @Test
  void updatePipelineWorkspaceIdNotAdminUser() throws Exception {
    doThrow(new ForbiddenException("error string")).when(samServiceMock).checkAdminAuthz(testUser);

    mockMvc
        .perform(
            patch(
                    String.format(
                        "/api/admin/v1/pipeline/%s", PipelinesEnum.IMPUTATION_MINIMAC4.getValue()))
                .contentType(MediaType.APPLICATION_JSON)
                .content(createTestJobPostBody(TEST_WORKSPACE_UUID)))
        .andExpect(status().isForbidden());
  }

  @Test
  void getAdminPipelineOk() throws Exception {
    when(pipelinesServiceMock.getPipeline(PipelinesEnum.IMPUTATION_MINIMAC4))
        .thenReturn(MockMvcUtils.getTestPipeline());
    MvcResult result =
        mockMvc
            .perform(
                get(
                    String.format(
                        "/api/admin/v1/pipeline/%s", PipelinesEnum.IMPUTATION_MINIMAC4.getValue())))
            .andExpect(status().isOk())
            .andExpect(content().contentType(MediaType.APPLICATION_JSON))
            .andReturn();

    ApiAdminPipeline response =
        new ObjectMapper()
            .readValue(result.getResponse().getContentAsString(), ApiAdminPipeline.class);

    // this is all mocked data so really not worth checking values, really just testing that it's a
    // 200 status with a properly formatted response
    assertEquals("pipelineName", response.getPipelineName());
  }

  @Test
  void getAdminPipelineNotAdminUser() throws Exception {
    doThrow(new ForbiddenException("error string")).when(samServiceMock).checkAdminAuthz(testUser);

    mockMvc
        .perform(
            get(
                String.format(
                    "/api/admin/v1/pipeline/%s", PipelinesEnum.IMPUTATION_MINIMAC4.getValue())))
        .andExpect(status().isForbidden());
  }

  private String createTestJobPostBody(UUID workspaceId) throws JsonProcessingException {
    ApiUpdatePipelineRequestBody apiUpdatePipelineRequestBody =
        new ApiUpdatePipelineRequestBody().workspaceId(workspaceId);
    return MockMvcUtils.convertToJsonString(apiUpdatePipelineRequestBody);
  }
}
