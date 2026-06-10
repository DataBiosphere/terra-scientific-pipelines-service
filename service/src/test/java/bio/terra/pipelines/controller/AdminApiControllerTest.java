package bio.terra.pipelines.controller;

import static bio.terra.pipelines.testutils.MockMvcUtils.*;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.mockito.ArgumentMatchers.any;
import static org.mockito.Mockito.*;
import static org.springframework.test.web.servlet.request.MockMvcRequestBuilders.get;
import static org.springframework.test.web.servlet.request.MockMvcRequestBuilders.patch;
import static org.springframework.test.web.servlet.result.MockMvcResultMatchers.content;
import static org.springframework.test.web.servlet.result.MockMvcResultMatchers.status;

import bio.terra.common.exception.ForbiddenException;
import bio.terra.common.exception.NotFoundException;
import bio.terra.common.iam.BearerTokenFactory;
import bio.terra.common.iam.SamUser;
import bio.terra.common.iam.SamUserFactory;
import bio.terra.pipelines.app.configuration.external.SamConfiguration;
import bio.terra.pipelines.app.controller.AdminApiController;
import bio.terra.pipelines.app.controller.GlobalExceptionHandler;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.db.entities.UserQuota;
import bio.terra.pipelines.dependencies.sam.SamService;
import bio.terra.pipelines.generated.model.ApiAdminPipeline;
import bio.terra.pipelines.generated.model.ApiAdminQuota;
import bio.terra.pipelines.generated.model.ApiAdminQuotaV2;
import bio.terra.pipelines.generated.model.ApiUpdatePipelineRequestBody;
import bio.terra.pipelines.generated.model.ApiUpdateQuotaLimitRequestBody;
import bio.terra.pipelines.model.Pipeline;
import bio.terra.pipelines.notifications.NotificationService;
import bio.terra.pipelines.service.PipelinesService;
import bio.terra.pipelines.service.QuotasService;
import bio.terra.pipelines.testutils.MockMvcUtils;
import bio.terra.pipelines.testutils.TestUtils;
import com.fasterxml.jackson.core.JsonProcessingException;
import com.fasterxml.jackson.databind.ObjectMapper;
import jakarta.servlet.http.HttpServletRequest;
import java.util.Optional;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Nested;
import org.junit.jupiter.api.Test;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.boot.test.autoconfigure.web.servlet.WebMvcTest;
import org.springframework.http.MediaType;
import org.springframework.test.context.ContextConfiguration;
import org.springframework.test.context.bean.override.mockito.MockitoBean;
import org.springframework.test.web.servlet.MockMvc;
import org.springframework.test.web.servlet.MvcResult;

@ContextConfiguration(classes = {AdminApiController.class, GlobalExceptionHandler.class})
@WebMvcTest()
class AdminApiControllerTest {
  @MockitoBean PipelinesService pipelinesServiceMock;
  @MockitoBean QuotasService quotasServiceMock;
  @MockitoBean SamUserFactory samUserFactoryMock;
  @MockitoBean BearerTokenFactory bearerTokenFactory;
  @MockitoBean SamConfiguration samConfiguration;
  @MockitoBean SamService samServiceMock;
  @MockitoBean NotificationService notificationService;

  @Autowired private MockMvc mockMvc;
  private final SamUser testUser = MockMvcUtils.TEST_SAM_USER;

  @BeforeEach
  void beforeEach() {
    when(samConfiguration.baseUri()).thenReturn("baseSamUri");
    when(samUserFactoryMock.from(any(HttpServletRequest.class), eq("baseSamUri")))
        .thenReturn(testUser);
    doNothing().when(samServiceMock).checkAdminAuthz(testUser);
    when(samServiceMock.isAdmin(testUser)).thenReturn(true);
  }

  @Test
  void updatePipelineWorkspaceOk() throws Exception {
    when(pipelinesServiceMock.adminUpdatePipelineWorkspace(
            PipelinesEnum.ARRAY_IMPUTATION,
            TestUtils.TEST_PIPELINE_VERSION_1,
            null,
            TEST_WORKSPACE_BILLING_PROJECT,
            TEST_WORKSPACE_NAME,
            TEST_TOOL_VERSION))
        .thenReturn(MockMvcUtils.getTestPipeline());
    MvcResult result =
        mockMvc
            .perform(
                patch(
                        String.format(
                            "/api/admin/v1/pipelines/%s/%s",
                            PipelinesEnum.ARRAY_IMPUTATION.getLowerCaseValue(),
                            TestUtils.TEST_PIPELINE_VERSION_1))
                    .contentType(MediaType.APPLICATION_JSON)
                    .content(
                        createTestJobPostBody(
                            TEST_WORKSPACE_BILLING_PROJECT,
                            TEST_WORKSPACE_NAME,
                            TEST_TOOL_VERSION,
                            null)))
            .andExpect(status().isOk())
            .andExpect(content().contentType(MediaType.APPLICATION_JSON))
            .andReturn();

    ApiAdminPipeline response =
        new ObjectMapper()
            .readValue(result.getResponse().getContentAsString(), ApiAdminPipeline.class);

    // this is all mocked data so really not worth checking values, really just testing that it's a
    // 200 status with a properly formatted response
    assertEquals(TEST_WORKSPACE_BILLING_PROJECT, response.getWorkspaceBillingProject());
    assertEquals(TEST_WORKSPACE_NAME, response.getWorkspaceName());
    assertEquals(
        TEST_WORKSPACE_STORAGE_CONTAINER_NAME, response.getWorkspaceStorageContainerName());
    assertEquals(TEST_WORKSPACE_GOOGLE_PROJECT, response.getWorkspaceGoogleProject());
    assertEquals(TEST_TOOL_VERSION, response.getToolVersion());
    assertEquals(false, response.isIsHidden());
    assertNotNull(response.getUpdated());
  }

  @Test
  void updatePipelineWorkspaceOkWithHidden() throws Exception {
    Pipeline hiddenPipeline = getTestPipeline().toBuilder().hidden(true).build();
    when(pipelinesServiceMock.adminUpdatePipelineWorkspace(
            PipelinesEnum.ARRAY_IMPUTATION,
            TestUtils.TEST_PIPELINE_VERSION_1,
            true,
            TEST_WORKSPACE_BILLING_PROJECT,
            TEST_WORKSPACE_NAME,
            TEST_TOOL_VERSION))
        .thenReturn(hiddenPipeline);
    MvcResult result =
        mockMvc
            .perform(
                patch(
                        String.format(
                            "/api/admin/v1/pipelines/%s/%s",
                            PipelinesEnum.ARRAY_IMPUTATION.getLowerCaseValue(),
                            TestUtils.TEST_PIPELINE_VERSION_1))
                    .contentType(MediaType.APPLICATION_JSON)
                    .content(
                        createTestJobPostBody(
                            TEST_WORKSPACE_BILLING_PROJECT,
                            TEST_WORKSPACE_NAME,
                            TEST_TOOL_VERSION,
                            true)))
            .andExpect(status().isOk())
            .andExpect(content().contentType(MediaType.APPLICATION_JSON))
            .andReturn();

    ApiAdminPipeline response =
        new ObjectMapper()
            .readValue(result.getResponse().getContentAsString(), ApiAdminPipeline.class);

    // this is all mocked data so really not worth checking values, really just testing that it's a
    // 200 status with a properly formatted response
    assertEquals(TEST_WORKSPACE_BILLING_PROJECT, response.getWorkspaceBillingProject());
    assertEquals(TEST_WORKSPACE_NAME, response.getWorkspaceName());
    assertEquals(
        TEST_WORKSPACE_STORAGE_CONTAINER_NAME, response.getWorkspaceStorageContainerName());
    assertEquals(TEST_WORKSPACE_GOOGLE_PROJECT, response.getWorkspaceGoogleProject());
    assertEquals(TEST_TOOL_VERSION, response.getToolVersion());
    assertTrue(response.isIsHidden());
    assertNotNull(response.getUpdated());
  }

  @Test
  void updatePipelineWorkspaceIdRequireWorkspaceName() throws Exception {
    mockMvc
        .perform(
            patch(
                    String.format(
                        "/api/admin/v1/pipelines/%s/%s",
                        PipelinesEnum.ARRAY_IMPUTATION.getLowerCaseValue(),
                        TestUtils.TEST_PIPELINE_VERSION_1))
                .contentType(MediaType.APPLICATION_JSON)
                .content(
                    createTestJobPostBody(
                        TEST_WORKSPACE_BILLING_PROJECT, null, TEST_TOOL_VERSION, null)))
        .andExpect(status().isBadRequest());
  }

  @Test
  void updatePipelineWorkspaceIdRequireWorkspaceProject() throws Exception {
    mockMvc
        .perform(
            patch(
                    String.format(
                        "/api/admin/v1/pipelines/%s/%s",
                        PipelinesEnum.ARRAY_IMPUTATION.getLowerCaseValue(),
                        TestUtils.TEST_PIPELINE_VERSION_1))
                .contentType(MediaType.APPLICATION_JSON)
                .content(createTestJobPostBody(null, TEST_WORKSPACE_NAME, TEST_TOOL_VERSION, null)))
        .andExpect(status().isBadRequest());
  }

  @Test
  void updatePipelineWorkspaceIdRequireToolVersion() throws Exception {
    mockMvc
        .perform(
            patch(
                    String.format(
                        "/api/admin/v1/pipelines/%s/%s",
                        PipelinesEnum.ARRAY_IMPUTATION.getLowerCaseValue(),
                        TestUtils.TEST_PIPELINE_VERSION_1))
                .contentType(MediaType.APPLICATION_JSON)
                .content(
                    createTestJobPostBody(
                        TEST_WORKSPACE_BILLING_PROJECT, TEST_WORKSPACE_NAME, null, null)))
        .andExpect(status().isBadRequest());
  }

  @Test
  void updatePipelineWorkspaceIdNotAdminUser() throws Exception {
    doThrow(new ForbiddenException("error string")).when(samServiceMock).checkAdminAuthz(testUser);

    mockMvc
        .perform(
            patch(
                    String.format(
                        "/api/admin/v1/pipelines/%s/%s",
                        PipelinesEnum.ARRAY_IMPUTATION.getLowerCaseValue(),
                        TestUtils.TEST_PIPELINE_VERSION_1))
                .contentType(MediaType.APPLICATION_JSON)
                .content(
                    createTestJobPostBody(
                        TEST_WORKSPACE_BILLING_PROJECT,
                        TEST_WORKSPACE_NAME,
                        TEST_TOOL_VERSION,
                        null)))
        .andExpect(status().isForbidden());
  }

  @Test
  void updatePipelineWorkspaceIdBadVersion() throws Exception {
    when(pipelinesServiceMock.adminUpdatePipelineWorkspace(
            PipelinesEnum.ARRAY_IMPUTATION,
            TestUtils.TEST_PIPELINE_VERSION_1,
            null,
            TEST_WORKSPACE_BILLING_PROJECT,
            TEST_WORKSPACE_NAME,
            TEST_TOOL_VERSION))
        .thenThrow(new NotFoundException("badversion"));

    mockMvc
        .perform(
            patch(
                    String.format(
                        "/api/admin/v1/pipelines/%s/%s",
                        PipelinesEnum.ARRAY_IMPUTATION.getLowerCaseValue(),
                        TestUtils.TEST_PIPELINE_VERSION_1))
                .contentType(MediaType.APPLICATION_JSON)
                .content(
                    createTestJobPostBody(
                        TEST_WORKSPACE_BILLING_PROJECT,
                        TEST_WORKSPACE_NAME,
                        TEST_TOOL_VERSION,
                        null)))
        .andExpect(status().isNotFound());
  }

  @Test
  void getAdminPipelineOk() throws Exception {
    when(pipelinesServiceMock.getPipeline(
            PipelinesEnum.ARRAY_IMPUTATION, TestUtils.TEST_PIPELINE_VERSION_1, true))
        .thenReturn(MockMvcUtils.getTestPipeline());
    MvcResult result =
        mockMvc
            .perform(
                get(
                    String.format(
                        "/api/admin/v1/pipelines/%s/%s",
                        PipelinesEnum.ARRAY_IMPUTATION.getLowerCaseValue(),
                        TestUtils.TEST_PIPELINE_VERSION_1)))
            .andExpect(status().isOk())
            .andExpect(content().contentType(MediaType.APPLICATION_JSON))
            .andReturn();

    ApiAdminPipeline response =
        new ObjectMapper()
            .readValue(result.getResponse().getContentAsString(), ApiAdminPipeline.class);

    // this is all mocked data so really not worth checking values, really just testing that it's a
    // 200 status with a properly formatted response
    assertEquals(PipelinesEnum.ARRAY_IMPUTATION.getLowerCaseValue(), response.getPipelineName());
  }

  @Test
  void getAdminPipelineBadVersion() throws Exception {
    when(pipelinesServiceMock.getPipeline(
            PipelinesEnum.ARRAY_IMPUTATION, TestUtils.TEST_PIPELINE_VERSION_1, true))
        .thenThrow(new NotFoundException("badversion"));

    mockMvc
        .perform(
            get(
                String.format(
                    "/api/admin/v1/pipelines/%s/%s",
                    PipelinesEnum.ARRAY_IMPUTATION.getLowerCaseValue(),
                    TestUtils.TEST_PIPELINE_VERSION_1)))
        .andExpect(status().isNotFound());
  }

  @Test
  void getAdminPipelineNotAdminUser() throws Exception {
    doThrow(new ForbiddenException("error string")).when(samServiceMock).checkAdminAuthz(testUser);

    mockMvc
        .perform(
            get(
                String.format(
                    "/api/admin/v1/pipelines/%s/%s",
                    PipelinesEnum.ARRAY_IMPUTATION.getLowerCaseValue(),
                    TestUtils.TEST_PIPELINE_VERSION_1)))
        .andExpect(status().isForbidden());
  }

  @Nested
  @Deprecated
  @DisplayName("getQuotaForPipelineAndUser V1 tests")
  class GetAdminQuotaV1Tests {

    @Test
    void getUserQuotaOk() throws Exception {
      when(quotasServiceMock.getQuotaForUserAndPipeline(
              TEST_SAM_USER.getSubjectId(), PipelinesEnum.ARRAY_IMPUTATION))
          .thenReturn(Optional.of(TEST_USER_QUOTA_1));
      MvcResult result =
          mockMvc
              .perform(
                  get(
                      String.format(
                          "/api/admin/v1/quotas/%s/%s",
                          PipelinesEnum.ARRAY_IMPUTATION.getLowerCaseValue(),
                          TEST_SAM_USER.getSubjectId())))
              .andExpect(status().isOk())
              .andExpect(content().contentType(MediaType.APPLICATION_JSON))
              .andReturn();

      ApiAdminQuota response =
          new ObjectMapper()
              .readValue(result.getResponse().getContentAsString(), ApiAdminQuota.class);

      // this is all mocked data so really not worth checking values, really just testing that it's
      // a 200 status with a properly formatted response
      assertEquals(PipelinesEnum.ARRAY_IMPUTATION.getLowerCaseValue(), response.getPipelineName());
      assertEquals(TEST_USER_QUOTA_1.getUserId(), response.getUserId());
      assertEquals(TEST_USER_QUOTA_1.getQuota(), response.getQuotaLimit());
      assertEquals(TEST_USER_QUOTA_1.getQuotaConsumed(), response.getQuotaConsumed());
    }

    @Test
    void getAdminQuotaNotAdminUser() throws Exception {
      doThrow(new ForbiddenException("error string"))
          .when(samServiceMock)
          .checkAdminAuthz(testUser);

      mockMvc
          .perform(
              get(
                  String.format(
                      "/api/admin/v1/quotas/%s/%s",
                      PipelinesEnum.ARRAY_IMPUTATION.getLowerCaseValue(),
                      TEST_SAM_USER.getSubjectId())))
          .andExpect(status().isForbidden());
    }

    @Test
    void getAdminQuotaUserQuotaDoesntExist() throws Exception {
      when(quotasServiceMock.getQuotaForUserAndPipeline(
              TEST_SAM_USER.getSubjectId(), PipelinesEnum.ARRAY_IMPUTATION))
          .thenReturn(Optional.empty());
      mockMvc
          .perform(
              get(
                  String.format(
                      "/api/admin/v1/quotas/%s/%s",
                      PipelinesEnum.ARRAY_IMPUTATION.getLowerCaseValue(),
                      TEST_SAM_USER.getSubjectId())))
          .andExpect(status().isBadRequest());
    }
  }

  @Nested
  @Deprecated
  @DisplayName("updateQuotaLimitForPipelineAndUser V1 tests")
  class UpdateAdminQuotaV1Tests {

    @Test
    void updateAdminQuotaOk() throws Exception {
      UserQuota updatedUserQuota =
          new UserQuota(PipelinesEnum.ARRAY_IMPUTATION, TEST_SAM_USER.getSubjectId(), 800, 0);
      when(quotasServiceMock.getQuotaForUserAndPipeline(
              TEST_SAM_USER.getSubjectId(), PipelinesEnum.ARRAY_IMPUTATION))
          .thenReturn(Optional.of(TEST_USER_QUOTA_1));
      when(quotasServiceMock.adminUpdateQuotaLimit(TEST_USER_QUOTA_1, 800))
          .thenReturn(updatedUserQuota);
      when(pipelinesServiceMock.getLatestPipeline(PipelinesEnum.ARRAY_IMPUTATION))
          .thenReturn(MockMvcUtils.getTestPipeline());

      MvcResult result =
          mockMvc
              .perform(
                  patch(
                          String.format(
                              "/api/admin/v1/quotas/%s/%s",
                              PipelinesEnum.ARRAY_IMPUTATION.getLowerCaseValue(),
                              TEST_SAM_USER.getSubjectId()))
                      .contentType(MediaType.APPLICATION_JSON)
                      .content(createTestJobPostBody(800)))
              .andExpect(status().isOk())
              .andExpect(content().contentType(MediaType.APPLICATION_JSON))
              .andReturn();

      ApiAdminQuota response =
          new ObjectMapper()
              .readValue(result.getResponse().getContentAsString(), ApiAdminQuota.class);

      // this is all mocked data so really not worth checking values, really just testing that it's
      // a 200 status with a properly formatted response
      assertEquals(PipelinesEnum.ARRAY_IMPUTATION.getLowerCaseValue(), response.getPipelineName());
      assertEquals(TEST_SAM_USER.getSubjectId(), response.getUserId());
      assertEquals(800, response.getQuotaLimit());

      // assert that NotificationService was called with right parameters
      verify(notificationService)
          .configureAndSendUserQuotaChangedNotification(
              TEST_SAM_USER.getSubjectId(), "displayName", 1000, 800, 10, 790);
    }

    @Test
    void updateAdminQuotaUserQuotaDoesntExist() throws Exception {
      when(quotasServiceMock.getQuotaForUserAndPipeline(
              TEST_SAM_USER.getSubjectId(), PipelinesEnum.ARRAY_IMPUTATION))
          .thenReturn(Optional.empty());
      mockMvc
          .perform(
              patch(
                      String.format(
                          "/api/admin/v1/quotas/%s/%s",
                          PipelinesEnum.ARRAY_IMPUTATION.getLowerCaseValue(),
                          TEST_SAM_USER.getSubjectId()))
                  .contentType(MediaType.APPLICATION_JSON)
                  .content(createTestJobPostBody(800)))
          .andExpect(status().isBadRequest());

      // assert NotificationService is never called
      verify(notificationService, never())
          .configureAndSendUserQuotaChangedNotification(
              any(), any(), anyInt(), anyInt(), anyInt(), anyInt());
    }

    @Test
    void updateAdminQuotaRequireQuotaLimit() throws Exception {
      mockMvc
          .perform(
              patch(
                      String.format(
                          "/api/admin/v1/quotas/%s/%s",
                          PipelinesEnum.ARRAY_IMPUTATION.getLowerCaseValue(),
                          TEST_SAM_USER.getSubjectId()))
                  .contentType(MediaType.APPLICATION_JSON)
                  .content("{}"))
          .andExpect(status().isBadRequest());

      // assert NotificationService is never called
      verify(notificationService, never())
          .configureAndSendUserQuotaChangedNotification(
              any(), any(), anyInt(), anyInt(), anyInt(), anyInt());
    }

    @Test
    void updateAdminQuotaIdNotAdminUser() throws Exception {
      doThrow(new ForbiddenException("error string"))
          .when(samServiceMock)
          .checkAdminAuthz(testUser);

      mockMvc
          .perform(
              patch(
                      String.format(
                          "/api/admin/v1/quotas/%s/%s",
                          PipelinesEnum.ARRAY_IMPUTATION.getLowerCaseValue(),
                          TEST_SAM_USER.getSubjectId()))
                  .contentType(MediaType.APPLICATION_JSON)
                  .content(createTestJobPostBody(500)))
          .andExpect(status().isForbidden());

      // assert NotificationService is never called
      verify(notificationService, never())
          .configureAndSendUserQuotaChangedNotification(
              any(), any(), anyInt(), anyInt(), anyInt(), anyInt());
    }
  }

  @Nested
  @DisplayName("getQuotaForPipelineAndUser V2 tests")
  class GetAdminQuotaV2Tests {

    @Test
    void getUserQuotaOk() throws Exception {
      String userEmail = TEST_SAM_USER.getEmail();
      when(samServiceMock.getUserIdFromEmail(testUser, userEmail))
          .thenReturn(TEST_SAM_USER.getSubjectId());
      when(quotasServiceMock.getOrCreateQuotaForUserAndPipeline(
              TEST_SAM_USER.getSubjectId(), PipelinesEnum.ARRAY_IMPUTATION))
          .thenReturn(TEST_USER_QUOTA_1);
      MvcResult result =
          mockMvc
              .perform(
                  get(
                      String.format(
                          "/api/admin/v2/quotas/%s/%s",
                          PipelinesEnum.ARRAY_IMPUTATION.getLowerCaseValue(), userEmail)))
              .andExpect(status().isOk())
              .andExpect(content().contentType(MediaType.APPLICATION_JSON))
              .andReturn();

      ApiAdminQuotaV2 response =
          new ObjectMapper()
              .readValue(result.getResponse().getContentAsString(), ApiAdminQuotaV2.class);

      assertEquals(PipelinesEnum.ARRAY_IMPUTATION.getLowerCaseValue(), response.getPipelineName());
      assertEquals(TEST_SAM_USER.getSubjectId(), response.getUserId());
      assertEquals(TEST_SAM_USER.getEmail(), response.getUserEmail());
      assertEquals(TEST_USER_QUOTA_1.getQuota(), response.getQuotaLimit());
      assertEquals(TEST_USER_QUOTA_1.getQuotaConsumed(), response.getQuotaConsumed());
    }

    @Test
    void getAdminQuotaNotAdminUser() throws Exception {
      String userEmail = TEST_SAM_USER.getEmail();
      doThrow(new ForbiddenException("Not an admin error"))
          .when(samServiceMock)
          .checkAdminAuthz(testUser);

      mockMvc
          .perform(
              get(
                  String.format(
                      "/api/admin/v2/quotas/%s/%s",
                      PipelinesEnum.ARRAY_IMPUTATION.getLowerCaseValue(), userEmail)))
          .andExpect(status().isForbidden());
    }

    @Test
    void getAdminQuotaUserNotFoundInSam() throws Exception {
      String userEmail = "nonexistent@example.com";
      when(samServiceMock.getUserIdFromEmail(testUser, userEmail))
          .thenThrow(new NotFoundException("User not found in SAM"));

      mockMvc
          .perform(
              get(
                  String.format(
                      "/api/admin/v2/quotas/%s/%s",
                      PipelinesEnum.ARRAY_IMPUTATION.getLowerCaseValue(), userEmail)))
          .andExpect(status().isNotFound());
    }

    @Test
    void getAdminQuotaUserQuotaDoesntExist() throws Exception {
      String userEmail = TEST_SAM_USER.getEmail();
      when(samServiceMock.getUserIdFromEmail(testUser, userEmail))
          .thenReturn(TEST_SAM_USER.getSubjectId());

      // this is to simulate creation of new row in database for a user
      // that doesn't have a quota yet
      UserQuota newlyCreatedQuota =
          new UserQuota(
              PipelinesEnum.ARRAY_IMPUTATION,
              TEST_SAM_USER.getSubjectId(),
              2500, // default quota limit
              0); // zero consumed since it's a new user
      when(quotasServiceMock.getOrCreateQuotaForUserAndPipeline(
              TEST_SAM_USER.getSubjectId(), PipelinesEnum.ARRAY_IMPUTATION))
          .thenReturn(newlyCreatedQuota);

      MvcResult result =
          mockMvc
              .perform(
                  get(
                      String.format(
                          "/api/admin/v2/quotas/%s/%s",
                          PipelinesEnum.ARRAY_IMPUTATION.getLowerCaseValue(), userEmail)))
              .andExpect(status().isOk())
              .andExpect(content().contentType(MediaType.APPLICATION_JSON))
              .andReturn();

      ApiAdminQuotaV2 response =
          new ObjectMapper()
              .readValue(result.getResponse().getContentAsString(), ApiAdminQuotaV2.class);

      assertEquals(PipelinesEnum.ARRAY_IMPUTATION.getLowerCaseValue(), response.getPipelineName());
      assertEquals(TEST_SAM_USER.getSubjectId(), response.getUserId());
      assertEquals(userEmail, response.getUserEmail());
      assertEquals(2500, response.getQuotaLimit()); // default quota
      assertEquals(0, response.getQuotaConsumed()); // new user has consumed nothing
    }
  }

  @Nested
  @DisplayName("updateQuotaLimitForPipelineAndUser V2 tests")
  class UpdateAdminQuotaV2Tests {

    @Test
    void updateAdminQuotaOk() throws Exception {
      String userEmail = TEST_SAM_USER.getEmail();
      UserQuota updatedUserQuota =
          new UserQuota(PipelinesEnum.ARRAY_IMPUTATION, TEST_SAM_USER.getSubjectId(), 8000, 10);

      when(samServiceMock.getUserIdFromEmail(testUser, userEmail))
          .thenReturn(TEST_SAM_USER.getSubjectId());
      when(quotasServiceMock.getQuotaForUserAndPipeline(
              TEST_SAM_USER.getSubjectId(), PipelinesEnum.ARRAY_IMPUTATION))
          .thenReturn(Optional.of(TEST_USER_QUOTA_1));
      when(quotasServiceMock.adminUpdateQuotaLimit(TEST_USER_QUOTA_1, 8000))
          .thenReturn(updatedUserQuota);
      when(pipelinesServiceMock.getLatestPipeline(PipelinesEnum.ARRAY_IMPUTATION))
          .thenReturn(MockMvcUtils.getTestPipeline());

      MvcResult result =
          mockMvc
              .perform(
                  patch(
                          String.format(
                              "/api/admin/v2/quotas/%s/%s",
                              PipelinesEnum.ARRAY_IMPUTATION.getLowerCaseValue(),
                              TEST_SAM_USER.getEmail()))
                      .contentType(MediaType.APPLICATION_JSON)
                      .content(createTestJobPostBody(8000)))
              .andExpect(status().isOk())
              .andExpect(content().contentType(MediaType.APPLICATION_JSON))
              .andReturn();

      ApiAdminQuotaV2 response =
          new ObjectMapper()
              .readValue(result.getResponse().getContentAsString(), ApiAdminQuotaV2.class);

      assertEquals(PipelinesEnum.ARRAY_IMPUTATION.getLowerCaseValue(), response.getPipelineName());
      assertEquals(TEST_SAM_USER.getSubjectId(), response.getUserId());
      assertEquals(TEST_SAM_USER.getEmail(), response.getUserEmail());
      assertEquals(8000, response.getQuotaLimit());
      assertEquals(10, response.getQuotaConsumed());

      // assert that NotificationService was called with right parameters
      verify(notificationService)
          .configureAndSendUserQuotaChangedNotification(
              TEST_SAM_USER.getSubjectId(), "displayName", 1000, 8000, 10, 7990);
    }

    @Test
    void updateAdminQuotaRequireQuotaLimit() throws Exception {
      mockMvc
          .perform(
              patch(
                      String.format(
                          "/api/admin/v2/quotas/%s/%s",
                          PipelinesEnum.ARRAY_IMPUTATION.getLowerCaseValue(),
                          TEST_SAM_USER.getEmail()))
                  .contentType(MediaType.APPLICATION_JSON)
                  .content("{}"))
          .andExpect(status().isBadRequest());

      // assert NotificationService is never called
      verify(notificationService, never())
          .configureAndSendUserQuotaChangedNotification(
              any(), any(), anyInt(), anyInt(), anyInt(), anyInt());
    }

    @Test
    void updateAdminQuotaNotAdminUser() throws Exception {
      doThrow(new ForbiddenException("Not an admin error"))
          .when(samServiceMock)
          .checkAdminAuthz(testUser);

      mockMvc
          .perform(
              patch(
                      String.format(
                          "/api/admin/v2/quotas/%s/%s",
                          PipelinesEnum.ARRAY_IMPUTATION.getLowerCaseValue(),
                          TEST_SAM_USER.getEmail()))
                  .contentType(MediaType.APPLICATION_JSON)
                  .content(createTestJobPostBody(500)))
          .andExpect(status().isForbidden());

      // assert NotificationService is never called
      verify(notificationService, never())
          .configureAndSendUserQuotaChangedNotification(
              any(), any(), anyInt(), anyInt(), anyInt(), anyInt());
    }

    @Test
    void updateAdminQuotaUserNotFoundInSam() throws Exception {
      String userEmail = "nonexistent@example.com";
      when(samServiceMock.getUserIdFromEmail(testUser, userEmail))
          .thenThrow(new NotFoundException("User not found in SAM"));

      mockMvc
          .perform(
              patch(
                      String.format(
                          "/api/admin/v2/quotas/%s/%s",
                          PipelinesEnum.ARRAY_IMPUTATION.getLowerCaseValue(), userEmail))
                  .contentType(MediaType.APPLICATION_JSON)
                  .content(createTestJobPostBody(500)))
          .andExpect(status().isNotFound());

      // assert NotificationService is never called
      verify(notificationService, never())
          .configureAndSendUserQuotaChangedNotification(
              any(), any(), anyInt(), anyInt(), anyInt(), anyInt());
    }

    @Test
    void updateAdminQuotaUserQuotaDoesntExist() throws Exception {
      String userEmail = TEST_SAM_USER.getEmail();
      when(samServiceMock.getUserIdFromEmail(testUser, userEmail))
          .thenReturn(TEST_SAM_USER.getSubjectId());

      when(quotasServiceMock.getQuotaForUserAndPipeline(
              TEST_SAM_USER.getSubjectId(), PipelinesEnum.ARRAY_IMPUTATION))
          .thenReturn(Optional.empty());

      UserQuota newlyCreatedQuota =
          new UserQuota(PipelinesEnum.ARRAY_IMPUTATION, TEST_SAM_USER.getSubjectId(), 5000, 0);

      when(quotasServiceMock.createQuotaForUserAndPipeline(
              TEST_SAM_USER.getSubjectId(), PipelinesEnum.ARRAY_IMPUTATION, 5000))
          .thenReturn(newlyCreatedQuota);

      when(pipelinesServiceMock.getLatestPipeline(PipelinesEnum.ARRAY_IMPUTATION))
          .thenReturn(MockMvcUtils.getTestPipeline());

      MvcResult result =
          mockMvc
              .perform(
                  patch(
                          String.format(
                              "/api/admin/v2/quotas/%s/%s",
                              PipelinesEnum.ARRAY_IMPUTATION.getLowerCaseValue(),
                              TEST_SAM_USER.getEmail()))
                      .contentType(MediaType.APPLICATION_JSON)
                      .content(createTestJobPostBody(5000)))
              .andExpect(status().isOk())
              .andExpect(content().contentType(MediaType.APPLICATION_JSON))
              .andReturn();

      ApiAdminQuotaV2 response =
          new ObjectMapper()
              .readValue(result.getResponse().getContentAsString(), ApiAdminQuotaV2.class);

      assertEquals(PipelinesEnum.ARRAY_IMPUTATION.getLowerCaseValue(), response.getPipelineName());
      assertEquals(TEST_SAM_USER.getSubjectId(), response.getUserId());
      assertEquals(TEST_SAM_USER.getEmail(), response.getUserEmail());
      assertEquals(5000, response.getQuotaLimit());
      assertEquals(0, response.getQuotaConsumed());

      // assert update was not called
      verify(quotasServiceMock, never()).adminUpdateQuotaLimit(any(), anyInt());

      // assert that NotificationService was called with right parameters
      verify(notificationService)
          .configureAndSendUserQuotaChangedNotification(
              TEST_SAM_USER.getSubjectId(), "displayName", 0, 5000, 0, 5000);
    }
  }

  private String createTestJobPostBody(
      String workspaceBillingProject, String workspaceName, String toolVersion, Boolean isHidden)
      throws JsonProcessingException {
    ApiUpdatePipelineRequestBody apiUpdatePipelineRequestBody =
        new ApiUpdatePipelineRequestBody()
            .workspaceBillingProject(workspaceBillingProject)
            .workspaceName(workspaceName)
            .toolVersion(toolVersion)
            .isHidden(isHidden);
    return MockMvcUtils.convertToJsonString(apiUpdatePipelineRequestBody);
  }

  private String createTestJobPostBody(int quotaLimit) throws JsonProcessingException {
    ApiUpdateQuotaLimitRequestBody apiUpdateQuotaLimitRequestBody =
        new ApiUpdateQuotaLimitRequestBody().quotaLimit(quotaLimit);
    return MockMvcUtils.convertToJsonString(apiUpdateQuotaLimitRequestBody);
  }
}
