package bio.terra.pipelines.service;

import static org.junit.jupiter.api.Assertions.*;
import static org.mockito.ArgumentMatchers.any;
import static org.mockito.ArgumentMatchers.anyBoolean;
import static org.mockito.Mockito.when;

import bio.terra.cbas.model.MethodDetails;
import bio.terra.cbas.model.MethodListResponse;
import bio.terra.pipelines.app.configuration.internal.ImputationConfiguration;
import bio.terra.pipelines.dependencies.cbas.CbasService;
import bio.terra.pipelines.dependencies.cbas.CbasServiceApiException;
import bio.terra.pipelines.dependencies.leonardo.LeonardoService;
import bio.terra.pipelines.dependencies.leonardo.LeonardoServiceApiException;
import bio.terra.pipelines.dependencies.sam.SamService;
import bio.terra.pipelines.dependencies.stairway.JobBuilder;
import bio.terra.pipelines.dependencies.stairway.JobService;
import bio.terra.pipelines.dependencies.wds.WdsService;
import bio.terra.pipelines.dependencies.wds.WdsServiceApiException;
import bio.terra.pipelines.dependencies.wds.WdsServiceException;
import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import bio.terra.pipelines.testutils.TestUtils;
import java.util.List;
import java.util.UUID;
import org.broadinstitute.dsde.workbench.client.leonardo.ApiException;
import org.broadinstitute.dsde.workbench.client.leonardo.model.ListAppResponse;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.mockito.InjectMocks;
import org.mockito.Mock;

class ImputationServiceMockTest extends BaseEmbeddedDbTest {
  @InjectMocks private ImputationService imputationService;
  @Mock private SamService samService;
  @Mock private LeonardoService leonardoService;
  @Mock private WdsService wdsService;
  @Mock private CbasService cbasService;
  @Mock private JobService mockJobService;
  @Mock private JobBuilder mockJobBuilder;

  private final String workspaceId = "workspaceId";

  // parameters used repeatedly by various tests, and things we'll want mocks to respond to
  // universally
  private final String testUserId = TestUtils.TEST_USER_ID_1;

  private final Object testPipelineInputs = TestUtils.TEST_PIPELINE_INPUTS;
  private final UUID testUUID = TestUtils.TEST_NEW_UUID;
  @Mock ImputationConfiguration imputationConfiguration = new ImputationConfiguration(workspaceId);

  @BeforeEach
  void initMocks() {
    // stairway submit method returns a good flightId
    when(mockJobService.newJob()).thenReturn(mockJobBuilder);
    when(mockJobBuilder.jobId(any())).thenReturn(mockJobBuilder);
    when(mockJobBuilder.flightClass(any())).thenReturn(mockJobBuilder);
    when(mockJobBuilder.addParameter(any(), any())).thenReturn(mockJobBuilder);
    when(mockJobBuilder.submit()).thenReturn(testUUID);
  }

  @Test
  void queryForWorkspaceApps() {
    when(samService.getTspsServiceAccountToken()).thenReturn("saToken");
    when(leonardoService.getApps(any(), any(), anyBoolean()))
        .thenReturn(
            List.of(
                new ListAppResponse().workspaceId(workspaceId),
                new ListAppResponse().workspaceId(workspaceId)));
    when(cbasService.getAllMethods(any(), any()))
        .thenReturn(
            new MethodListResponse().addMethodsItem(new MethodDetails().name("testMethod")));

    List<ListAppResponse> response =
        imputationService.queryForWorkspaceApps(TestUtils.TEST_PIPELINE_1);
    assertEquals(2, response.size());
    response.stream().forEach(r -> assertEquals(workspaceId, r.getWorkspaceId()));
  }

  @Test
  void queryForWorkspaceAppsIsEmptyWhenLeonardoException() {
    when(samService.getTspsServiceAccountToken()).thenReturn("saToken");
    when(leonardoService.getApps(any(), any(), anyBoolean()))
        .thenThrow(new LeonardoServiceApiException(new ApiException()));
    assertTrue(imputationService.queryForWorkspaceApps(TestUtils.TEST_PIPELINE_1).isEmpty());
  }

  @Test
  void queryForWorkspaceAppsIsEmptyWhenWdsException() throws WdsServiceException {
    when(samService.getTspsServiceAccountToken()).thenReturn("saToken");
    when(wdsService.querySchema(any(), any(), any()))
        .thenThrow(
            new WdsServiceApiException(new org.databiosphere.workspacedata.client.ApiException()));
    assertTrue(imputationService.queryForWorkspaceApps(TestUtils.TEST_PIPELINE_1).isEmpty());
  }

  @Test
  void queryForWorkspaceAppsIsEmptyWhenCbasException() throws WdsServiceException {
    when(samService.getTspsServiceAccountToken()).thenReturn("saToken");
    when(cbasService.getAllMethods(any(), any()))
        .thenThrow(new CbasServiceApiException(new bio.terra.cbas.client.ApiException()));
    assertTrue(imputationService.queryForWorkspaceApps(TestUtils.TEST_PIPELINE_1).isEmpty());
  }

  @Test
  void createJobSuccess() {
    // note this doesn't actually kick off a job
    UUID writtenUUID =
        imputationService.createImputationJob(
            testUUID,
            testUserId,
            "test description",
            TestUtils.TEST_PIPELINE_1,
            testPipelineInputs,
            TestUtils.TEST_RESULT_URL);
    assertEquals(testUUID, writtenUUID);
  }
}
