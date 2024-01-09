package bio.terra.pipelines.service;

import static org.junit.jupiter.api.Assertions.*;
import static org.mockito.ArgumentMatchers.any;
import static org.mockito.ArgumentMatchers.anyBoolean;
import static org.mockito.Mockito.when;

import bio.terra.pipelines.app.configuration.internal.ImputationConfiguration;
import bio.terra.pipelines.dependencies.leonardo.LeonardoService;
import bio.terra.pipelines.dependencies.leonardo.LeonardoServiceApiException;
import bio.terra.pipelines.dependencies.sam.SamService;
import bio.terra.pipelines.dependencies.stairway.StairwayJobBuilder;
import bio.terra.pipelines.dependencies.stairway.StairwayJobService;
import bio.terra.pipelines.dependencies.wds.WdsService;
import bio.terra.pipelines.dependencies.wds.WdsServiceApiException;
import bio.terra.pipelines.dependencies.wds.WdsServiceException;
import bio.terra.pipelines.testutils.BaseContainerTest;
import bio.terra.pipelines.testutils.TestUtils;
import java.util.List;
import java.util.UUID;
import org.broadinstitute.dsde.workbench.client.leonardo.ApiException;
import org.broadinstitute.dsde.workbench.client.leonardo.model.ListAppResponse;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.mockito.InjectMocks;
import org.mockito.Mock;

class ImputationServiceMockTest extends BaseContainerTest {
  @InjectMocks private ImputationService imputationService;
  @Mock private SamService samService;
  @Mock private LeonardoService leonardoService;
  @Mock private WdsService wdsService;
  @Mock private StairwayJobService mockStairwayJobService;
  @Mock private StairwayJobBuilder mockStairwayJobBuilder;

  private final String workspaceId = "workspaceId";

  // parameters used repeatedly by various tests, and things we'll want mocks to respond to
  // universally
  private final String testUserId = TestUtils.TEST_USER_ID_1;
  private final String testPipelineVersion = TestUtils.TEST_PIPELINE_VERSION_1;

  private final Object testPipelineInputs = TestUtils.TEST_PIPELINE_INPUTS;
  private final UUID testUUID = TestUtils.TEST_NEW_UUID;
  @Mock ImputationConfiguration imputationConfiguration = new ImputationConfiguration(workspaceId);

  @BeforeEach
  void initMocks() {
    // stairway submit method returns a good flightId
    when(mockStairwayJobService.newJob()).thenReturn(mockStairwayJobBuilder);
    when(mockStairwayJobBuilder.jobId(any())).thenReturn(mockStairwayJobBuilder);
    when(mockStairwayJobBuilder.flightClass(any())).thenReturn(mockStairwayJobBuilder);
    when(mockStairwayJobBuilder.pipelineId(any())).thenReturn(mockStairwayJobBuilder);
    when(mockStairwayJobBuilder.pipelineVersion(any())).thenReturn(mockStairwayJobBuilder);
    when(mockStairwayJobBuilder.userId(any())).thenReturn(mockStairwayJobBuilder);
    when(mockStairwayJobBuilder.pipelineInputs(any())).thenReturn(mockStairwayJobBuilder);
    when(mockStairwayJobBuilder.submit()).thenReturn(testUUID);
  }

  @Test
  void queryForWorkspaceApps() {
    when(samService.getTspsServiceAccountToken()).thenReturn("saToken");
    when(leonardoService.getApps(any(), any(), anyBoolean()))
        .thenReturn(
            List.of(
                new ListAppResponse().workspaceId(workspaceId),
                new ListAppResponse().workspaceId(workspaceId)));

    List<ListAppResponse> response = imputationService.queryForWorkspaceApps();
    assertEquals(2, response.size());
    response.stream().forEach(r -> assertEquals(workspaceId, r.getWorkspaceId()));
  }

  @Test
  void queryForWorkspaceAppsIsEmptyWhenLeonardoException() {
    when(samService.getTspsServiceAccountToken()).thenReturn("saToken");
    when(leonardoService.getApps(any(), any(), anyBoolean()))
        .thenThrow(new LeonardoServiceApiException(new ApiException()));
    assertTrue(imputationService.queryForWorkspaceApps().isEmpty());
  }

  @Test
  void queryForWorkspaceAppsIsEmptyWhenWdsException() throws WdsServiceException {
    when(samService.getTspsServiceAccountToken()).thenReturn("saToken");
    when(wdsService.querySchema(any(), any(), any()))
        .thenThrow(
            new WdsServiceApiException(new org.databiosphere.workspacedata.client.ApiException()));
    assertTrue(imputationService.queryForWorkspaceApps().isEmpty());
  }

  @Test
  void testCreateJob_success() {
    // a job isn't actually kicked off
    UUID writtenUUID =
        imputationService.createImputationJob(testUserId, testPipelineVersion, testPipelineInputs);
    assertEquals(testUUID, writtenUUID);
  }
}
