package bio.terra.pipelines.service;

import static org.junit.jupiter.api.Assertions.*;
import static org.mockito.ArgumentMatchers.any;
import static org.mockito.ArgumentMatchers.anyBoolean;
import static org.mockito.Mockito.when;

import bio.terra.pipelines.app.configuration.internal.ImputationConfiguration;
import bio.terra.pipelines.dependencies.leonardo.LeonardoService;
import bio.terra.pipelines.dependencies.leonardo.LeonardoServiceApiException;
import bio.terra.pipelines.dependencies.sam.SamService;
import bio.terra.pipelines.dependencies.wds.WdsService;
import bio.terra.pipelines.dependencies.wds.WdsServiceApiException;
import bio.terra.pipelines.dependencies.wds.WdsServiceException;
import bio.terra.pipelines.testutils.BaseContainerTest;
import java.util.List;
import org.broadinstitute.dsde.workbench.client.leonardo.ApiException;
import org.broadinstitute.dsde.workbench.client.leonardo.model.ListAppResponse;
import org.junit.jupiter.api.Test;
import org.mockito.InjectMocks;
import org.mockito.Mock;

class ImputationServiceMockTest extends BaseContainerTest {
  @InjectMocks private ImputationService imputationService;
  @Mock private SamService samService;
  @Mock private LeonardoService leonardoService;
  @Mock private WdsService wdsService;

  private final String workspaceId = "workspaceId";
  @Mock ImputationConfiguration imputationConfiguration = new ImputationConfiguration(workspaceId);

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
}
