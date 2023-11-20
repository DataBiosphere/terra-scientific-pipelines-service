package bio.terra.pipelines.service;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.mockito.ArgumentMatchers.any;
import static org.mockito.ArgumentMatchers.anyBoolean;
import static org.mockito.Mockito.when;

import bio.terra.pipelines.app.configuration.internal.ImputationConfiguration;
import bio.terra.pipelines.dependencies.leonardo.LeonardoService;
import bio.terra.pipelines.dependencies.leonardo.LeonardoServiceApiException;
import bio.terra.pipelines.dependencies.sam.SamService;
import bio.terra.pipelines.testutils.BaseContainerTest;
import java.util.List;
import org.broadinstitute.dsde.workbench.client.leonardo.ApiException;
import org.broadinstitute.dsde.workbench.client.leonardo.model.ListAppResponse;
import org.junit.jupiter.api.Test;
import org.mockito.InjectMocks;
import org.mockito.Mock;

public class ImputationServiceMockTest extends BaseContainerTest {
  @InjectMocks private ImputationService imputationService;
  @Mock private SamService samService;
  @Mock private LeonardoService leonardoService;

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
  void queryForWorkspaceAppsIsNullWhenLeonardoException() {
    when(samService.getTspsServiceAccountToken()).thenReturn("saToken");
    when(leonardoService.getApps(any(), any(), anyBoolean()))
        .thenThrow(new LeonardoServiceApiException(new ApiException()));
    assertNull(imputationService.queryForWorkspaceApps());
  }
}
