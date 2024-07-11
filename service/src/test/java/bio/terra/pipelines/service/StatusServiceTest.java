package bio.terra.pipelines.service;

import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.mockito.ArgumentMatchers.any;
import static org.mockito.Mockito.doReturn;
import static org.mockito.Mockito.doThrow;
import static org.mockito.Mockito.when;

import bio.terra.pipelines.app.configuration.internal.StatusCheckConfiguration;
import bio.terra.pipelines.dependencies.sam.SamService;
import bio.terra.pipelines.dependencies.workspacemanager.WorkspaceManagerService;
import bio.terra.pipelines.generated.model.ApiSystemStatusSystems;
import bio.terra.pipelines.testutils.BaseTest;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.extension.ExtendWith;
import org.mockito.Mock;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.boot.test.mock.mockito.MockBean;
import org.springframework.jdbc.core.ConnectionCallback;
import org.springframework.jdbc.core.JdbcTemplate;
import org.springframework.jdbc.core.namedparam.NamedParameterJdbcTemplate;
import org.springframework.test.context.ContextConfiguration;
import org.springframework.test.context.junit.jupiter.SpringExtension;

@ExtendWith(SpringExtension.class)
@ContextConfiguration(classes = StatusService.class)
class StatusServiceTest extends BaseTest {
  @Autowired private StatusService statusService;

  @MockBean private NamedParameterJdbcTemplate jdbcTemplate;
  @MockBean private StatusCheckConfiguration configuration;
  @MockBean private SamService samService;
  @MockBean private WorkspaceManagerService workspaceManagerService;

  @Mock private final JdbcTemplate jdbcTemplateMock = new JdbcTemplate();

  @BeforeEach
  void setup() {
    when(jdbcTemplate.getJdbcTemplate()).thenReturn(jdbcTemplateMock);
    when(configuration.enabled()).thenReturn(true);
    when(configuration.stalenessThresholdSeconds()).thenReturn(60);
  }

  @Test
  void testStatusOk() {
    doReturn(true).when(jdbcTemplateMock).execute(any(ConnectionCallback.class));
    when(samService.checkHealthApiSystemStatus()).thenReturn(new ApiSystemStatusSystems().ok(true));
    when(workspaceManagerService.checkHealthApiSystemStatus())
        .thenReturn(new ApiSystemStatusSystems().ok(true));

    statusService.checkStatus();
    assertTrue(statusService.getCurrentStatus());
  }

  @Test
  void testStatusStale() {
    // make the staleness threshold 0 seconds
    when(configuration.stalenessThresholdSeconds()).thenReturn(0);

    doReturn(true).when(jdbcTemplateMock).execute(any(ConnectionCallback.class));
    when(samService.checkHealthApiSystemStatus()).thenReturn(new ApiSystemStatusSystems().ok(true));
    when(workspaceManagerService.checkHealthApiSystemStatus())
        .thenReturn(new ApiSystemStatusSystems().ok(true));

    statusService.checkStatus();
    assertFalse(statusService.getCurrentStatus());
  }

  @Test
  void testStatusBadDb() {
    doThrow(new RuntimeException()).when(jdbcTemplateMock).execute(any(ConnectionCallback.class));
    when(samService.checkHealthApiSystemStatus()).thenReturn(new ApiSystemStatusSystems().ok(true));
    when(workspaceManagerService.checkHealthApiSystemStatus())
        .thenReturn(new ApiSystemStatusSystems().ok(true));

    statusService.checkStatus();
    assertFalse(statusService.getCurrentStatus());
  }

  @Test
  void testStatusBadWsm() {
    doReturn(true).when(jdbcTemplateMock).execute(any(ConnectionCallback.class));
    when(samService.checkHealthApiSystemStatus()).thenReturn(new ApiSystemStatusSystems().ok(true));
    when(workspaceManagerService.checkHealthApiSystemStatus())
        .thenReturn(new ApiSystemStatusSystems().ok(false));

    statusService.checkStatus();
    assertFalse(statusService.getCurrentStatus());
  }

  @Test
  void testStatusBadSam() {
    doReturn(true).when(jdbcTemplateMock).execute(any(ConnectionCallback.class));
    when(samService.checkHealthApiSystemStatus())
        .thenReturn(new ApiSystemStatusSystems().ok(false));
    when(workspaceManagerService.checkHealthApiSystemStatus())
        .thenReturn(new ApiSystemStatusSystems().ok(true));

    statusService.checkStatus();
    assertFalse(statusService.getCurrentStatus());
  }
}
