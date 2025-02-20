package bio.terra.pipelines.service;

import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.mockito.ArgumentMatchers.any;
import static org.mockito.Mockito.doReturn;
import static org.mockito.Mockito.doThrow;
import static org.mockito.Mockito.when;

import bio.terra.pipelines.app.configuration.internal.StatusCheckConfiguration;
import bio.terra.pipelines.dependencies.sam.SamService;
import bio.terra.pipelines.generated.model.ApiSystemStatusSystems;
import bio.terra.pipelines.testutils.BaseTest;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.extension.ExtendWith;
import org.mockito.Mock;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.jdbc.core.ConnectionCallback;
import org.springframework.jdbc.core.JdbcTemplate;
import org.springframework.jdbc.core.namedparam.NamedParameterJdbcTemplate;
import org.springframework.test.context.ContextConfiguration;
import org.springframework.test.context.bean.override.mockito.MockitoBean;
import org.springframework.test.context.junit.jupiter.SpringExtension;

@ExtendWith(SpringExtension.class)
@ContextConfiguration(classes = StatusService.class)
class StatusServiceTest extends BaseTest {
  @Autowired private StatusService statusService;

  @MockitoBean private NamedParameterJdbcTemplate jdbcTemplate;
  @MockitoBean private StatusCheckConfiguration configuration;
  @MockitoBean private SamService samService;

  @Mock private final JdbcTemplate jdbcTemplateMock = new JdbcTemplate();

  @BeforeEach
  void setup() {
    // by default, everything should work. we can override some of these in specific tests
    when(jdbcTemplate.getJdbcTemplate()).thenReturn(jdbcTemplateMock);
    when(configuration.enabled()).thenReturn(true);
    when(configuration.stalenessThresholdSeconds()).thenReturn(60);

    doReturn(true).when(jdbcTemplateMock).execute(any(ConnectionCallback.class));
    when(samService.checkHealthApiSystemStatus()).thenReturn(new ApiSystemStatusSystems().ok(true));
  }

  @Test
  void testStatusOk() {
    statusService.checkStatus();
    assertTrue(statusService.getCurrentStatus());
  }

  @Test
  void testStatusConfigurationNotEnabled() {
    // when status checking is disabled, the status defaults to true (ok)
    when(configuration.enabled()).thenReturn(false);

    statusService.checkStatus();
    assertTrue(statusService.getCurrentStatus());
  }

  @Test
  void testStatusStale() {
    // make the staleness threshold 0 seconds
    when(configuration.stalenessThresholdSeconds()).thenReturn(0);

    statusService.checkStatus();
    assertFalse(statusService.getCurrentStatus());
  }

  @Test
  void testStatusBadDb() {
    doThrow(new RuntimeException()).when(jdbcTemplateMock).execute(any(ConnectionCallback.class));

    statusService.checkStatus();
    assertFalse(statusService.getCurrentStatus());
  }

  @Test
  void testStatusBadSam() {
    when(samService.checkHealthApiSystemStatus())
        .thenReturn(new ApiSystemStatusSystems().ok(false));

    statusService.checkStatus();
    assertFalse(statusService.getCurrentStatus());
  }
}
