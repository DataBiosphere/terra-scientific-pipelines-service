package bio.terra.pipelines.dependencies.sam;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.mockito.Mockito.*;

import bio.terra.pipelines.dependencies.common.HealthCheck;
import bio.terra.pipelines.generated.model.ApiSystemStatusSystems;
import org.broadinstitute.dsde.workbench.client.sam.ApiException;
import org.broadinstitute.dsde.workbench.client.sam.api.StatusApi;
import org.broadinstitute.dsde.workbench.client.sam.model.SystemStatus;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.extension.ExtendWith;
import org.mockito.junit.jupiter.MockitoExtension;

@ExtendWith(MockitoExtension.class)
class SamServiceTest {

  @Test
  void checkHealth() throws ApiException {
    SamClient samClient = mock(SamClient.class);
    StatusApi statusApi = mock(StatusApi.class);

    SystemStatus expectedSystemStatus = new SystemStatus();
    expectedSystemStatus.setOk(true);

    when(samClient.statusApi()).thenReturn(statusApi);
    when(statusApi.getSystemStatus()).thenReturn(expectedSystemStatus);

    SamService samService = spy(new SamService(samClient));

    HealthCheck.Result actualResult = samService.checkHealth();

    assertEquals(
        new HealthCheck.Result(expectedSystemStatus.getOk(), expectedSystemStatus.toString()),
        actualResult);

    ApiSystemStatusSystems apiSystemStatusSystems = samService.checkHealthApiSystemStatus();
    assertEquals(
        new ApiSystemStatusSystems()
            .ok(expectedSystemStatus.getOk())
            .addMessagesItem(expectedSystemStatus.toString()),
        apiSystemStatusSystems);
  }

  @Test
  void checkHealthWithException() throws ApiException {
    SamClient samClient = mock(SamClient.class);
    StatusApi statusApi = mock(StatusApi.class);

    SystemStatus systemStatus = new SystemStatus();
    systemStatus.setOk(true);

    String exceptionMessage = "this is my exception message";
    ApiException apiException = new ApiException(exceptionMessage);
    when(samClient.statusApi()).thenReturn(statusApi);
    when(statusApi.getSystemStatus()).thenThrow(apiException);

    HealthCheck.Result expectedResultOnFail =
        new HealthCheck.Result(false, apiException.getMessage());
    SamService samService = spy(new SamService(samClient));

    HealthCheck.Result actualResult = samService.checkHealth();

    assertEquals(expectedResultOnFail, actualResult);

    ApiSystemStatusSystems apiSystemStatusSystems = samService.checkHealthApiSystemStatus();
    assertEquals(
        new ApiSystemStatusSystems()
            .ok(expectedResultOnFail.isOk())
            .addMessagesItem(expectedResultOnFail.message()),
        apiSystemStatusSystems);
  }
}
