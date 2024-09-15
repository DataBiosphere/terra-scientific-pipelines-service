package bio.terra.pipelines.dependencies.sam;

import static org.junit.jupiter.api.Assertions.*;
import static org.mockito.Mockito.*;

import bio.terra.common.exception.ForbiddenException;
import bio.terra.common.exception.InternalServerErrorException;
import bio.terra.common.iam.BearerToken;
import bio.terra.common.iam.SamUser;
import bio.terra.pipelines.dependencies.common.HealthCheck;
import bio.terra.pipelines.generated.model.ApiSystemStatusSystems;
import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import com.google.auth.oauth2.AccessToken;
import com.google.auth.oauth2.GoogleCredentials;
import java.io.IOException;
import org.broadinstitute.dsde.workbench.client.sam.ApiException;
import org.broadinstitute.dsde.workbench.client.sam.api.AdminApi;
import org.broadinstitute.dsde.workbench.client.sam.api.StatusApi;
import org.broadinstitute.dsde.workbench.client.sam.model.SystemStatus;
import org.broadinstitute.dsde.workbench.client.sam.model.UserStatus;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.extension.ExtendWith;
import org.mockito.InjectMocks;
import org.mockito.MockedStatic;
import org.mockito.Mockito;
import org.mockito.junit.jupiter.MockitoExtension;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.boot.test.mock.mockito.MockBean;

@ExtendWith(MockitoExtension.class)
class SamServiceTest extends BaseEmbeddedDbTest {
  @Autowired @InjectMocks SamService samService;
  @MockBean SamClient samClient;

  @Test
  void checkHealth() throws ApiException {
    StatusApi statusApi = mock(StatusApi.class);

    SystemStatus expectedSystemStatus = new SystemStatus();
    expectedSystemStatus.setOk(true);

    when(samClient.statusApi()).thenReturn(statusApi);
    when(statusApi.getSystemStatus()).thenReturn(expectedSystemStatus);

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
    StatusApi statusApi = mock(StatusApi.class);

    SystemStatus systemStatus = new SystemStatus();
    systemStatus.setOk(true);

    String exceptionMessage = "this is my exception message";
    ApiException apiException = new ApiException(exceptionMessage);
    when(samClient.statusApi()).thenReturn(statusApi);
    when(statusApi.getSystemStatus()).thenThrow(apiException);

    HealthCheck.Result expectedResultOnFail =
        new HealthCheck.Result(false, apiException.getMessage());

    HealthCheck.Result actualResult = samService.checkHealth();

    assertEquals(expectedResultOnFail, actualResult);

    ApiSystemStatusSystems apiSystemStatusSystems = samService.checkHealthApiSystemStatus();
    assertEquals(
        new ApiSystemStatusSystems()
            .ok(expectedResultOnFail.isOk())
            .addMessagesItem(expectedResultOnFail.message()),
        apiSystemStatusSystems);
  }

  @Test
  void getServiceAccountToken() throws IOException {
    GoogleCredentials mockCredentials = GoogleCredentials.create(new AccessToken("hi", null));

    try (MockedStatic<GoogleCredentials> utilities = Mockito.mockStatic(GoogleCredentials.class)) {
      SamClient samClient = mock(SamClient.class);
      utilities.when(GoogleCredentials::getApplicationDefault).thenReturn(mockCredentials);
      SamService samService = new SamService(samClient);
      String token = samService.getTeaspoonsServiceAccountToken();
      assertEquals("hi", token);
    }
  }

  /**
   * NOTE: if this test fails when you run it locally, it's likely because you have credentials
   * available locally; try running `rm $HOME/.config/gcloud/application_default_credentials.json`
   * to clear them.
   */
  @Test
  void getServiceAccountTokenThrows() {
    // this should throw an exception because there are no credentials available by default
    SamService samService = new SamService(samClient);
    assertThrows(InternalServerErrorException.class, samService::getTeaspoonsServiceAccountToken);
  }

  @Test
  void isAdmin() throws ApiException {
    AdminApi adminApi = mock(AdminApi.class);

    when(adminApi.adminGetUserByEmail("blahEmail")).thenReturn(new UserStatus());
    when(samClient.adminApi("blahToken")).thenReturn(adminApi);

    samService.checkAdminAuthz(
        new SamUser("blahEmail", "blahSubjectId", new BearerToken("blahToken")));
  }

  @Test
  void isAdminNotAdminForbiddenException() throws ApiException {
    SamClient samClient = mock(SamClient.class);
    AdminApi adminApi = mock(AdminApi.class);

    when(adminApi.adminGetUserByEmail("doesnt matter")).thenThrow(ApiException.class);
    when(samClient.adminApi("blah")).thenReturn(adminApi);

    SamService samService = spy(new SamService(samClient));

    SamUser samUser = new SamUser("doesnt matter", "really doesnt matter", new BearerToken("blah"));
    assertThrows(
        ForbiddenException.class,
        () -> {
          samService.checkAdminAuthz(samUser);
        });
  }
}
