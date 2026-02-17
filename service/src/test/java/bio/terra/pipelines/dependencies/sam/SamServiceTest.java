package bio.terra.pipelines.dependencies.sam;

import static org.junit.jupiter.api.Assertions.*;
import static org.mockito.Mockito.*;

import bio.terra.common.exception.ForbiddenException;
import bio.terra.common.exception.InternalServerErrorException;
import bio.terra.common.iam.BearerToken;
import bio.terra.common.iam.SamUser;
import bio.terra.common.sam.exception.SamBadRequestException;
import bio.terra.pipelines.dependencies.common.HealthCheck;
import bio.terra.pipelines.generated.model.ApiSystemStatusSystems;
import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import bio.terra.pipelines.testutils.TestUtils;
import com.google.auth.oauth2.AccessToken;
import com.google.auth.oauth2.GoogleCredentials;
import java.util.List;
import org.broadinstitute.dsde.workbench.client.sam.ApiException;
import org.broadinstitute.dsde.workbench.client.sam.api.AdminApi;
import org.broadinstitute.dsde.workbench.client.sam.api.GoogleApi;
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
import org.springframework.test.context.bean.override.mockito.MockitoBean;

@ExtendWith(MockitoExtension.class)
class SamServiceTest extends BaseEmbeddedDbTest {
  @Autowired @InjectMocks SamService samService;
  @MockitoBean SamClient samClient;

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
  void getProxyGroupForUser() throws ApiException {
    SamUser testUser = TestUtils.TEST_SAM_USER_1;
    String expectedProxyGroup = "proxyGroup";
    GoogleApi googleApi = mock(GoogleApi.class);
    when(googleApi.getProxyGroup(testUser.getEmail())).thenReturn(expectedProxyGroup);
    when(samClient.googleApi(testUser.getBearerToken().getToken())).thenReturn(googleApi);

    assertEquals(expectedProxyGroup, samService.getProxyGroupForUser(testUser));
  }

  @Test
  void getUserPetServiceAccountTokenReadOnly() throws ApiException {
    List<String> readOnlyScopes =
        List.of(
            "https://www.googleapis.com/auth/userinfo.email",
            "https://www.googleapis.com/auth/userinfo.profile",
            "https://www.googleapis.com/auth/devstorage.read_only");
    SamUser testUser = TestUtils.TEST_SAM_USER_1;
    String expectedTokenString = "petToken";
    GoogleApi googleApi = mock(GoogleApi.class);
    when(googleApi.getArbitraryPetServiceAccountToken(readOnlyScopes)).thenReturn("petToken");
    when(samClient.googleApi(testUser.getBearerToken().getToken())).thenReturn(googleApi);

    assertEquals(
        new BearerToken(expectedTokenString),
        samService.getUserPetServiceAccountTokenReadOnly(testUser));
  }

  @Test
  void getUserPetServiceAccountTokenReadOnlyApiException() throws ApiException {
    List<String> readOnlyScopes =
        List.of(
            "https://www.googleapis.com/auth/userinfo.email",
            "https://www.googleapis.com/auth/userinfo.profile",
            "https://www.googleapis.com/auth/devstorage.read_only");
    SamUser testUser = TestUtils.TEST_SAM_USER_1;
    GoogleApi googleApi = mock(GoogleApi.class);
    ApiException apiException = new ApiException(400, "this is an api exception");
    when(googleApi.getArbitraryPetServiceAccountToken(readOnlyScopes)).thenThrow(apiException);
    when(samClient.googleApi(testUser.getBearerToken().getToken())).thenReturn(googleApi);

    // SamService (via TCL) translates the 400 ApiException into a SamBadRequestException
    assertThrows(
        SamBadRequestException.class,
        () -> samService.getUserPetServiceAccountTokenReadOnly(testUser));
  }

  @Test
  void getServiceAccountToken() {
    GoogleCredentials mockCredentials = GoogleCredentials.create(new AccessToken("hi", null));

    try (MockedStatic<GoogleCredentials> utilities = Mockito.mockStatic(GoogleCredentials.class)) {
      utilities.when(GoogleCredentials::getApplicationDefault).thenReturn(mockCredentials);
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
    assertThrows(InternalServerErrorException.class, samService::getTeaspoonsServiceAccountToken);
  }

  @Test
  void isAdmin() throws ApiException {
    AdminApi adminApi = mock(AdminApi.class);

    when(adminApi.adminGetUserByEmail("blahEmail")).thenReturn(new UserStatus());
    when(samClient.adminApi("blahToken")).thenReturn(adminApi);

    samService.checkAdminAuthz(
        new SamUser("blahEmail", "doesnt matter", new BearerToken("blahToken")));
  }

  @Test
  void isAdminNotAdminForbiddenException() throws ApiException {
    AdminApi adminApi = mock(AdminApi.class);

    String userEmail = "doesnt matter";
    when(adminApi.adminGetUserByEmail(userEmail)).thenThrow(ApiException.class);
    when(samClient.adminApi("blah")).thenReturn(adminApi);

    SamUser samUser = new SamUser(userEmail, "really doesnt matter", new BearerToken("blah"));
    assertThrows(
        ForbiddenException.class,
        () -> {
          samService.checkAdminAuthz(samUser);
        });
  }
}
