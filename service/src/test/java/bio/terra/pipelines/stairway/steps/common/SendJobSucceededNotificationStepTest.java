package bio.terra.pipelines.stairway.steps.common;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.mockito.Mockito.doNothing;
import static org.mockito.Mockito.doThrow;
import static org.mockito.Mockito.times;
import static org.mockito.Mockito.verify;
import static org.mockito.Mockito.when;

import bio.terra.pipelines.dependencies.stairway.JobMapKeys;
import bio.terra.pipelines.notifications.NotificationService;
import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import bio.terra.pipelines.testutils.TestUtils;
import bio.terra.stairway.FlightContext;
import bio.terra.stairway.FlightMap;
import bio.terra.stairway.StepStatus;
import java.util.UUID;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.mockito.Mock;

class SendJobSucceededNotificationStepTest extends BaseEmbeddedDbTest {
  @Mock private NotificationService notificationService;
  @Mock private FlightContext flightContext;

  private final UUID testJobId = TestUtils.TEST_NEW_UUID;
  private final String testUserId = TestUtils.TEST_USER_1_ID;

  @BeforeEach
  void setup() {
    var inputParameters = new FlightMap();
    var workingMap = new FlightMap();

    inputParameters.put(JobMapKeys.USER_ID, testUserId);

    when(flightContext.getFlightId()).thenReturn(testJobId.toString());
    when(flightContext.getInputParameters()).thenReturn(inputParameters);
    when(flightContext.getWorkingMap()).thenReturn(workingMap);
  }

  @Test
  void doStepSuccess() {
    // assume success
    doNothing()
        .when(notificationService)
        .configureAndSendPipelineRunSucceededNotification(testJobId, testUserId);

    var sendJobSucceededNotificationStep =
        new SendJobSucceededNotificationStep(notificationService);
    var result = sendJobSucceededNotificationStep.doStep(flightContext);

    assertEquals(StepStatus.STEP_RESULT_SUCCESS, result.getStepStatus());
    verify(notificationService, times(1))
        .configureAndSendPipelineRunSucceededNotification(testJobId, testUserId);
  }

  @Test
  void doStepFailureStillSuccess() {
    // throw exception
    doThrow(new RuntimeException())
        .when(notificationService)
        .configureAndSendPipelineRunSucceededNotification(testJobId, testUserId);

    var sendJobSucceededNotificationStep =
        new SendJobSucceededNotificationStep(notificationService);
    var result = sendJobSucceededNotificationStep.doStep(flightContext);

    // step should catch exception and just log the error rather than throwing
    assertEquals(StepStatus.STEP_RESULT_SUCCESS, result.getStepStatus());
    verify(notificationService, times(1))
        .configureAndSendPipelineRunSucceededNotification(testJobId, testUserId);
  }

  @Test
  void undoStepSuccess() {
    var sendJobSucceededNotificationStep =
        new SendJobSucceededNotificationStep(notificationService);
    var result = sendJobSucceededNotificationStep.undoStep(flightContext);

    assertEquals(StepStatus.STEP_RESULT_SUCCESS, result.getStepStatus());
  }
}
