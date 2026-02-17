package bio.terra.pipelines.common.utils;

import static org.junit.jupiter.params.provider.Arguments.arguments;
import static org.mockito.Mockito.doNothing;
import static org.mockito.Mockito.verify;
import static org.mockito.Mockito.verifyNoInteractions;

import bio.terra.pipelines.dependencies.stairway.JobMapKeys;
import bio.terra.pipelines.notifications.NotificationService;
import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import bio.terra.pipelines.testutils.StairwayTestUtils;
import bio.terra.pipelines.testutils.TestFlightContext;
import bio.terra.pipelines.testutils.TestUtils;
import bio.terra.stairway.FlightContext;
import bio.terra.stairway.FlightStatus;
import java.util.UUID;
import java.util.stream.Stream;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.Arguments;
import org.junit.jupiter.params.provider.MethodSource;
import org.mockito.InjectMocks;
import org.mockito.Mock;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.test.context.bean.override.mockito.MockitoBean;

class StairwaySendFailedJobNotificationHookTest extends BaseEmbeddedDbTest {
  @Autowired @InjectMocks
  StairwaySendFailedJobNotificationHook stairwaySendFailedJobNotificationHook;

  @MockitoBean NotificationService notificationService;

  @Mock private FlightContext flightContext;

  private final UUID testJobId = TestUtils.TEST_NEW_UUID;

  @BeforeEach
  void setup() {
    doNothing()
        .when(notificationService)
        .configureAndSendPipelineRunFailedNotification(
            testJobId, TestUtils.TEST_USER_1_ID, flightContext);
    stairwaySendFailedJobNotificationHook =
        new StairwaySendFailedJobNotificationHook(notificationService);
  }

  private static Stream<Arguments> flightContexts() {

    return Stream.of(
        // arguments: whether to include hook key in input Params,
        // hook key value, flight status at endFlight time, whether failed job notification should
        // be sent

        arguments(
            true,
            true,
            FlightStatus.SUCCESS,
            false), // flight was successful, so the pipelineRun status should not be set to FAILED
        arguments(
            true,
            true,
            FlightStatus.ERROR,
            true), // flight failed, so the pipelineRun status should be set to FAILED
        arguments(
            true,
            true,
            FlightStatus.FATAL,
            true), // flight failed dismally, so the pipelineRun status should be set to FAILED
        arguments(
            true,
            true,
            FlightStatus.QUEUED,
            true), // flight in an unexpected status for end of flight, so the pipelineRun status
        // should be set to FAILED
        arguments(
            false,
            false, // doesn't matter
            FlightStatus.ERROR,
            false), // flight failed, but the hook key was not included in the input parameters
        arguments(
            true,
            false,
            FlightStatus.ERROR,
            false)); // flight failed, but the hook key value is false
  }

  @ParameterizedTest
  @MethodSource("flightContexts")
  <T> void endFlight(
      boolean includeHookKey,
      boolean hookKeyValue,
      FlightStatus endFlightStatus,
      boolean shouldSendFailedJobNotification)
      throws InterruptedException {
    var context = new TestFlightContext();

    if (includeHookKey) {
      // this includes setting the DO_SET_PIPELINE_RUN_STATUS_FAILED_HOOK key to true
      StairwayTestUtils.constructCreateJobInputs(context.getInputParameters());
      context
          .getInputParameters()
          .put(JobMapKeys.DO_SEND_JOB_FAILURE_NOTIFICATION_HOOK, hookKeyValue);
    }

    stairwaySendFailedJobNotificationHook.startFlight(context);

    // set the end flight status
    context.flightStatus(endFlightStatus);

    stairwaySendFailedJobNotificationHook.endFlight(context);

    if (shouldSendFailedJobNotification) {
      verify(notificationService)
          .configureAndSendPipelineRunFailedNotification(
              testJobId, TestUtils.TEST_USER_1_ID, context);
    } else {
      verifyNoInteractions(notificationService);
    }
  }
}
