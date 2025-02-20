package bio.terra.pipelines.notifications;

import static org.junit.jupiter.api.Assertions.assertDoesNotThrow;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.mockito.ArgumentMatchers.any;
import static org.mockito.Mockito.doNothing;
import static org.mockito.Mockito.doThrow;
import static org.mockito.Mockito.times;
import static org.mockito.Mockito.verify;
import static org.mockito.Mockito.when;

import bio.terra.pipelines.app.configuration.internal.NotificationConfiguration;
import bio.terra.pipelines.app.configuration.internal.PipelinesCommonConfiguration;
import bio.terra.pipelines.common.utils.CommonPipelineRunStatusEnum;
import bio.terra.pipelines.db.entities.Pipeline;
import bio.terra.pipelines.db.entities.PipelineRun;
import bio.terra.pipelines.db.entities.UserQuota;
import bio.terra.pipelines.db.repositories.PipelineRunsRepository;
import bio.terra.pipelines.dependencies.rawls.RawlsServiceApiException;
import bio.terra.pipelines.service.PipelineRunsService;
import bio.terra.pipelines.service.PipelinesService;
import bio.terra.pipelines.service.QuotasService;
import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import bio.terra.pipelines.testutils.TestUtils;
import bio.terra.stairway.FlightContext;
import bio.terra.stairway.StepResult;
import bio.terra.stairway.StepStatus;
import com.fasterxml.jackson.databind.ObjectMapper;
import java.io.IOException;
import java.time.Instant;
import java.util.UUID;
import org.junit.jupiter.api.Test;
import org.mockito.InjectMocks;
import org.mockito.Mock;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.test.context.bean.override.mockito.MockitoBean;

class NotificationServiceTest extends BaseEmbeddedDbTest {
  @InjectMocks @Autowired NotificationService notificationService;
  @Autowired PipelineRunsService pipelineRunsService;
  @Autowired PipelineRunsRepository pipelineRunsRepository;
  @Autowired PipelinesService pipelinesService;
  @Autowired QuotasService quotasService;
  @Autowired NotificationConfiguration notificationConfiguration;
  @Autowired PipelinesCommonConfiguration pipelinesCommonConfiguration;
  @Autowired ObjectMapper objectMapper;
  @MockitoBean PubsubService pubsubService;
  @Mock private FlightContext flightContext;

  UUID testJobId = TestUtils.TEST_NEW_UUID;
  String testUserId = TestUtils.TEST_USER_ID_1;
  Integer testQuotaConsumedByJob = 1000;
  String testUserDescription = TestUtils.TEST_USER_PROVIDED_DESCRIPTION;
  String testErrorMessage = "test error message";

  @Test
  void formatDateTime() {
    Instant instant = Instant.parse("2021-08-25T12:34:56.789Z");
    String formattedDateTime = notificationService.formatInstantToReadableString(instant);
    assertEquals("Wed, 25 Aug 2021 12:34:56 GMT", formattedDateTime);
  }

  @Test
  void configureAndSendPipelineRunSucceededNotification() throws IOException {
    Pipeline pipeline = pipelinesService.getPipelineById(1L);
    PipelineRun writtenPipelineRun =
        createCompletedPipelineRunInDb(pipeline, CommonPipelineRunStatusEnum.SUCCEEDED);

    // initialize and set user quota
    UserQuota userQuota =
        quotasService.getOrCreateQuotaForUserAndPipeline(testUserId, pipeline.getName());
    UserQuota updatedUserQuota =
        quotasService.updateQuotaConsumed(userQuota, testQuotaConsumedByJob);
    int expectedQuotaRemaining = updatedUserQuota.getQuota() - updatedUserQuota.getQuotaConsumed();

    String stringifiedJobSucceededNotification =
        objectMapper.writeValueAsString(
            new TeaspoonsJobSucceededNotification(
                testUserId,
                pipeline.getDisplayName(),
                testJobId.toString(),
                notificationService.formatInstantToReadableString(writtenPipelineRun.getCreated()),
                notificationService.formatInstantToReadableString(writtenPipelineRun.getUpdated()),
                testQuotaConsumedByJob.toString(),
                String.valueOf(expectedQuotaRemaining),
                testUserDescription,
                pipelinesCommonConfiguration.getUserDataTtlDays().toString()));
    // success is a void method
    doNothing()
        .when(pubsubService)
        .publishMessage(
            notificationConfiguration.projectId(),
            notificationConfiguration.topicId(),
            stringifiedJobSucceededNotification);

    notificationService.configureAndSendPipelineRunSucceededNotification(testJobId, testUserId);

    // verify that the pubsub method was called
    verify(pubsubService, times(1))
        .publishMessage(
            notificationConfiguration.projectId(),
            notificationConfiguration.topicId(),
            stringifiedJobSucceededNotification);
  }

  @Test
  void configureAndSendPipelineRunSucceededNotificationIOException() throws IOException {
    Pipeline pipeline = pipelinesService.getPipelineById(1L);
    createCompletedPipelineRunInDb(pipeline, CommonPipelineRunStatusEnum.SUCCEEDED);

    doThrow(new IOException()).when(pubsubService).publishMessage(any(), any(), any());

    // exception should be caught
    assertDoesNotThrow(
        () ->
            notificationService.configureAndSendPipelineRunSucceededNotification(
                testJobId, testUserId));
  }

  @Test
  void configureAndSendPipelineRunFailedNotification() throws IOException {
    Pipeline pipeline = pipelinesService.getPipelineById(1L);
    PipelineRun writtenPipelineRun =
        createCompletedPipelineRunInDb(pipeline, CommonPipelineRunStatusEnum.FAILED);

    // initialize and set a custom user quota value. this is not quota consumed by the job.
    int customUserQuota = 2000;
    UserQuota userQuota =
        quotasService.getOrCreateQuotaForUserAndPipeline(testUserId, pipeline.getName());
    UserQuota updatedUserQuota = quotasService.updateQuotaConsumed(userQuota, customUserQuota);
    int expectedQuotaRemaining = updatedUserQuota.getQuota() - updatedUserQuota.getQuotaConsumed();

    when(flightContext.getFlightId()).thenReturn(testJobId.toString());
    RawlsServiceApiException rawlsServiceApiException =
        new RawlsServiceApiException(testErrorMessage);
    StepResult stepResultFailedWithException =
        new StepResult(StepStatus.STEP_RESULT_FAILURE_FATAL, rawlsServiceApiException);
    when(flightContext.getResult()).thenReturn(stepResultFailedWithException);

    String stringifiedJobFailedNotification =
        objectMapper.writeValueAsString(
            new TeaspoonsJobFailedNotification(
                testUserId,
                pipeline.getDisplayName(),
                testJobId.toString(),
                testErrorMessage,
                notificationService.formatInstantToReadableString(writtenPipelineRun.getCreated()),
                notificationService.formatInstantToReadableString(writtenPipelineRun.getUpdated()),
                String.valueOf(expectedQuotaRemaining),
                testUserDescription));
    // success is a void method
    doNothing()
        .when(pubsubService)
        .publishMessage(
            notificationConfiguration.projectId(),
            notificationConfiguration.topicId(),
            stringifiedJobFailedNotification);

    notificationService.configureAndSendPipelineRunFailedNotification(
        testJobId, testUserId, flightContext);

    // verify that the pubsub method was called
    verify(pubsubService, times(1))
        .publishMessage(
            notificationConfiguration.projectId(),
            notificationConfiguration.topicId(),
            stringifiedJobFailedNotification);
  }

  @Test
  void configureAndSendPipelineRunFailedNotificationNoUserQuota() throws IOException {
    Pipeline pipeline = pipelinesService.getPipelineById(1L);
    PipelineRun writtenPipelineRun =
        createCompletedPipelineRunInDb(pipeline, CommonPipelineRunStatusEnum.FAILED);

    // don't initialize user in user_quota table
    int expectedQuotaRemaining =
        quotasService.getPipelineQuota(pipeline.getName()).getDefaultQuota();

    when(flightContext.getFlightId()).thenReturn(testJobId.toString());
    RawlsServiceApiException rawlsServiceApiException =
        new RawlsServiceApiException(testErrorMessage);
    StepResult stepResultFailedWithException =
        new StepResult(StepStatus.STEP_RESULT_FAILURE_FATAL, rawlsServiceApiException);
    when(flightContext.getResult()).thenReturn(stepResultFailedWithException);

    String stringifiedJobFailedNotification =
        objectMapper.writeValueAsString(
            new TeaspoonsJobFailedNotification(
                testUserId,
                pipeline.getDisplayName(),
                testJobId.toString(),
                testErrorMessage,
                notificationService.formatInstantToReadableString(writtenPipelineRun.getCreated()),
                notificationService.formatInstantToReadableString(writtenPipelineRun.getUpdated()),
                String.valueOf(expectedQuotaRemaining),
                testUserDescription));
    // success is a void method
    doNothing()
        .when(pubsubService)
        .publishMessage(
            notificationConfiguration.projectId(),
            notificationConfiguration.topicId(),
            stringifiedJobFailedNotification);

    notificationService.configureAndSendPipelineRunFailedNotification(
        testJobId, testUserId, flightContext);

    // verify that the pubsub method was called
    verify(pubsubService, times(1))
        .publishMessage(
            notificationConfiguration.projectId(),
            notificationConfiguration.topicId(),
            stringifiedJobFailedNotification);
  }

  @Test
  void configureAndSendPipelineRunFailedNotificationWithoutException() throws IOException {
    Pipeline pipeline = pipelinesService.getPipelineById(1L);
    PipelineRun writtenPipelineRun =
        createCompletedPipelineRunInDb(pipeline, CommonPipelineRunStatusEnum.FAILED);

    // don't initialize user in user_quota table
    int expectedQuotaRemaining =
        quotasService.getPipelineQuota(pipeline.getName()).getDefaultQuota();

    when(flightContext.getFlightId()).thenReturn(testJobId.toString());
    StepResult stepResultFailedWithoutException =
        new StepResult(StepStatus.STEP_RESULT_FAILURE_FATAL);
    when(flightContext.getResult()).thenReturn(stepResultFailedWithoutException);

    String stringifiedJobFailedNotification =
        objectMapper.writeValueAsString(
            new TeaspoonsJobFailedNotification(
                testUserId,
                pipeline.getDisplayName(),
                testJobId.toString(),
                "Unknown error",
                notificationService.formatInstantToReadableString(writtenPipelineRun.getCreated()),
                notificationService.formatInstantToReadableString(writtenPipelineRun.getUpdated()),
                String.valueOf(expectedQuotaRemaining),
                testUserDescription));
    // success is a void method
    doNothing()
        .when(pubsubService)
        .publishMessage(
            notificationConfiguration.projectId(),
            notificationConfiguration.topicId(),
            stringifiedJobFailedNotification);

    notificationService.configureAndSendPipelineRunFailedNotification(
        testJobId, testUserId, flightContext);

    // verify that the pubsub method was called
    verify(pubsubService, times(1))
        .publishMessage(
            notificationConfiguration.projectId(),
            notificationConfiguration.topicId(),
            stringifiedJobFailedNotification);
  }

  @Test
  void configureAndSendPipelineRunFailedNotificationIOException() throws IOException {
    Pipeline pipeline = pipelinesService.getPipelineById(1L);
    createCompletedPipelineRunInDb(pipeline, CommonPipelineRunStatusEnum.FAILED);

    when(flightContext.getFlightId()).thenReturn(testJobId.toString());
    RawlsServiceApiException rawlsServiceApiException =
        new RawlsServiceApiException(testErrorMessage);
    StepResult stepResultFailedWithException =
        new StepResult(StepStatus.STEP_RESULT_FAILURE_FATAL, rawlsServiceApiException);
    when(flightContext.getResult()).thenReturn(stepResultFailedWithException);

    doThrow(new IOException()).when(pubsubService).publishMessage(any(), any(), any());

    // exception should be caught
    assertDoesNotThrow(
        () ->
            notificationService.configureAndSendPipelineRunFailedNotification(
                testJobId, testUserId, flightContext));
  }

  /**
   * Helper method for tests to create a completed pipeline run in the database.
   *
   * @param pipeline the pipeline to create the run for
   * @param statusEnum the status of the pipeline run
   * @return the completed pipeline run
   */
  private PipelineRun createCompletedPipelineRunInDb(
      Pipeline pipeline, CommonPipelineRunStatusEnum statusEnum) {
    PipelineRun completedPipelineRun =
        new PipelineRun(
            testJobId,
            testUserId,
            pipeline.getId(),
            pipeline.getToolVersion(),
            pipeline.getWorkspaceId(),
            pipeline.getWorkspaceBillingProject(),
            pipeline.getWorkspaceName(),
            pipeline.getWorkspaceStorageContainerName(),
            pipeline.getWorkspaceGoogleProject(),
            null, // timestamps auto generated by db
            null,
            statusEnum,
            testUserDescription,
            testQuotaConsumedByJob);

    return pipelineRunsRepository.save(completedPipelineRun);
  }
}
