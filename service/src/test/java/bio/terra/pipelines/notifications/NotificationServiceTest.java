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
import org.springframework.boot.test.mock.mockito.MockBean;

class NotificationServiceTest extends BaseEmbeddedDbTest {
  @InjectMocks @Autowired NotificationService notificationService;
  @Autowired PipelineRunsService pipelineRunsService;
  @Autowired PipelineRunsRepository pipelineRunsRepository;
  @Autowired PipelinesService pipelinesService;
  @Autowired QuotasService quotasService;
  @Autowired NotificationConfiguration notificationConfiguration;
  @Autowired ObjectMapper objectMapper;
  @MockBean PubsubService pubsubService;
  @Mock private FlightContext flightContext;

  UUID testJobId = TestUtils.TEST_NEW_UUID;
  String testUserId = TestUtils.TEST_USER_ID_1;
  Integer testQuotaConsumedByJob = 1000;
  String testUserDescription = TestUtils.TEST_USER_PROVIDED_DESCRIPTION;
  String testErrorMessage = "test error message";

  @Test
  void createBaseTeaspoonsJobNotification() {
    Pipeline pipeline = pipelinesService.getPipelineById(1L);
    PipelineRun writtenPipelineRun =
        createCompletedPipelineRunInDb(pipeline, CommonPipelineRunStatusEnum.SUCCEEDED);

    // initialize and set user quota
    UserQuota userQuota =
        quotasService.getOrCreateQuotaForUserAndPipeline(testUserId, pipeline.getName());
    UserQuota updatedUserQuota =
        quotasService.updateQuotaConsumed(userQuota, testQuotaConsumedByJob);
    int expectedQuotaRemaining = updatedUserQuota.getQuota() - updatedUserQuota.getQuotaConsumed();

    BaseTeaspoonsJobNotification baseTeaspoonsJobNotification =
        notificationService.createBaseTeaspoonsJobNotification(testJobId, testUserId);

    assertEquals(testUserId, baseTeaspoonsJobNotification.getRecipientUserId());
    assertEquals(pipeline.getDisplayName(), baseTeaspoonsJobNotification.getPipelineDisplayName());
    assertEquals(testJobId.toString(), baseTeaspoonsJobNotification.getJobId());
    assertEquals(
        notificationService.formatInstantToReadableString(writtenPipelineRun.getCreated()),
        baseTeaspoonsJobNotification.getTimeSubmitted());
    assertEquals(
        notificationService.formatInstantToReadableString(writtenPipelineRun.getUpdated()),
        baseTeaspoonsJobNotification.getTimeCompleted());
    assertEquals(
        String.valueOf(expectedQuotaRemaining), baseTeaspoonsJobNotification.getQuotaRemaining());
    assertEquals(testUserDescription, baseTeaspoonsJobNotification.getUserDescription());
  }

  @Test
  void createBaseTeaspoonsJobNotificationNoQuota() {
    Pipeline pipeline = pipelinesService.getPipelineById(1L);
    createCompletedPipelineRunInDb(pipeline, CommonPipelineRunStatusEnum.SUCCEEDED);

    // don't initialize user in user_quota table
    int expectedQuotaRemaining =
        quotasService.getPipelineQuota(pipeline.getName()).getDefaultQuota();

    BaseTeaspoonsJobNotification baseTeaspoonsJobNotification =
        notificationService.createBaseTeaspoonsJobNotification(testJobId, testUserId);

    assertEquals(testUserId, baseTeaspoonsJobNotification.getRecipientUserId());
    assertEquals(
        String.valueOf(expectedQuotaRemaining), baseTeaspoonsJobNotification.getQuotaRemaining());
  }

  @Test
  void formatDateTime() {
    Instant instant = Instant.parse("2021-08-25T12:34:56.789Z");
    String formattedDateTime = notificationService.formatInstantToReadableString(instant);
    assertEquals("Wed, 25 Aug 2021 12:34:56 GMT", formattedDateTime);
  }

  @Test
  void createTeaspoonsJobFailedNotification() {
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

    TeaspoonsJobFailedNotification teaspoonsJobFailedNotification =
        notificationService.createTeaspoonsJobFailedNotification(
            testJobId, testUserId, flightContext);

    assertEquals(testUserId, teaspoonsJobFailedNotification.getRecipientUserId());
    assertEquals(
        pipeline.getDisplayName(), teaspoonsJobFailedNotification.getPipelineDisplayName());
    assertEquals(testJobId.toString(), teaspoonsJobFailedNotification.getJobId());
    assertEquals(
        notificationService.formatInstantToReadableString(writtenPipelineRun.getCreated()),
        teaspoonsJobFailedNotification.getTimeSubmitted());
    assertEquals(
        notificationService.formatInstantToReadableString(writtenPipelineRun.getUpdated()),
        teaspoonsJobFailedNotification.getTimeCompleted());
    assertEquals(
        String.valueOf(expectedQuotaRemaining), teaspoonsJobFailedNotification.getQuotaRemaining());
    assertEquals(testUserDescription, teaspoonsJobFailedNotification.getUserDescription());

    // fields specific to failed job notification
    assertEquals(
        "TeaspoonsJobFailedNotification", teaspoonsJobFailedNotification.getNotificationType());
    assertEquals(testErrorMessage, teaspoonsJobFailedNotification.getErrorMessage());
    assertEquals("0", teaspoonsJobFailedNotification.getQuotaConsumedByJob());
  }

  @Test
  void createTeaspoonsJobFailedNotificationNoMessage() {
    Pipeline pipeline = pipelinesService.getPipelineById(1L);
    createCompletedPipelineRunInDb(pipeline, CommonPipelineRunStatusEnum.FAILED);

    when(flightContext.getFlightId()).thenReturn(testJobId.toString());
    StepResult stepResultFailedNoException = new StepResult(StepStatus.STEP_RESULT_FAILURE_FATAL);
    when(flightContext.getResult()).thenReturn(stepResultFailedNoException);

    TeaspoonsJobFailedNotification teaspoonsJobFailedNotification =
        notificationService.createTeaspoonsJobFailedNotification(
            testJobId, testUserId, flightContext);

    assertEquals(testUserId, teaspoonsJobFailedNotification.getRecipientUserId());

    // test for message when no exception is present
    assertEquals(
        "TeaspoonsJobFailedNotification", teaspoonsJobFailedNotification.getNotificationType());
    assertEquals("Unknown error", teaspoonsJobFailedNotification.getErrorMessage());
  }

  @Test
  void createTeaspoonsJobSucceededNotification() {
    Pipeline pipeline = pipelinesService.getPipelineById(1L);
    PipelineRun writtenPipelineRun =
        createCompletedPipelineRunInDb(pipeline, CommonPipelineRunStatusEnum.SUCCEEDED);

    // initialize and set user quota
    UserQuota userQuota =
        quotasService.getOrCreateQuotaForUserAndPipeline(testUserId, pipeline.getName());
    UserQuota updatedUserQuota =
        quotasService.updateQuotaConsumed(userQuota, testQuotaConsumedByJob);
    int expectedQuotaRemaining = updatedUserQuota.getQuota() - updatedUserQuota.getQuotaConsumed();

    TeaspoonsJobSucceededNotification teaspoonsJobSucceededNotification =
        notificationService.createTeaspoonsJobSucceededNotification(testJobId, testUserId);

    assertEquals(testUserId, teaspoonsJobSucceededNotification.getRecipientUserId());
    assertEquals(
        pipeline.getDisplayName(), teaspoonsJobSucceededNotification.getPipelineDisplayName());
    assertEquals(testJobId.toString(), teaspoonsJobSucceededNotification.getJobId());
    assertEquals(
        notificationService.formatInstantToReadableString(writtenPipelineRun.getCreated()),
        teaspoonsJobSucceededNotification.getTimeSubmitted());
    assertEquals(
        notificationService.formatInstantToReadableString(writtenPipelineRun.getUpdated()),
        teaspoonsJobSucceededNotification.getTimeCompleted());
    assertEquals(
        String.valueOf(expectedQuotaRemaining),
        teaspoonsJobSucceededNotification.getQuotaRemaining());
    assertEquals(testUserDescription, teaspoonsJobSucceededNotification.getUserDescription());

    // fields specific to succeeded job notification
    assertEquals(
        "TeaspoonsJobSucceededNotification",
        teaspoonsJobSucceededNotification.getNotificationType());
    assertEquals(
        testQuotaConsumedByJob.toString(),
        teaspoonsJobSucceededNotification.getQuotaConsumedByJob());
  }

  @Test
  void configureAndSendPipelineRunSucceededNotification() throws IOException, InterruptedException {
    Pipeline pipeline = pipelinesService.getPipelineById(1L);
    createCompletedPipelineRunInDb(pipeline, CommonPipelineRunStatusEnum.SUCCEEDED);

    // success is a void method
    doNothing().when(pubsubService).publishMessage(any(), any(), any());

    notificationService.configureAndSendPipelineRunSucceededNotification(testJobId, testUserId);

    String stringifiedJobSucceededNotification =
        objectMapper.writeValueAsString(
            notificationService.createTeaspoonsJobSucceededNotification(testJobId, testUserId));
    // verify that the pubsub method was called
    verify(pubsubService, times(1))
        .publishMessage(
            notificationConfiguration.projectId(),
            notificationConfiguration.topicId(),
            stringifiedJobSucceededNotification);
  }

  @Test
  void configureAndSendPipelineRunSucceededNotificationIOException()
      throws IOException, InterruptedException {
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
  void configureAndSendPipelineRunSucceededNotificationInterruptedException()
      throws IOException, InterruptedException {
    Pipeline pipeline = pipelinesService.getPipelineById(1L);
    createCompletedPipelineRunInDb(pipeline, CommonPipelineRunStatusEnum.SUCCEEDED);

    doThrow(new InterruptedException()).when(pubsubService).publishMessage(any(), any(), any());

    // exception should be caught
    assertDoesNotThrow(
        () ->
            notificationService.configureAndSendPipelineRunSucceededNotification(
                testJobId, testUserId));
  }

  @Test
  void configureAndSendPipelineRunFailedNotification() throws IOException, InterruptedException {
    Pipeline pipeline = pipelinesService.getPipelineById(1L);
    createCompletedPipelineRunInDb(pipeline, CommonPipelineRunStatusEnum.FAILED);

    when(flightContext.getFlightId()).thenReturn(testJobId.toString());
    RawlsServiceApiException rawlsServiceApiException =
        new RawlsServiceApiException(testErrorMessage);
    StepResult stepResultFailedWithException =
        new StepResult(StepStatus.STEP_RESULT_FAILURE_FATAL, rawlsServiceApiException);
    when(flightContext.getResult()).thenReturn(stepResultFailedWithException);

    // success is a void method
    doNothing().when(pubsubService).publishMessage(any(), any(), any());

    notificationService.configureAndSendPipelineRunFailedNotification(
        testJobId, testUserId, flightContext);

    String stringifiedJobFailedNotification =
        objectMapper.writeValueAsString(
            notificationService.createTeaspoonsJobFailedNotification(
                testJobId, testUserId, flightContext));
    // verify that the pubsub method was called
    verify(pubsubService, times(1))
        .publishMessage(
            notificationConfiguration.projectId(),
            notificationConfiguration.topicId(),
            stringifiedJobFailedNotification);
  }

  @Test
  void configureAndSendPipelineRunFailedNotificationIOException()
      throws IOException, InterruptedException {
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

  @Test
  void configureAndSendPipelineRunFailedNotificationInterruptedException()
      throws IOException, InterruptedException {
    Pipeline pipeline = pipelinesService.getPipelineById(1L);
    createCompletedPipelineRunInDb(pipeline, CommonPipelineRunStatusEnum.FAILED);

    when(flightContext.getFlightId()).thenReturn(testJobId.toString());
    RawlsServiceApiException rawlsServiceApiException =
        new RawlsServiceApiException(testErrorMessage);
    StepResult stepResultFailedWithException =
        new StepResult(StepStatus.STEP_RESULT_FAILURE_FATAL, rawlsServiceApiException);
    when(flightContext.getResult()).thenReturn(stepResultFailedWithException);

    doThrow(new InterruptedException()).when(pubsubService).publishMessage(any(), any(), any());

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
            pipeline.getWdlMethodVersion(),
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
