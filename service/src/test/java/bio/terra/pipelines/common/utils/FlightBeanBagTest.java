package bio.terra.pipelines.common.utils;

import static org.junit.jupiter.api.Assertions.assertEquals;

import bio.terra.pipelines.app.configuration.internal.PipelineConfigurations;
import bio.terra.pipelines.dependencies.rawls.RawlsService;
import bio.terra.pipelines.dependencies.sam.SamService;
import bio.terra.pipelines.notifications.NotificationService;
import bio.terra.pipelines.service.*;
import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import org.junit.jupiter.api.Test;
import org.mockito.Mock;
import org.springframework.beans.factory.annotation.Autowired;

class FlightBeanBagTest extends BaseEmbeddedDbTest {

  @Autowired private PipelinesService pipelinesService;
  @Autowired private PipelineRunsService pipelineRunsService;
  @Autowired private PipelineInputsOutputsService pipelineInputsOutputsService;
  @Autowired private SamService samService;
  @Autowired private RawlsService rawlsService;
  @Autowired private QuotasService quotasService;
  @Autowired private DownloadCallCounterService downloadCallCounterService;

  @Mock
  private NotificationService
      notificationService; // mock because at startup tries to auto-create a topic

  @Autowired private PipelineConfigurations pipelineConfigurations;
  @Autowired private DataDeliveryService dataDeliveryService;

  @Test
  void testFlightBeanBag() {
    FlightBeanBag flightBeanBag =
        new FlightBeanBag(
            pipelinesService,
            pipelineRunsService,
            pipelineInputsOutputsService,
            dataDeliveryService,
            samService,
            rawlsService,
            quotasService,
            notificationService,
            downloadCallCounterService,
            pipelineConfigurations);
    assertEquals(pipelinesService, flightBeanBag.getPipelinesService());
    assertEquals(pipelineRunsService, flightBeanBag.getPipelineRunsService());
    assertEquals(pipelineInputsOutputsService, flightBeanBag.getPipelineInputsOutputsService());
    assertEquals(samService, flightBeanBag.getSamService());
    assertEquals(rawlsService, flightBeanBag.getRawlsService());
    assertEquals(quotasService, flightBeanBag.getQuotasService());
    assertEquals(notificationService, flightBeanBag.getNotificationService());
    assertEquals(downloadCallCounterService, flightBeanBag.getDownloadCallCounterService());
    assertEquals(pipelineConfigurations, flightBeanBag.getPipelineConfigurations());
  }
}
