package bio.terra.pipelines.common.utils;

import static org.junit.jupiter.api.Assertions.assertEquals;

import bio.terra.pipelines.app.configuration.external.CbasConfiguration;
import bio.terra.pipelines.app.configuration.internal.ImputationConfiguration;
import bio.terra.pipelines.app.configuration.internal.PipelinesCommonConfiguration;
import bio.terra.pipelines.dependencies.cbas.CbasService;
import bio.terra.pipelines.dependencies.leonardo.LeonardoService;
import bio.terra.pipelines.dependencies.rawls.RawlsService;
import bio.terra.pipelines.dependencies.sam.SamService;
import bio.terra.pipelines.dependencies.wds.WdsService;
import bio.terra.pipelines.dependencies.workspacemanager.WorkspaceManagerService;
import bio.terra.pipelines.notifications.NotificationService;
import bio.terra.pipelines.service.PipelineInputsOutputsService;
import bio.terra.pipelines.service.PipelineRunsService;
import bio.terra.pipelines.service.PipelinesService;
import bio.terra.pipelines.service.QuotasService;
import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import org.junit.jupiter.api.Test;
import org.mockito.Mock;
import org.springframework.beans.factory.annotation.Autowired;

class FlightBeanBagTest extends BaseEmbeddedDbTest {

  @Autowired private PipelinesService pipelinesService;
  @Autowired private PipelineRunsService pipelineRunsService;
  @Autowired private PipelineInputsOutputsService pipelineInputsOutputsService;
  @Autowired private SamService samService;
  @Autowired private LeonardoService leonardoService;
  @Autowired private WdsService wdsService;
  @Autowired private CbasService cbasService;
  @Autowired private WorkspaceManagerService workspaceManagerService;
  @Autowired private RawlsService rawlsService;
  @Autowired private QuotasService quotasService;

  @Mock
  private NotificationService
      notificationService; // mock because at startup tries to auto-create a topic

  @Autowired private ImputationConfiguration imputationConfiguration;
  @Autowired private CbasConfiguration cbasConfiguration;
  @Autowired private PipelinesCommonConfiguration pipelinesCommonConfiguration;

  @Test
  void testFlightBeanBag() {
    FlightBeanBag flightBeanBag =
        new FlightBeanBag(
            pipelinesService,
            pipelineRunsService,
            pipelineInputsOutputsService,
            samService,
            leonardoService,
            wdsService,
            cbasService,
            rawlsService,
            quotasService,
            notificationService,
            workspaceManagerService,
            imputationConfiguration,
            cbasConfiguration,
            pipelinesCommonConfiguration);
    assertEquals(pipelinesService, flightBeanBag.getPipelinesService());
    assertEquals(pipelineRunsService, flightBeanBag.getPipelineRunsService());
    assertEquals(pipelineInputsOutputsService, flightBeanBag.getPipelineInputsOutputsService());
    assertEquals(samService, flightBeanBag.getSamService());
    assertEquals(leonardoService, flightBeanBag.getLeonardoService());
    assertEquals(wdsService, flightBeanBag.getWdsService());
    assertEquals(cbasService, flightBeanBag.getCbasService());
    assertEquals(workspaceManagerService, flightBeanBag.getWorkspaceManagerService());
    assertEquals(rawlsService, flightBeanBag.getRawlsService());
    assertEquals(quotasService, flightBeanBag.getQuotasService());
    assertEquals(notificationService, flightBeanBag.getNotificationService());
    assertEquals(imputationConfiguration, flightBeanBag.getImputationConfiguration());
    assertEquals(cbasConfiguration, flightBeanBag.getCbasConfiguration());
    assertEquals(pipelinesCommonConfiguration, flightBeanBag.getPipelinesCommonConfiguration());
  }
}
