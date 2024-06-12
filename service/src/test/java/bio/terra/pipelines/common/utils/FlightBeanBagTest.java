package bio.terra.pipelines.common.utils;

import static org.junit.jupiter.api.Assertions.assertEquals;

import bio.terra.pipelines.app.configuration.external.CbasConfiguration;
import bio.terra.pipelines.app.configuration.internal.ImputationConfiguration;
import bio.terra.pipelines.dependencies.cbas.CbasService;
import bio.terra.pipelines.dependencies.leonardo.LeonardoService;
import bio.terra.pipelines.dependencies.sam.SamService;
import bio.terra.pipelines.dependencies.wds.WdsService;
import bio.terra.pipelines.dependencies.workspacemanager.WorkspaceService;
import bio.terra.pipelines.service.PipelineRunsService;
import bio.terra.pipelines.service.PipelinesService;
import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import org.junit.jupiter.api.Test;
import org.springframework.beans.factory.annotation.Autowired;

class FlightBeanBagTest extends BaseEmbeddedDbTest {

  @Autowired private PipelinesService pipelinesService;
  @Autowired private PipelineRunsService pipelineRunsService;
  @Autowired private SamService samService;
  @Autowired private LeonardoService leonardoService;
  @Autowired private WdsService wdsService;
  @Autowired private CbasService cbasService;
  @Autowired private WorkspaceService workspaceService;
  @Autowired private ImputationConfiguration imputationConfiguration;
  @Autowired private CbasConfiguration cbasConfiguration;

  @Test
  void testFlightBeanBag() {
    FlightBeanBag flightBeanBag =
        new FlightBeanBag(
            pipelinesService,
            pipelineRunsService,
            samService,
            leonardoService,
            wdsService,
            cbasService,
            workspaceService,
            imputationConfiguration,
            cbasConfiguration);
    assertEquals(pipelinesService, flightBeanBag.getPipelinesService());
    assertEquals(pipelineRunsService, flightBeanBag.getPipelineRunsService());
    assertEquals(samService, flightBeanBag.getSamService());
    assertEquals(leonardoService, flightBeanBag.getLeonardoService());
    assertEquals(wdsService, flightBeanBag.getWdsService());
    assertEquals(cbasService, flightBeanBag.getCbasService());
    assertEquals(workspaceService, flightBeanBag.getWorkspaceService());
    assertEquals(imputationConfiguration, flightBeanBag.getImputationConfiguration());
    assertEquals(cbasConfiguration, flightBeanBag.getCbasConfiguration());
  }
}
