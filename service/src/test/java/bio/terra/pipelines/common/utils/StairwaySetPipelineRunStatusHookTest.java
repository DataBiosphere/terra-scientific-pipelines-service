package bio.terra.pipelines.common.utils;

import static bio.terra.pipelines.testutils.TestUtils.createNewPipelineRunWithJobId;
import static org.junit.jupiter.api.Assertions.assertEquals;

import bio.terra.pipelines.db.entities.PipelineRun;
import bio.terra.pipelines.db.repositories.PipelineRunsRepository;
import bio.terra.pipelines.service.PipelineRunsService;
import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import bio.terra.pipelines.testutils.TestFlightContext;
import bio.terra.pipelines.testutils.TestUtils;
import bio.terra.stairway.FlightStatus;
import java.util.UUID;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.springframework.beans.factory.annotation.Autowired;

class StairwaySetPipelineRunStatusHookTest extends BaseEmbeddedDbTest {
  StairwaySetPipelineRunStatusHook stairwaySetPipelineRunStatusHook;
  @Autowired PipelineRunsService pipelineRunsService;

  @Autowired PipelineRunsRepository pipelineRunsRepository;

  private final UUID testJobId = TestUtils.TEST_NEW_UUID;

  @BeforeEach
  void setup() {
    stairwaySetPipelineRunStatusHook = new StairwaySetPipelineRunStatusHook(pipelineRunsService);
  }

  @Test
  void endFlight_success() throws InterruptedException {

    var context =
        new TestFlightContext()
            .flightClassName("bio.terra.testing.flight.TestFlight")
            .stepClassName("bio.terra.testing.StepClass");
    // write pipelineRun to the db
    PipelineRun pipelineRun = createNewPipelineRunWithJobId(UUID.fromString(context.getFlightId()));
    pipelineRunsRepository.save(pipelineRun);

    context.getInputParameters().put("user_id", TestUtils.TEST_USER_ID_1);

    stairwaySetPipelineRunStatusHook.startFlight(context);

    context.flightStatus(FlightStatus.SUCCESS);

    stairwaySetPipelineRunStatusHook.endFlight(context);

    PipelineRun writtenPipelineRun =
        pipelineRunsRepository.findByJobIdAndUserId(testJobId, TestUtils.TEST_USER_ID_1).get();
    assertEquals(CommonPipelineRunStatusEnum.SUCCEEDED, writtenPipelineRun.getStatus());
  }

  @Test
  void endFlight_error() throws InterruptedException {

    var context =
        new TestFlightContext()
            .flightClassName("bio.terra.testing.flight.TestFlight")
            .stepClassName("bio.terra.testing.StepClass");
    // write pipelineRun to the db
    PipelineRun pipelineRun = createNewPipelineRunWithJobId(UUID.fromString(context.getFlightId()));
    pipelineRunsRepository.save(pipelineRun);

    context.getInputParameters().put("user_id", TestUtils.TEST_USER_ID_1);

    stairwaySetPipelineRunStatusHook.startFlight(context);

    // make the flight fail
    context.flightStatus(FlightStatus.ERROR);

    stairwaySetPipelineRunStatusHook.endFlight(context);

    PipelineRun writtenPipelineRun =
        pipelineRunsRepository.findByJobIdAndUserId(testJobId, TestUtils.TEST_USER_ID_1).get();
    assertEquals(CommonPipelineRunStatusEnum.FAILED, writtenPipelineRun.getStatus());
  }
}
