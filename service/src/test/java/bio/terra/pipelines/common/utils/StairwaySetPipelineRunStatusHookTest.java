package bio.terra.pipelines.common.utils;

import static bio.terra.pipelines.testutils.TestUtils.createNewPipelineRunWithJobId;
import static org.junit.jupiter.api.Assertions.assertDoesNotThrow;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotEquals;

import bio.terra.pipelines.db.entities.PipelineRun;
import bio.terra.pipelines.db.repositories.PipelineRunsRepository;
import bio.terra.pipelines.service.PipelineRunsService;
import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import bio.terra.pipelines.testutils.StairwayTestUtils;
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
  void endFlight_notPipelineRunTypeFlight_success() throws InterruptedException {
    var context =
        new TestFlightContext()
            .flightClassName("bio.terra.testing.flight.TestFlight")
            .stepClassName("bio.terra.testing.StepClass"); // stepClassName doesn't matter

    stairwaySetPipelineRunStatusHook.startFlight(context);
    // logic should not be executed because this TestFlight does not contain the
    // do_set_pipeline_run_status_failed_hook key in the inputParameters; if it tried to run the
    // logic on this flight, it would throw an exception because there is no "user_id" key in the
    // inputParameters
    assertDoesNotThrow(() -> stairwaySetPipelineRunStatusHook.endFlight(context));
  }

  @Test
  void endFlight_pipelineRunTypeFlight_success() throws InterruptedException {

    var context =
        new TestFlightContext()
            .flightClassName("bio.terra.pipelines.stairway.imputation.RunImputationGcpJobFlight")
            .stepClassName("bio.terra.testing.StepClass"); // stepClassName doesn't matter
    // write pipelineRun to the db
    PipelineRun pipelineRun = createNewPipelineRunWithJobId(UUID.fromString(context.getFlightId()));
    pipelineRunsRepository.save(pipelineRun);

    // this includes setting the DO_SET_PIPELINE_RUN_STATUS_FAILED_HOOK key to true
    StairwayTestUtils.constructCreateJobInputs(context.getInputParameters());

    stairwaySetPipelineRunStatusHook.startFlight(context);

    context.flightStatus(FlightStatus.SUCCESS);

    stairwaySetPipelineRunStatusHook.endFlight(context);

    // the flight did not fail, so the pipelineRun status should not be set to FAILED
    PipelineRun writtenPipelineRun =
        pipelineRunsRepository.findByJobIdAndUserId(testJobId, TestUtils.TEST_USER_ID_1).get();
    assertNotEquals(CommonPipelineRunStatusEnum.FAILED, writtenPipelineRun.getStatus());
  }

  @Test
  void endFlight_notPipelineRunTypeFlight_error() throws InterruptedException {
    var context =
        new TestFlightContext()
            .flightClassName("bio.terra.testing.flight.TestFlight")
            .stepClassName("bio.terra.testing.StepClass"); // stepClassName doesn't matter

    stairwaySetPipelineRunStatusHook.startFlight(context);
    // make the flight fail
    context.flightStatus(FlightStatus.ERROR);

    // logic should not be executed because this TestFlight does not contain the
    // do_set_pipeline_run_status_failed_hook key in the inputParameters; if it tried to run the
    // logic on this flight, it would throw an exception because there is no "user_id" key in the
    // inputParameters
    assertDoesNotThrow(() -> stairwaySetPipelineRunStatusHook.endFlight(context));
  }

  @Test
  void endFlight_PipelineRunTypeFlight_error() throws InterruptedException {

    var context =
        new TestFlightContext()
            .flightClassName("bio.terra.pipelines.stairway.imputation.RunImputationGcpJobFlight")
            .stepClassName("bio.terra.testing.StepClass"); // stepClassName doesn't matter
    // write pipelineRun to the db
    PipelineRun pipelineRun = createNewPipelineRunWithJobId(UUID.fromString(context.getFlightId()));
    pipelineRunsRepository.save(pipelineRun);

    // this includes setting the DO_SET_PIPELINE_RUN_STATUS_FAILED_HOOK key to true
    StairwayTestUtils.constructCreateJobInputs(context.getInputParameters());

    stairwaySetPipelineRunStatusHook.startFlight(context);

    // make the flight fail
    context.flightStatus(FlightStatus.ERROR);
    assertEquals(FlightStatus.ERROR, context.getFlightStatus());

    stairwaySetPipelineRunStatusHook.endFlight(context);

    PipelineRun writtenPipelineRun =
        pipelineRunsRepository.findByJobIdAndUserId(testJobId, TestUtils.TEST_USER_ID_1).get();
    assertEquals(CommonPipelineRunStatusEnum.FAILED, writtenPipelineRun.getStatus());
  }
}
