package bio.terra.pipelines.stairway.flights.datadelivery;

import static org.junit.jupiter.api.Assertions.*;

import bio.terra.pipelines.common.utils.FlightBeanBag;
import bio.terra.pipelines.dependencies.stairway.JobMapKeys;
import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import bio.terra.pipelines.testutils.TestUtils;
import bio.terra.stairway.FlightMap;
import java.util.UUID;
import org.junit.jupiter.api.Test;
import org.springframework.beans.factory.annotation.Autowired;

class DeliverDataToGcsFlightTest extends BaseEmbeddedDbTest {

  @Autowired FlightBeanBag flightBeanBag;

  @Test
  void testCreateDeliverDataToGcsFlight() {
    FlightMap inputParameters = new FlightMap();
    inputParameters.put(JobMapKeys.USER_ID, TestUtils.TEST_USER_ID_1);
    inputParameters.put(JobMapKeys.PIPELINE_ID, TestUtils.TEST_PIPELINE_ID_1);
    inputParameters.put(
        JobMapKeys.PIPELINE_NAME, TestUtils.TEST_PIPELINE_1_IMPUTATION_ENUM.getValue());
    inputParameters.put(
        DataDeliveryJobMapKeys.DESTINATION_GCS_PATH, "gs://test-bucket/test-destination");
    inputParameters.put(DataDeliveryJobMapKeys.PIPELINE_RUN_ID, UUID.randomUUID());

    DeliverDataToGcsFlight flight = new DeliverDataToGcsFlight(inputParameters, flightBeanBag);

    assertNotNull(flight);

    // Verify it has the one expected step (will need to update this later on)
    assertEquals(1, flight.getSteps().size());
    assertEquals(
        "DeliverOutputFilesToGcsStep",
        flight.getSteps().get(0).getClass().getSimpleName(),
        "Flight should have DeliverOutputFilesToGcsStep");
  }
}
