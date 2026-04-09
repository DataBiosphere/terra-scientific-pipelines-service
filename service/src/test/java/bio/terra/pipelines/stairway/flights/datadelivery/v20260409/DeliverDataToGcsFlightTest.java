package bio.terra.pipelines.stairway.flights.datadelivery.v20260409;

import static org.junit.jupiter.api.Assertions.*;

import bio.terra.pipelines.common.utils.FlightBeanBag;
import bio.terra.pipelines.dependencies.stairway.JobMapKeys;
import bio.terra.pipelines.stairway.flights.datadelivery.DataDeliveryJobMapKeys;
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
    inputParameters.put(JobMapKeys.USER_ID, TestUtils.TEST_USER_1_ID);
    inputParameters.put(
        DataDeliveryJobMapKeys.DESTINATION_GCS_PATH, "gs://test-bucket/test-destination");
    inputParameters.put(DataDeliveryJobMapKeys.PIPELINE_RUN_ID, UUID.randomUUID());

    DeliverDataToGcsFlight flight = new DeliverDataToGcsFlight(inputParameters, flightBeanBag);

    assertNotNull(flight);

    // Verify it has the expected steps
    assertEquals(4, flight.getSteps().size());
    assertEquals(
        "CreateDataDeliveryRecordStep",
        flight.getSteps().get(0).getClass().getSimpleName(),
        "Flight should have CreateDataDeliveryRecordStep");

    assertEquals(
        "DeliverOutputFilesToGcsStep",
        flight.getSteps().get(1).getClass().getSimpleName(),
        "Flight should have DeliverOutputFilesToGcsStep");

    assertEquals(
        "UpdateDataDeliveryStatusStep",
        flight.getSteps().get(2).getClass().getSimpleName(),
        "Flight should have UpdateDataDeliveryStatusStep");

    assertEquals(
        "DeleteOutputSourceFilesStep",
        flight.getSteps().get(3).getClass().getSimpleName(),
        "Flight should have DeleteOutputSourceFilesStep");
  }
}
