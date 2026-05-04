package bio.terra.pipelines.service;

import static org.junit.jupiter.api.Assertions.*;

import bio.terra.common.exception.NotFoundException;
import bio.terra.pipelines.common.GcsFile;
import bio.terra.pipelines.common.utils.CommonPipelineRunStatusEnum;
import bio.terra.pipelines.common.utils.DataDeliveryStatusEnum;
import bio.terra.pipelines.db.entities.DataDelivery;
import bio.terra.pipelines.db.entities.PipelineRun;
import bio.terra.pipelines.db.repositories.DataDeliveryRepository;
import bio.terra.pipelines.db.repositories.PipelineRunsRepository;
import bio.terra.pipelines.testutils.BaseEmbeddedDbTest;
import bio.terra.pipelines.testutils.TestUtils;
import java.util.UUID;
import org.junit.jupiter.api.Test;
import org.springframework.beans.factory.annotation.Autowired;

class DataDeliveryServiceTest extends BaseEmbeddedDbTest {
  @Autowired DataDeliveryService dataDeliveryService;
  @Autowired DataDeliveryRepository dataDeliveryRepository;
  @Autowired PipelineRunsRepository pipelineRunsRepository;

  @Test
  void createDataDelivery() {
    PipelineRun pipelineRun = createTestPipelineRun();

    GcsFile gcsPath = new GcsFile("gs://test-bucket/test-path");

    DataDelivery dataDelivery =
        dataDeliveryService.createDataDelivery(
            pipelineRun.getId(), UUID.randomUUID(), DataDeliveryStatusEnum.RUNNING, gcsPath);

    DataDelivery savedDelivery = dataDeliveryRepository.findById(dataDelivery.getId()).orElse(null);
    assertNotNull(savedDelivery);
    assertEquals(dataDelivery.getId(), savedDelivery.getId());
    assertEquals(dataDelivery.getStatus(), savedDelivery.getStatus());
  }

  @Test
  void getLatestDataDeliveryByPipelineRunId() {
    PipelineRun pipelineRun = createTestPipelineRun();

    // Create multiple data delivery records for the same pipeline run
    dataDeliveryService.createDataDelivery(
        pipelineRun.getId(),
        UUID.randomUUID(),
        DataDeliveryStatusEnum.RUNNING,
        new GcsFile("gs://test-bucket/path1"));

    DataDelivery delivery2 =
        dataDeliveryService.createDataDelivery(
            pipelineRun.getId(),
            UUID.randomUUID(),
            DataDeliveryStatusEnum.SUCCEEDED,
            new GcsFile("gs://test-bucket/path2"));

    // Get the latest delivery
    DataDelivery latestDelivery =
        dataDeliveryService.getLatestDataDeliveryByPipelineRunId(pipelineRun.getId());

    // Verify we got the most recent one
    assertNotNull(latestDelivery);
    assertEquals(delivery2.getId(), latestDelivery.getId());
    assertEquals(DataDeliveryStatusEnum.SUCCEEDED, latestDelivery.getStatus());
    assertEquals("gs://test-bucket/path2", latestDelivery.getGcsDestinationPath());
  }

  @Test
  void getLatestDataDeliveryByPipelineRunIdReturnsNullWhenNoneExist() {
    PipelineRun pipelineRun = createTestPipelineRun();

    DataDelivery result =
        dataDeliveryService.getLatestDataDeliveryByPipelineRunId(pipelineRun.getId());

    assertNull(result);
  }

  @Test
  void updateDataDeliveryStatus() {
    PipelineRun pipelineRun = createTestPipelineRun();

    dataDeliveryService.createDataDelivery(
        pipelineRun.getId(),
        UUID.randomUUID(),
        DataDeliveryStatusEnum.RUNNING,
        new GcsFile("gs://test-bucket/test-path"));

    DataDelivery updatedDelivery =
        dataDeliveryService.updateDataDeliveryStatus(
            pipelineRun.getId(), DataDeliveryStatusEnum.SUCCEEDED);

    DataDelivery persistedDelivery =
        dataDeliveryRepository.findById(updatedDelivery.getId()).orElse(null);
    assertNotNull(persistedDelivery);
    assertEquals(DataDeliveryStatusEnum.SUCCEEDED, persistedDelivery.getStatus());
  }

  @Test
  void updateDataDeliveryStatusThrowsNotFoundExceptionWhenNoRecordExists() {
    PipelineRun pipelineRun = createTestPipelineRun();

    // Attempt to update status when no delivery exists
    Long pipelineRunId = pipelineRun.getId();
    NotFoundException exception =
        assertThrows(
            NotFoundException.class,
            () ->
                dataDeliveryService.updateDataDeliveryStatus(
                    pipelineRunId, DataDeliveryStatusEnum.SUCCEEDED));

    // Verify the exception message
    assertTrue(exception.getMessage().contains("No data delivery record found"));
    assertTrue(exception.getMessage().contains(pipelineRun.getId().toString()));
  }

  private PipelineRun createTestPipelineRun() {
    PipelineRun pipelineRun = new PipelineRun();
    pipelineRun.setJobId(UUID.randomUUID());
    pipelineRun.setUserId(TestUtils.TEST_USER_1_ID);
    // Keep pipeline id aligned with seeded metadata + pipelines rows used by FK constraints.
    pipelineRun.setPipelineId(TestUtils.TEST_PIPELINE_ID_1);
    pipelineRun.setStatus(CommonPipelineRunStatusEnum.SUCCEEDED);
    pipelineRun.setWorkspaceBillingProject(TestUtils.CONTROL_WORKSPACE_BILLING_PROJECT);
    pipelineRun.setWorkspaceName(TestUtils.CONTROL_WORKSPACE_NAME);
    pipelineRun.setWorkspaceStorageContainerName(TestUtils.CONTROL_WORKSPACE_CONTAINER_NAME);
    pipelineRun.setWorkspaceGoogleProject(TestUtils.CONTROL_WORKSPACE_GOOGLE_PROJECT);
    return pipelineRunsRepository.save(pipelineRun);
  }
}
