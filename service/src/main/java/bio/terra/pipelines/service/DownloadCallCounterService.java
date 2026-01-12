package bio.terra.pipelines.service;

import bio.terra.pipelines.db.entities.DownloadCallCount;
import bio.terra.pipelines.db.repositories.DownloadCallCounterRepository;
import java.util.Optional;
import java.util.UUID;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.stereotype.Service;

/* Service to track and increment the number of times output has been downloaded for a given job.
 * We use calls to the /output/signed-urls endpoint as a proxy for download attempts.
 */
@Service
public class DownloadCallCounterService {

  private static final Logger logger = LoggerFactory.getLogger(DownloadCallCounterService.class);
  private final DownloadCallCounterRepository downloadCallCounterRepository;

  public DownloadCallCounterService(DownloadCallCounterRepository downloadCallCounterRepository) {
    this.downloadCallCounterRepository = downloadCallCounterRepository;
  }

  /**
   * Increment and return the download call count for a given jobId. If the jobId doesn't exist in
   * the table, insert a new row for that job.
   */
  public void incrementDownloadCallCount(UUID jobId) {
    Optional<DownloadCallCount> existingRow = downloadCallCounterRepository.findByJobId(jobId);
    DownloadCallCount rowToIncrement = existingRow.orElseGet(() -> new DownloadCallCount(jobId, 0));
    int incrementedCount = rowToIncrement.getCount() + 1;
    rowToIncrement.setCount(incrementedCount);
    downloadCallCounterRepository.save(rowToIncrement);

    logger.info("JobId {} download call count incremented to {}", jobId, incrementedCount);
  }
}
