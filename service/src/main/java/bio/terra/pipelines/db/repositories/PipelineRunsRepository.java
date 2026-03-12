package bio.terra.pipelines.db.repositories;

import bio.terra.pipelines.db.entities.PipelineRun;
import java.util.List;
import java.util.Optional;
import java.util.UUID;
import org.springframework.data.jpa.repository.JpaSpecificationExecutor;
import org.springframework.data.jpa.repository.Query;
import org.springframework.data.repository.CrudRepository;
import org.springframework.data.repository.PagingAndSortingRepository;
import org.springframework.data.repository.query.Param;

public interface PipelineRunsRepository
    extends CrudRepository<PipelineRun, Long>,
        PagingAndSortingRepository<PipelineRun, Long>,
        JpaSpecificationExecutor<PipelineRun> {
  List<PipelineRun> findAllByUserId(String userId);

  boolean existsByJobId(UUID jobId);

  Optional<PipelineRun> findByJobIdAndUserId(UUID jobId, String userId);

  /*
   Fetches PipelineRun with its Pipeline and PipelineOutputDefinitions in a single query
   to avoid LazyInitializationException when accessing pipeline.getPipelineOutputDefinitions()
   in getPipelineOutputsFileSize(). The LEFT JOIN FETCH ensures output definitions are loaded
   within the active Hibernate session before the transaction closes.
  */
  @Query(
      "SELECT pr FROM PipelineRun pr "
          + "JOIN FETCH pr.pipeline p "
          + "LEFT JOIN FETCH p.pipelineOutputDefinitions "
          + "WHERE pr.jobId = :jobId AND pr.userId = :userId")
  Optional<PipelineRun> findByJobIdAndUserIdWithPipeline(
      @Param("jobId") UUID jobId, @Param("userId") String userId);

  boolean existsByIdGreaterThan(Long id);

  int countByUserId(String userId);
}
