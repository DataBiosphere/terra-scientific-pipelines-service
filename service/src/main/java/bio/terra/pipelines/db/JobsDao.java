package bio.terra.pipelines.db;

import static java.sql.Types.TIMESTAMP;
import static java.sql.Types.VARCHAR;

import bio.terra.common.db.ReadTransaction;
import bio.terra.common.db.WriteTransaction;
import bio.terra.pipelines.app.configuration.TspsDatabaseConfiguration;
import bio.terra.pipelines.db.exception.DuplicateObjectException;
import bio.terra.pipelines.db.exception.JobNotFoundException;
import bio.terra.pipelines.service.model.Job;
import java.util.*;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.dao.DuplicateKeyException;
import org.springframework.dao.EmptyResultDataAccessException;
import org.springframework.dao.support.DataAccessUtils;
import org.springframework.jdbc.core.RowMapper;
import org.springframework.jdbc.core.namedparam.MapSqlParameterSource;
import org.springframework.jdbc.core.namedparam.NamedParameterJdbcTemplate;
import org.springframework.stereotype.Component;

@Component
public class JobsDao {
  private static final RowMapper<DbJob> DB_JOB_ROW_MAPPER =
      (rs, rowNum) ->
          new DbJob(
              UUID.fromString(rs.getString("job_id")),
              rs.getString("user_id"),
              rs.getString("pipeline_id"),
              rs.getString("pipeline_version"),
              rs.getTimestamp("time_submitted"),
              Optional.ofNullable(rs.getTimestamp("time_completed")),
              rs.getString("status"));

  private final NamedParameterJdbcTemplate jdbcTemplate;
  private final Logger logger = LoggerFactory.getLogger(JobsDao.class);

  @Autowired
  public JobsDao(TspsDatabaseConfiguration tspsDatabaseConfiguration) {
    this.jdbcTemplate = new NamedParameterJdbcTemplate(tspsDatabaseConfiguration.getDataSource());
  }

  /**
   * Writes a job to DB. Returns job UUID on success.
   *
   * @param job all properties of the job to create
   * @return jobUuid
   */
  @WriteTransaction
  public UUID createJob(Job job) {
    final String sql =
        """
                        INSERT INTO jobs (job_id, user_id, pipeline_id, pipeline_version, time_submitted, status)
                        VALUES (:jobId, :userId, :pipelineId, :pipelineVersion, :timeSubmitted, :status)
                        """;

    final UUID jobUuid = job.getJobId();

    MapSqlParameterSource params =
        new MapSqlParameterSource()
            .addValue("jobId", jobUuid.toString(), VARCHAR)
            .addValue("userId", job.getUserId(), VARCHAR)
            .addValue("pipelineId", job.getPipelineId(), VARCHAR)
            .addValue("pipelineVersion", job.getPipelineVersion(), VARCHAR)
            .addValue("timeSubmitted", job.getTimeSubmitted(), TIMESTAMP)
            .addValue("status", job.getStatus(), VARCHAR);
    try {
      jdbcTemplate.update(sql, params);
      logger.info("Inserted record for job {}", jobUuid);
    } catch (DuplicateKeyException e) {
      String message = e.getMessage();
      if (message != null
          && message.contains("duplicate key value violates unique constraint \"jobs_pkey\"")) {
        // Job with job_id already exists.
        // TODO if this happens, JobsService should retry with a new UUID instead. see TSPS-19
        throw new DuplicateObjectException(
            String.format(
                "Job with id %s already exists - pipelineId %s submitted on %s",
                jobUuid, job.getPipelineId(), job.getTimeSubmitted()),
            e);
      } else {
        throw e;
      }
    }
    return jobUuid;
  }

  /**
   * Return specified job submitted by the user
   *
   * @param userId unique identifier of the calling user
   * @param pipelineId identifier of the pipeline
   * @param jobId UUID identifier of the job
   * @return Job object
   */
  public Job getJob(String userId, String pipelineId, String jobId) {
    DbJob dbJob =
        getDbJobIfExists(userId, pipelineId, jobId)
            .orElseThrow(() -> new JobNotFoundException(String.format("Job {} not found.", jobId)));

    return Job.fromDb(dbJob);
  }

  /**
   * Return all jobs submitted by the user
   *
   * @param userId unique identifier of the calling user
   * @param pipelineId identifier of the pipeline
   * @return List of Job objects
   */
  @ReadTransaction
  public List<Job> getJobs(String userId, String pipelineId) {
    List<Job> jobList = new ArrayList<>();

    List<DbJob> dbJobList = getDbJobs(userId, pipelineId);

    for (DbJob dbJob : dbJobList) {
      jobList.add(Job.fromDb(dbJob));
    }

    return jobList;
  }

  @ReadTransaction
  private Optional<DbJob> getDbJobIfExists(String userId, String pipelineId, String jobId) {
    final String sql =
        """
                        SELECT job_id, user_id, pipeline_id, pipeline_version, time_submitted, time_completed, status
                        FROM jobs
                        WHERE user_id = :userId AND pipeline_id = :pipelineId AND job_id = :jobId
                        """;

    MapSqlParameterSource params =
        new MapSqlParameterSource()
            .addValue("userId", userId)
            .addValue("pipelineId", pipelineId)
            .addValue("jobId", jobId);

    try {
      DbJob result =
          DataAccessUtils.requiredSingleResult(jdbcTemplate.query(sql, params, DB_JOB_ROW_MAPPER));
      return Optional.of(result);
    } catch (EmptyResultDataAccessException e) {
      return Optional.empty();
    }
  }

  private List<DbJob> getDbJobs(String userId, String pipelineId) {
    final String sql =
        """
                    SELECT job_id, user_id, pipeline_id, pipeline_version, time_submitted, time_completed, status
                    FROM jobs
                    WHERE user_id = :userId AND pipeline_id = :pipelineId
                    """;

    MapSqlParameterSource params =
        new MapSqlParameterSource().addValue("userId", userId).addValue("pipelineId", pipelineId);

    return jdbcTemplate.query(sql, params, DB_JOB_ROW_MAPPER);
  }
}