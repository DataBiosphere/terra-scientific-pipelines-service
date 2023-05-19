package bio.terra.pipelines.db;

import static java.sql.Types.*;

import bio.terra.common.db.ReadTransaction;
import bio.terra.common.db.WriteTransaction;
import bio.terra.pipelines.app.configuration.TspsDatabaseConfiguration;
import bio.terra.pipelines.db.exception.JobNotFoundException;
import bio.terra.pipelines.service.model.Job;
import java.sql.Timestamp;
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

  private static final String USER_ID_PARAM = "userId";
  private static final String JOB_ID_PARAM = "jobId";
  private static final String PIPELINE_ID_PARAM = "pipelineId";
  private static final String PIPELINE_VERSION_PARAM = "pipelineVersion";
  private static final String TIME_SUBMITTED_PARAM = "timeSubmitted";
  private static final String STATUS_PARAM = "status";
  private static final String PIPELINE_INPUTS_PARAM = "pipelineInputs";
  private static final RowMapper<DbJob> DB_JOB_ROW_MAPPER =
      (rs, rowNum) ->
          new DbJob(
              UUID.fromString(rs.getString("job_id")),
              rs.getString("user_id"),
              rs.getString("pipeline_id"),
              rs.getString("pipeline_version"),
              rs.getTimestamp("time_submitted").toInstant(),
              Optional.ofNullable(rs.getTimestamp("time_completed")).map(Timestamp::toInstant),
              rs.getString("status"),
              rs.getString("pipeline_inputs"));

  private final NamedParameterJdbcTemplate jdbcTemplate;
  private final Logger logger = LoggerFactory.getLogger(JobsDao.class);

  @Autowired
  public JobsDao(TspsDatabaseConfiguration tspsDatabaseConfiguration) {
    this.jdbcTemplate = new NamedParameterJdbcTemplate(tspsDatabaseConfiguration.getDataSource());
  }

  /**
   * Writes a job to DB. Returns job UUID on success, or null if we violate the primary key
   * uniqueness constraint.
   *
   * @param job all properties of the job to create
   * @return jobUuid or null
   */
  @WriteTransaction
  public UUID createJob(Job job) {
    final String sql =
        """
                        INSERT INTO jobs (job_id, user_id, pipeline_id, pipeline_version, time_submitted, status, pipeline_inputs)
                        VALUES (:jobId, :userId, :pipelineId, :pipelineVersion, :timeSubmitted, :status, :pipelineInputs)
                        """;

    final UUID jobUuid = job.getJobId();

    MapSqlParameterSource params =
        new MapSqlParameterSource()
            .addValue(JOB_ID_PARAM, jobUuid.toString(), VARCHAR)
            .addValue(USER_ID_PARAM, job.getUserId(), VARCHAR)
            .addValue(PIPELINE_ID_PARAM, job.getPipelineId(), VARCHAR)
            .addValue(PIPELINE_VERSION_PARAM, job.getPipelineVersion(), VARCHAR)
            .addValue(TIME_SUBMITTED_PARAM, Timestamp.from(job.getTimeSubmitted()), TIMESTAMP)
            .addValue(STATUS_PARAM, job.getStatus(), VARCHAR)
            .addValue(PIPELINE_INPUTS_PARAM, job.getPipelineInputs(), VARCHAR);

    try {
      jdbcTemplate.update(sql, params);
      logger.info("Inserted record for job {}", jobUuid);
    } catch (DuplicateKeyException e) {
      String message = e.getMessage();

      // Check to see if the message contains a reference to a constraint on just the job_id field.
      // If so, it's the primary key and we can retry it
      if (message != null && message.toLowerCase().contains("(job_id)")) {
        // Job with job_id already exists.
        logger.warn("Duplicate jobId {} unable to be written to database", jobUuid);
        return null;
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
            .orElseThrow(() -> new JobNotFoundException(String.format("Job %s not found.", jobId)));

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
                        SELECT job_id, user_id, pipeline_id, pipeline_version, time_submitted, time_completed, status, pipeline_inputs
                        FROM jobs
                        WHERE user_id = :userId AND pipeline_id = :pipelineId AND job_id = :jobId
                        """;

    MapSqlParameterSource params =
        new MapSqlParameterSource()
            .addValue(USER_ID_PARAM, userId)
            .addValue(PIPELINE_ID_PARAM, pipelineId)
            .addValue(JOB_ID_PARAM, jobId);

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
                    SELECT job_id, user_id, pipeline_id, pipeline_version, time_submitted, time_completed, status, pipeline_inputs
                    FROM jobs
                    WHERE user_id = :userId AND pipeline_id = :pipelineId
                    """;

    MapSqlParameterSource params =
        new MapSqlParameterSource()
            .addValue(USER_ID_PARAM, userId)
            .addValue(PIPELINE_ID_PARAM, pipelineId);

    return jdbcTemplate.query(sql, params, DB_JOB_ROW_MAPPER);
  }
}
