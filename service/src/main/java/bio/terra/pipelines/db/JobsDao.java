package bio.terra.pipelines.db;

import static java.sql.Types.TIMESTAMP;
import static java.sql.Types.VARCHAR;

import bio.terra.common.db.WriteTransaction;
import bio.terra.pipelines.app.configuration.TspsDatabaseConfiguration;
import bio.terra.pipelines.db.exception.DuplicateObjectException;
import bio.terra.pipelines.service.model.Job;
import java.util.*;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.dao.DuplicateKeyException;
import org.springframework.jdbc.core.JdbcTemplate;
import org.springframework.jdbc.core.PreparedStatementSetter;
import org.springframework.jdbc.core.RowMapper;
import org.springframework.jdbc.core.namedparam.MapSqlParameterSource;
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
              rs.getTimestamp("timeSubmitted"),
              Optional.ofNullable(rs.getTimestamp("timeCompleted")).orElse(null),
              rs.getString("status"));

  private final JdbcTemplate tspsJdbcTemplate;
  private final Logger logger = LoggerFactory.getLogger(JobsDao.class);

  @Autowired
  public JobsDao(TspsDatabaseConfiguration tspsDatabaseConfiguration) {
    this.tspsJdbcTemplate = new JdbcTemplate(tspsDatabaseConfiguration.getDataSource());
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
      tspsJdbcTemplate.update(sql, params);
      logger.info("Inserted record for job {}", jobUuid);
    } catch (DuplicateKeyException e) {
      if (e.getMessage().contains("duplicate key value violates unique constraint \"jobs_pkey\"")) {
        // Job with job_id already exists.
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
   * Return all jobs submitted by the user
   *
   * @param userId unique identifier of the calling user
   * @return List of Pipeline objects
   */
  public List<Job> getJobs(String userId, String pipelineId) {
    List<Job> jobList = new ArrayList<>();

    List<DbJob> dbJobList = getDbJobs(userId, pipelineId);

    for (DbJob dbJob : dbJobList) {
      jobList.add(Job.fromDb(dbJob));
    }

    return jobList;
  }

  private List<DbJob> getDbJobs(String userId, String pipelineId) {
    final String sql =
        """
                    SELECT job_id, user_id, pipeline_id, timeSubmitted, timeCompleted, status
                    FROM jobs
                    WHERE user_id = :userId AND pipeline_id = :pipelineId
                    """;

    MapSqlParameterSource params =
        new MapSqlParameterSource().addValue("userId", userId).addValue("pipelineId", pipelineId);

    return tspsJdbcTemplate.query(sql, (PreparedStatementSetter) params, DB_JOB_ROW_MAPPER);
  }
}
