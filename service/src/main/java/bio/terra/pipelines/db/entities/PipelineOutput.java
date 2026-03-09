package bio.terra.pipelines.db.entities;

import jakarta.persistence.*;
import lombok.Getter;
import lombok.NoArgsConstructor;
import lombok.Setter;

/**
 * JPA entity for pipeline outputs table.
 *
 * <p>Note: The uniqueConstraints and indexes defined here document the database schema but do not
 * enforce or create these constraints. The actual constraints and indexes are managed via Liquibase
 * changesets.
 */
@Entity
@Getter
@Setter
@NoArgsConstructor
@Table(
    name = "pipeline_outputs",
    uniqueConstraints = {
      @UniqueConstraint(
          name = "pipeline_run_id_output_name_uk",
          columnNames = {"pipeline_run_id", "output_name"})
    },
    indexes = {
      @Index(name = "pipeline_outputs_pipeline_run_id_idx", columnList = "pipeline_run_id")
    })
public class PipelineOutput {
  @Id
  @Column(name = "id", nullable = false)
  @GeneratedValue(strategy = GenerationType.IDENTITY)
  private Long id;

  @Column(name = "pipeline_run_id", nullable = false)
  private Long pipelineRunId;

  @Column(name = "output_name", nullable = false, columnDefinition = "TEXT")
  private String outputName;

  @Column(name = "output_value", nullable = false, columnDefinition = "TEXT")
  private String outputValue;
}
