package bio.terra.pipelines.db.entities;

import jakarta.persistence.Column;
import jakarta.persistence.Entity;
import jakarta.persistence.GeneratedValue;
import jakarta.persistence.GenerationType;
import jakarta.persistence.Id;
import jakarta.persistence.Table;
import java.time.Instant;
import java.util.UUID;
import lombok.Getter;
import lombok.NoArgsConstructor;
import lombok.Setter;
import org.hibernate.annotations.CreationTimestamp;
import org.hibernate.annotations.SourceType;
import org.hibernate.annotations.UpdateTimestamp;

@Entity
@Getter
@Setter
@NoArgsConstructor
@Table(name = "download_output_calls")
public class DownloadCallCount {
  @Id
  @Column(name = "id", nullable = false)
  @GeneratedValue(strategy = GenerationType.IDENTITY)
  private Long id;

  @Column(name = "job_id", nullable = false, unique = true)
  private UUID jobId;

  @Column(name = "count")
  private int count;

  @Column(name = "first_call", insertable = false)
  @CreationTimestamp(source = SourceType.DB)
  private Instant firstCall;

  @Column(name = "latest_call", insertable = false)
  @UpdateTimestamp(source = SourceType.DB)
  private Instant latestCall;

  public DownloadCallCount(UUID jobId, int count) {
    this.jobId = jobId;
    this.count = count;
  }
}
