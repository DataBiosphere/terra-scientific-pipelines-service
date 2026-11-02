package bio.terra.pipelines.db.entities;

import jakarta.persistence.Column;
import jakarta.persistence.Entity;
import jakarta.persistence.Id;
import jakarta.persistence.Table;
import java.time.Instant;
import lombok.Getter;
import lombok.NoArgsConstructor;
import lombok.Setter;
import org.hibernate.annotations.SourceType;
import org.hibernate.annotations.UpdateTimestamp;

/** Single-row entity (id is always 1) holding service-wide state. */
@Entity
@Getter
@Setter
@NoArgsConstructor
@Table(name = "service_state")
public class ServiceState {
  public static final int SINGLETON_ID = 1;

  @Id
  @Column(name = "id", nullable = false)
  private Integer id = SINGLETON_ID;

  @Column(name = "paused", nullable = false)
  private boolean paused;

  @Column(name = "updated", insertable = false)
  @UpdateTimestamp(source = SourceType.DB)
  private Instant updated;
}
