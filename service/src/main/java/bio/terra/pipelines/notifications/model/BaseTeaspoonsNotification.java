package bio.terra.pipelines.notifications.model;

import lombok.Getter;

/**
 * Base class for all Teaspoons notifications. Contains common fields for all notification types.
 */
@Getter
public abstract class BaseTeaspoonsNotification {
  protected String notificationType;
  protected String pipelineDisplayName;
}
