package bio.terra.pipelines.service;

import bio.terra.common.db.WriteTransaction;
import bio.terra.common.exception.InternalServerErrorException;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.db.entities.PipelineQuota;
import bio.terra.pipelines.db.entities.UserQuota;
import bio.terra.pipelines.db.repositories.PipelineQuotasRepository;
import bio.terra.pipelines.db.repositories.UserQuotasRepository;
import java.util.Optional;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Service;

/** This Service handles logic around user quotas and pipeline quotas* */
@Service
public class QuotasService {
  private static final Logger logger = LoggerFactory.getLogger(QuotasService.class);
  private final UserQuotasRepository userQuotasRepository;
  private final PipelineQuotasRepository pipelineQuotasRepository;

  @Autowired
  QuotasService(
      UserQuotasRepository userQuotasRepository,
      PipelineQuotasRepository pipelineQuotasRepository) {
    this.userQuotasRepository = userQuotasRepository;
    this.pipelineQuotasRepository = pipelineQuotasRepository;
  }

  /**
   * This method gets the quota for a given user and pipeline. If the user quota does not exist, it
   * will create a new row in the user quotas table with the default quota for the pipeline.
   */
  @WriteTransaction
  public UserQuota getQuotaForUserAndPipeline(String userId, PipelinesEnum pipelineName) {
    // try to get the user quota
    Optional<UserQuota> userQuota =
        userQuotasRepository.findByUserIdAndPipelineName(userId, pipelineName);
    // if the user quota is not found, grab the default pipeline quota and make a new row in user
    // quotas table
    if (userQuota.isEmpty()) {
      logger.debug(
          "Couldn't find user quota for user {} and pipeline {}. Creating a new row",
          userId,
          pipelineName);
      PipelineQuota pipelineQuota = pipelineQuotasRepository.findByPipelineName(pipelineName);
      UserQuota newUserQuota = new UserQuota();
      newUserQuota.setUserId(userId);
      newUserQuota.setPipelineName(pipelineName);
      newUserQuota.setQuota(pipelineQuota.getDefaultQuota());
      userQuotasRepository.save(newUserQuota);
      return newUserQuota;
    }
    return userQuota.get();
  }

  /**
   * This method updates the quota consumed for a given user and pipeline. It will return the
   * updated user quota object.
   *
   * @param userQuota - the user quota object to update
   * @param newQuotaConsumed - the quota consumed
   * @return - the updated user quota object
   */
  public UserQuota updateQuotaConsumed(UserQuota userQuota, int newQuotaConsumed) {
    if (newQuotaConsumed < 0) {
      throw new InternalServerErrorException("Quota consumed must be greater than or equal to 0");
    }
    if (newQuotaConsumed > userQuota.getQuota()) {
      throw new InternalServerErrorException("Quota consumed cannot be greater than user quota");
    }
    userQuota.setQuotaConsumed(newQuotaConsumed);
    return userQuotasRepository.save(userQuota);
  }
}
