package bio.terra.pipelines.service.pao;

import bio.terra.pipelines.common.exception.PolicyNotImplementedException;
import bio.terra.pipelines.common.model.PolicyInputs;
import bio.terra.pipelines.db.PaoDao;
import bio.terra.pipelines.service.pao.model.Pao;
import bio.terra.pipelines.service.pao.model.PaoComponent;
import bio.terra.pipelines.service.pao.model.PaoObjectType;
import java.util.UUID;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Component;

/**
 * The PAO Service does all processing of the Policy Attribute Objects. Those are the the objects
 * that represent policies on objects in other components of Terra. It manages the graph (DAG) of
 * objects that depend on the policies of other objects.
 */
@Component
public class PaoService {
  private static final Logger logger = LoggerFactory.getLogger(PaoService.class);

  private final PaoDao paoDao;

  @Autowired
  public PaoService(PaoDao paoDao) {
    this.paoDao = paoDao;
  }

  public void clonePao(UUID sourceObjectId, UUID destinationObjectId) {
    throw new PolicyNotImplementedException("Deprecated method. Here until we change the WSM call");
  }

  /**
   * Create a policy attribute object
   *
   * @param objectId UUID of the object - client-relevant
   * @param component identity of the component, so we know what component owns UUID
   * @param objectType type of object in the component, so the component knows where to look up the
   *     UUID
   * @param inputs policy attributes
   */
  public void createPao(
      UUID objectId, PaoComponent component, PaoObjectType objectType, PolicyInputs inputs) {
    logger.info(
        "Create PAO id {} component {} object type {}",
        objectId,
        component.name(),
        objectType.name());

    // TODO: Validate policy inputs against the policy descriptions when those are available.

    // The DAO does the heavy lifting.
    paoDao.createPao(objectId, component, objectType, inputs);
  }

  public void deletePao(UUID objectId) {
    logger.info("Delete PAO id {}", objectId);
    paoDao.markPaoDeleted(objectId);
  }

  public Pao getPao(UUID objectId) {
    logger.info("Get PAO id {}", objectId);
    return paoDao.getPao(objectId);
  }
}
