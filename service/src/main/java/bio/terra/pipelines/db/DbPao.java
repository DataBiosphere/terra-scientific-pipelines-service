package bio.terra.pipelines.db;

import bio.terra.pipelines.service.pao.model.PaoComponent;
import bio.terra.pipelines.service.pao.model.PaoObjectType;
import java.util.Set;
import java.util.UUID;

/** Record to hold a PAO record when processing in the PaoDao */
public record DbPao(
    UUID objectId,
    PaoComponent component,
    PaoObjectType objectType,
    Set<String> sources,
    String attributeSetId,
    String effectiveSetId,
    boolean deleted) {}