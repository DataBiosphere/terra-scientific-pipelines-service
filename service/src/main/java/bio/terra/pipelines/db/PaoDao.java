package bio.terra.pipelines.db;

import bio.terra.common.db.ReadTransaction;
import bio.terra.common.db.WriteTransaction;
import bio.terra.pipelines.app.configuration.TspsDatabaseConfiguration;
import bio.terra.pipelines.common.exception.PolicyObjectNotFoundException;
import bio.terra.pipelines.common.model.PolicyInput;
import bio.terra.pipelines.common.model.PolicyInputs;
import bio.terra.pipelines.common.model.PolicyName;
import bio.terra.pipelines.db.exception.DuplicateObjectException;
import bio.terra.pipelines.service.pao.model.Pao;
import bio.terra.pipelines.service.pao.model.PaoComponent;
import bio.terra.pipelines.service.pao.model.PaoObjectType;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.UUID;
import java.util.stream.Collectors;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.dao.DuplicateKeyException;
import org.springframework.jdbc.core.RowMapper;
import org.springframework.jdbc.core.namedparam.MapSqlParameterSource;
import org.springframework.jdbc.core.namedparam.NamedParameterJdbcTemplate;
import org.springframework.stereotype.Component;

@Component
public class PaoDao {
  private static final RowMapper<DbPao> DB_PAO_ROW_MAPPER =
      (rs, rowNum) -> {
        String[] sourcesArray = (String[]) rs.getArray("sources").getArray();
        Set<String> sources = new HashSet<>(Arrays.asList(sourcesArray));
        return new DbPao(
            UUID.fromString(rs.getString("object_id")),
            PaoComponent.fromDb(rs.getString("component")),
            PaoObjectType.fromDb(rs.getString("object_type")),
            sources,
            rs.getString("attribute_set_id"),
            rs.getString("effective_set_id"),
            rs.getBoolean("deleted"));
      };

  private static final RowMapper<DbAttribute> DB_ATTRIBUTE_SET_ROW_MAPPER =
      (rs, rowNum) -> {
        String[] conflictsArray = (String[]) rs.getArray("conflicts").getArray();
        Set<UUID> conflicts =
            Arrays.stream(conflictsArray).map(UUID::fromString).collect(Collectors.toSet());

        return new DbAttribute(
            rs.getString("set_id"),
            new PolicyInput(
                new PolicyName(rs.getString("namespace"), rs.getString("name")),
                DbAdditionalData.fromDb(rs.getString("properties")),
                conflicts));
      };

  private final Logger logger = LoggerFactory.getLogger(PaoDao.class);
  private final NamedParameterJdbcTemplate tpsJdbcTemplate;

  @Autowired
  public PaoDao(TspsDatabaseConfiguration tspsDatabaseConfiguration) {
    this.tpsJdbcTemplate =
        new NamedParameterJdbcTemplate(tspsDatabaseConfiguration.getDataSource());
  }

  @WriteTransaction
  public void createPao(
      UUID objectId, PaoComponent component, PaoObjectType objectType, PolicyInputs inputs) {

    // Store the attribute set twice: once as the object's set and once as its effective set.
    // We could optimize this case, but the logic is cleaner if we treat them distinctly from the
    // outset.
    String attributeSetId = UUID.randomUUID().toString();
    String effectiveSetId = UUID.randomUUID().toString();
    createAttributeSet(attributeSetId, inputs);
    createAttributeSet(effectiveSetId, inputs);

    createDbPao(
        objectId,
        component.getDbComponent(),
        objectType.getDbObjectType(),
        "",
        attributeSetId,
        effectiveSetId);
  }

  @WriteTransaction
  public void deletePaos(Collection<DbPao> paos) {
    paos.forEach((DbPao pao) -> removeDbPao(pao));
  }

  /**
   * Set the 'deleted' field on the PAO to true.
   *
   * @param objectId The PAO to flag
   */
  public void markPaoDeleted(UUID objectId) {
    final String sql = "UPDATE policy_object SET deleted=true WHERE object_id=:object_id";
    MapSqlParameterSource params =
        new MapSqlParameterSource().addValue("object_id", objectId.toString());
    tpsJdbcTemplate.update(sql, params);
  }

  @ReadTransaction
  public Pao getPao(UUID objectId) {
    DbPao dbPao = getDbPao(objectId);
    Map<String, PolicyInputs> attributeSetMap =
        getAttributeSets(List.of(dbPao.attributeSetId(), dbPao.effectiveSetId()));
    return Pao.fromDb(dbPao, attributeSetMap);
  }

  // -- Graph Walk Methods --
  // These methods are intentionally without transaction annotations. They are used by the policy
  // update process. That process may do multiple reads of the database followed by a big update.
  // The transaction is managed outside of the DAO.

  /**
   * Given a list of PAO ids, return PAOs
   *
   * @param objectIdList UUIDs of Policy Attribute Objects
   * @return List of Pao objects
   */
  public List<Pao> getPaos(List<UUID> objectIdList) {
    List<Pao> paoList = new ArrayList<>();
    if (objectIdList.isEmpty()) {
      // Nothing to do
      return paoList;
    }

    List<DbPao> dbPaoList = getDbPaos(objectIdList);
    List<String> setIdList = new ArrayList<>();
    for (DbPao dbPao : dbPaoList) {
      setIdList.add(dbPao.attributeSetId());
      setIdList.add(dbPao.effectiveSetId());
    }

    // Gather all of the attribute sets
    Map<String, PolicyInputs> attributeSetMap = getAttributeSets(setIdList);

    for (DbPao dbPao : dbPaoList) {
      paoList.add(Pao.fromDb(dbPao, attributeSetMap));
    }

    return paoList;
  }

  /**
   * Given a source id, find the immediate dependents and return their ids
   *
   * @param sourceId source to hunt for
   * @return Set of dependent UUIDs that reference that source id
   */
  public Set<UUID> getDependentIds(UUID sourceId) {
    final String sql = "SELECT object_id FROM policy_object WHERE :source_id = ANY(sources)";

    MapSqlParameterSource params =
        new MapSqlParameterSource().addValue("source_id", sourceId.toString());

    return new HashSet<>(
        tpsJdbcTemplate.query(
            sql, params, (rs, rowNum) -> UUID.fromString(rs.getString("object_id"))));
  }

  /**
   * Given a source id, RECURSIVELY find all of the dependents and return their ids. Consult the
   * documentation on Postgres recursive queries in this document:
   *
   * <p>https://www.postgresql.org/docs/current/queries-with.html
   *
   * <p>And in this tutorial:
   *
   * <p>https://www.postgresqltutorial.com/postgresql-tutorial/postgresql-recursive-query/
   *
   * @param sourceId source to hunt for
   * @return Set of all dependent UUIDs
   */
  public Set<UUID> getAllDependentIds(UUID sourceId) {
    final String sql =
        """
        WITH RECURSIVE dependents AS (
          SELECT :source_id AS object_id
          UNION
          SELECT R.object_id FROM policy_object R INNER JOIN dependents D ON D.object_id = ANY(R.sources))
        SELECT * FROM dependents;
        """;

    MapSqlParameterSource params =
        new MapSqlParameterSource().addValue("source_id", sourceId.toString());

    return new HashSet<>(
        tpsJdbcTemplate.query(
            sql, params, (rs, rowNum) -> UUID.fromString(rs.getString("object_id"))));
  }

  private void removeDbPao(DbPao dbPao) {
    try {
      // Delete associated attribute set(s)
      deleteAttributeSet(dbPao.attributeSetId());
      deleteAttributeSet(dbPao.effectiveSetId());

      // Delete the policy object
      final String sql = "DELETE FROM policy_object WHERE object_id=:object_id";
      MapSqlParameterSource params =
          new MapSqlParameterSource().addValue("object_id", dbPao.objectId().toString());
      tpsJdbcTemplate.update(sql, params);
    } catch (PolicyObjectNotFoundException e) {
      // Delete throws no error on not found
    }
  }

  /**
   * The JdbcTemplate doesn't have a nice way to pass arrays into and out of Postgres. Since we use
   * UUIDs, we know there are no commas, so we can build the CSV and use the Postgres
   * string_to_array builtin to set the array value.
   *
   * @param uuidSet input set of UUIDs
   * @return csv of strings
   */
  private String makeCsvFromUuidSet(Set<UUID> uuidSet) {
    return String.join(",", uuidSet.stream().map(UUID::toString).toList());
  }

  private void createAttributeSet(String setId, PolicyInputs inputs) {
    final String setsql =
        """
        INSERT INTO attribute_set(set_id, namespace, name, properties, conflicts)
        VALUES(:set_id, :namespace, :name, cast(:properties AS jsonb), string_to_array(:conflicts,','))
        """;

    for (PolicyInput input : inputs.getInputs().values()) {
      String conflictCsv = makeCsvFromUuidSet(input.getConflicts());
      MapSqlParameterSource setparams =
          new MapSqlParameterSource()
              .addValue("set_id", setId)
              .addValue("namespace", input.getPolicyName().getNamespace())
              .addValue("name", input.getPolicyName().getName())
              .addValue("properties", DbAdditionalData.toDb(input.getAdditionalData()))
              .addValue("conflicts", conflictCsv);
      tpsJdbcTemplate.update(setsql, setparams);
      logger.info(
          "Inserted record for pao set id {}, policy {}, conflicts {}",
          setId,
          input.getPolicyName(),
          conflictCsv);
    }
  }

  private void createDbPao(
      UUID objectId,
      String component,
      String objectType,
      String sources,
      String attributeSetId,
      String effectiveSetId) {
    final String sql =
        """
        INSERT INTO policy_object
          (object_id, component, object_type, sources, attribute_set_id, effective_set_id)
        VALUES
          (:object_id, :component, :object_type, string_to_array(:sources, ','), :attribute_set_id, :effective_set_id)
        """;

    MapSqlParameterSource params =
        new MapSqlParameterSource()
            .addValue("object_id", objectId)
            .addValue("component", component)
            .addValue("object_type", objectType)
            .addValue("sources", sources)
            .addValue("attribute_set_id", attributeSetId)
            .addValue("effective_set_id", effectiveSetId);

    try {
      tpsJdbcTemplate.update(sql, params);
      logger.info("Inserted record for pao {}", objectId);
    } catch (DuplicateKeyException e) {
      throw new DuplicateObjectException(
          "Duplicate policy attributes object with objectId " + objectId);
    }
  }

  private void deleteAttributeSet(String setId) {
    final String sql = "DELETE FROM attribute_set WHERE set_id = :set_id";
    final var params = new MapSqlParameterSource().addValue("set_id", setId);
    tpsJdbcTemplate.update(sql, params);
  }

  public DbPao getDbPao(UUID objectId) {
    final String sql =
        """
        SELECT object_id, component, object_type, attribute_set_id, effective_set_id, sources, deleted
        FROM policy_object WHERE object_id = :object_id
        """;

    MapSqlParameterSource params =
        new MapSqlParameterSource().addValue("object_id", objectId.toString());

    List<DbPao> dbPao = tpsJdbcTemplate.query(sql, params, DB_PAO_ROW_MAPPER);
    if (dbPao.isEmpty()) {
      throw new PolicyObjectNotFoundException("Policy object not found: " + objectId);
    }
    return dbPao.get(0);
  }

  private List<DbPao> getDbPaos(List<UUID> objectIdList) {
    final String sql =
        """
        SELECT object_id, component, object_type, attribute_set_id, effective_set_id, sources, deleted
        FROM policy_object
        WHERE object_id IN (:object_id_list)
        """;

    if (objectIdList.isEmpty()) {
      return new ArrayList<>();
    }
    List<String> stringIdList = objectIdList.stream().map(UUID::toString).toList();

    MapSqlParameterSource params =
        new MapSqlParameterSource().addValue("object_id_list", stringIdList);

    return tpsJdbcTemplate.query(sql, params, DB_PAO_ROW_MAPPER);
  }

  /**
   * Build attribute sets from a list of set ids, possibly duplicate ones.
   *
   * @param setIdList list of set ids
   * @return map of set id to attribute set
   */
  private Map<String, PolicyInputs> getAttributeSets(List<String> setIdList) {
    final String sql =
        """
        SELECT set_id, namespace, name, properties, conflicts
        FROM attribute_set
        WHERE set_id IN (:set_id_list)
        """;

    var attributeSets = new HashMap<String, PolicyInputs>();
    if (setIdList.isEmpty()) {
      // Nothing to do - skip the query
      return attributeSets;
    }

    var uniqueSetIds = new HashSet<>(setIdList);
    var params = new MapSqlParameterSource().addValue("set_id_list", uniqueSetIds);

    List<DbAttribute> attributeList =
        tpsJdbcTemplate.query(sql, params, DB_ATTRIBUTE_SET_ROW_MAPPER);

    // Initialize the attributeSets with all input set ids and an empty PolicyInputs.
    // That is a valid return and covers the case where attribute sets are empty,
    // so do not have any rows in the attribute set table.
    for (String id : uniqueSetIds) {
      attributeSets.put(id, new PolicyInputs());
    }

    // For attribute rows we have, find their attribute set and add the policy
    for (DbAttribute attribute : attributeList) {
      PolicyInputs policyInputs = attributeSets.get(attribute.setId());
      policyInputs.addInput(attribute.policyInput());
    }

    return attributeSets;
  }
}