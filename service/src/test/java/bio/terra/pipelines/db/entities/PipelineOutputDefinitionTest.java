package bio.terra.pipelines.db.entities;

import static org.junit.jupiter.api.Assertions.*;

import bio.terra.pipelines.common.utils.PipelineVariableTypesEnum;
import bio.terra.pipelines.testutils.BaseTest;
import org.apache.commons.lang3.builder.HashCodeBuilder;
import org.junit.jupiter.api.Test;

class PipelineOutputDefinitionTest extends BaseTest {

  @Test
  void hashCodeEquals() {
    PipelineOutputDefinition pipelineOutputDefinition =
        new PipelineOutputDefinition(
            4L, "name", "wdlVariableName", PipelineVariableTypesEnum.FILE, true);
    pipelineOutputDefinition.setId(5L);
    assertEquals(
        new HashCodeBuilder(17, 31)
            .append(5L)
            .append(4L)
            .append("name")
            .append("wdlVariableName")
            .append(PipelineVariableTypesEnum.FILE)
            .append(true)
            .toHashCode(),
        pipelineOutputDefinition.hashCode());
  }

  @Test
  void equals() {
    PipelineOutputDefinition first =
        new PipelineOutputDefinition(
            4L, "name", "wdlVariableName", PipelineVariableTypesEnum.FILE, true);
    PipelineOutputDefinition sameAsFirst =
        new PipelineOutputDefinition(
            4L, "name", "wdlVariableName", PipelineVariableTypesEnum.FILE, true);
    PipelineOutputDefinition differentFromFirst =
        new PipelineOutputDefinition(6L, "name", "agadg", PipelineVariableTypesEnum.FILE, true);
    PipelineInputDefinition pipelineInputDefinition =
        new PipelineInputDefinition(
            4L,
            "name",
            "wdlVariableName",
            PipelineVariableTypesEnum.FILE,
            "suffix",
            true,
            false,
            false,
            "default");

    assertEquals(first, first);
    assertEquals(first, sameAsFirst);
    assertNotEquals(first, differentFromFirst);
    assertNotEquals(sameAsFirst, differentFromFirst);
    assertNotEquals(first, pipelineInputDefinition);
  }
}
