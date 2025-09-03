package bio.terra.pipelines.db.entities;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotEquals;

import bio.terra.pipelines.common.utils.PipelineVariableTypesEnum;
import bio.terra.pipelines.testutils.BaseTest;
import org.apache.commons.lang3.builder.HashCodeBuilder;
import org.junit.jupiter.api.Test;

class PipelineInputDefinitionTest extends BaseTest {

  @Test
  void hashCodeEquals() {
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
    pipelineInputDefinition.setId(5L);
    assertEquals(
        new HashCodeBuilder(17, 31)
            .append(5L)
            .append(4L)
            .append("name")
            .append("wdlVariableName")
            .append(PipelineVariableTypesEnum.FILE)
            .append("suffix")
            .append(Boolean.TRUE)
            .append(false)
            .append(false)
            .append("default")
            .toHashCode(),
        pipelineInputDefinition.hashCode());
  }

  @Test
  void equals() {
    PipelineInputDefinition first =
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
    PipelineInputDefinition sameAsFirst =
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
    PipelineInputDefinition differentFromFirst =
        new PipelineInputDefinition(
            6L,
            "name",
            "wdlVariableName",
            PipelineVariableTypesEnum.FILE,
            "deeeeefault",
            true,
            true,
            false,
            "default");
    PipelineOutputDefinition pipelineOutputDefinition =
        new PipelineOutputDefinition(
            4L, "name", "wdlVariableName", PipelineVariableTypesEnum.FILE, true);
    assertEquals(first, first);
    assertEquals(first, sameAsFirst);
    assertNotEquals(first, differentFromFirst);
    assertNotEquals(sameAsFirst, differentFromFirst);
    assertNotEquals(first, pipelineOutputDefinition);
  }
}
