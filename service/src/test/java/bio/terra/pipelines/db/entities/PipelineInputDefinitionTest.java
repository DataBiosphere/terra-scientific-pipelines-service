package bio.terra.pipelines.db.entities;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotEquals;

import bio.terra.pipelines.common.utils.PipelineVariableTypesEnum;
import bio.terra.pipelines.testutils.BaseTest;
import java.util.Optional;
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
            "default",
            null,
            null);
    pipelineInputDefinition.setId(5L);
    assertEquals(
        new HashCodeBuilder(17, 31)
            .append(5L)
            .append(4L)
            .append("name")
            .append("wdlVariableName")
            .append(PipelineVariableTypesEnum.FILE)
            .append("suffix")
            .append(true)
            .append(false)
            .append(false)
            .append("default")
            .append(Optional.empty())
            .append(Optional.empty())
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
            "default",
            1.0,
            2.0);
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
            "default",
            1.0,
            2.0);
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
            "default",
            3.0,
            4.0);
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
