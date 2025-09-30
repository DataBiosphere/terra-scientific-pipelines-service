package bio.terra.pipelines.common.utils;

import static bio.terra.pipelines.common.utils.PipelineVariableTypesEnum.BOOLEAN;
import static bio.terra.pipelines.common.utils.PipelineVariableTypesEnum.FILE;
import static bio.terra.pipelines.common.utils.PipelineVariableTypesEnum.FILE_ARRAY;
import static bio.terra.pipelines.common.utils.PipelineVariableTypesEnum.INTEGER;
import static bio.terra.pipelines.common.utils.PipelineVariableTypesEnum.STRING;
import static bio.terra.pipelines.common.utils.PipelineVariableTypesEnum.STRING_ARRAY;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.params.provider.Arguments.arguments;

import bio.terra.pipelines.db.entities.PipelineInputDefinition;
import bio.terra.pipelines.testutils.BaseTest;
import com.fasterxml.jackson.core.type.TypeReference;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.stream.Stream;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.Arguments;
import org.junit.jupiter.params.provider.MethodSource;

class PipelineVariableTypesEnumTest extends BaseTest {

  private static Stream<Arguments> castValidations() {
    String commonInputName = "inputName";
    // error messages
    String stringTypeErrorMessage = "%s must be a string".formatted(commonInputName);
    String stringPatternErrorMessage =
        "%s must only contain alphanumeric characters, dashes, underscores, and periods"
            .formatted(commonInputName);
    String integerTypeErrorMessage = "%s must be an integer".formatted(commonInputName);
    String booleanTypeErrorMessage = "%s must be a boolean".formatted(commonInputName);
    String vcfTypeErrorMessage =
        "%s must be a path to a file ending in .vcf.gz".formatted(commonInputName);
    String stringArrayTypeErrorMessage =
        "%s must be an array of strings".formatted(commonInputName);
    String vcfArrayTypeErrorMessage =
        "%s must be an array of paths to files ending in .vcf.gz".formatted(commonInputName);
    String notNullOrEmptyErrorMessage = "%s must not be null or empty".formatted(commonInputName);

    // the only used information in the input definitions is name, type, and fileSuffix
    PipelineInputDefinition integerInputDefinition =
        new PipelineInputDefinition(
            1L, commonInputName, "integer_input", INTEGER, null, true, true, false, null);
    PipelineInputDefinition stringInputDefinition =
        new PipelineInputDefinition(
            1L, commonInputName, "string_input", STRING, null, true, true, false, null);
    PipelineInputDefinition booleanInputDefinition =
        new PipelineInputDefinition(
            1L, commonInputName, "boolean_input", BOOLEAN, null, true, true, false, null);
    PipelineInputDefinition fileVcfInputDefinition =
        new PipelineInputDefinition(
            1L, commonInputName, "file_vcf_input", FILE, ".vcf.gz", true, true, false, null);
    PipelineInputDefinition fileBedInputDefinition =
        new PipelineInputDefinition(
            1L, commonInputName, "file_bed_input", FILE, ".bed", true, true, false, null);
    PipelineInputDefinition fileNoSuffixInputDefinition =
        new PipelineInputDefinition(
            1L, commonInputName, "file_no_suffix_input", FILE, "", true, true, false, null);
    PipelineInputDefinition stringArrayInputDefinition =
        new PipelineInputDefinition(
            1L, commonInputName, "string_array_input", STRING_ARRAY, null, true, true, false, null);
    PipelineInputDefinition fileArrayVcfInputDefinition =
        new PipelineInputDefinition(
            1L,
            commonInputName,
            "file_array_vcf_input",
            FILE_ARRAY,
            ".vcf.gz",
            true,
            true,
            false,
            null);
    PipelineInputDefinition fileArrayBedInputDefinition =
        new PipelineInputDefinition(
            1L,
            commonInputName,
            "file_array_bed_input",
            FILE_ARRAY,
            ".bed",
            true,
            true,
            false,
            null);

    return Stream.of(
        // arguments: input definition, input value to cast, expected cast value if
        // successful, error message if fails

        // INTEGER
        arguments(integerInputDefinition, 123, 123, null),
        arguments(integerInputDefinition, "123", 123, null),
        arguments(integerInputDefinition, "I am a string", null, integerTypeErrorMessage),
        arguments(integerInputDefinition, 2.3, null, integerTypeErrorMessage),
        arguments(integerInputDefinition, "2.3", null, integerTypeErrorMessage),
        arguments(integerInputDefinition, null, null, integerTypeErrorMessage),
        arguments(integerInputDefinition, "", null, integerTypeErrorMessage),

        // STRING
        arguments(stringInputDefinition, "IAmAString", "IAmAString", null),
        arguments(
            stringInputDefinition, "I am a string", "I am a string", stringPatternErrorMessage),
        arguments(
            stringInputDefinition, "\"IAmAString\"", "\"IAmAString\"", stringPatternErrorMessage),
        arguments(
            stringInputDefinition,
            "$tr1nG.w1th.0th3r.ch@r@ct3r$",
            "$tr1nG.w1th.0th3r.ch@r@ct3r$",
            stringPatternErrorMessage),
        arguments(
            stringInputDefinition,
            "    IAmAStringWithExtraWhitespace    ",
            "IAmAStringWithExtraWhitespace",
            null),
        arguments(
            stringInputDefinition,
            List.of("this", "is", "not", "a", "string"),
            null,
            stringTypeErrorMessage),
        arguments(stringInputDefinition, 123, null, stringTypeErrorMessage),
        arguments(stringInputDefinition, null, null, stringTypeErrorMessage),
        arguments(stringInputDefinition, "", null, stringTypeErrorMessage),

        // BOOLEAN
        arguments(booleanInputDefinition, true, true, null),
        arguments(booleanInputDefinition, false, false, null),
        arguments(booleanInputDefinition, "true", true, null),
        arguments(booleanInputDefinition, " false ", false, null),
        arguments(booleanInputDefinition, "TRUE", true, null),
        arguments(booleanInputDefinition, "FALSE", false, null),
        arguments(booleanInputDefinition, "foo", null, booleanTypeErrorMessage),
        arguments(booleanInputDefinition, 1, null, booleanTypeErrorMessage),
        arguments(booleanInputDefinition, null, null, booleanTypeErrorMessage),

        // FILE
        arguments(fileVcfInputDefinition, "path/to/file.vcf.gz", "path/to/file.vcf.gz", null),
        arguments(fileVcfInputDefinition, "   path/to/file.vcf.gz   ", "path/to/file.vcf.gz", null),
        arguments(
            fileVcfInputDefinition,
            "path/to/file.vcf",
            "path/to/file.vcf",
            vcfTypeErrorMessage), // cast is successful but validation fails
        arguments(fileVcfInputDefinition, 3, null, vcfTypeErrorMessage),
        arguments(fileVcfInputDefinition, null, null, vcfTypeErrorMessage),
        arguments(fileVcfInputDefinition, "", null, vcfTypeErrorMessage),
        arguments(
            fileBedInputDefinition,
            "path/to/not/a/bed/file",
            "path/to/not/a/bed/file",
            "%s must be a path to a file ending in .bed"
                .formatted(commonInputName)), // cast is successful but validation fails
        arguments(fileNoSuffixInputDefinition, "path/to/file", "path/to/file", null),

        // STRING_ARRAY
        arguments(
            stringArrayInputDefinition,
            List.of("this", "is", "a", "list", "of", "strings"),
            List.of("this", "is", "a", "list", "of", "strings"),
            null),
        arguments(
            stringArrayInputDefinition, List.of("singletonList"), List.of("singletonList"), null),
        arguments(
            stringArrayInputDefinition,
            List.of("this ", " is", " a ", "  list", "of  ", "  strings  ", "withExtraWhitespace"),
            List.of("this", "is", "a", "list", "of", "strings", "withExtraWhitespace"),
            null),
        arguments(
            stringArrayInputDefinition,
            "[\"this\", \"is\", \"a\", \"stringy\", \"list\",  \"of\", \"strings\"]",
            List.of("this", "is", "a", "stringy", "list", "of", "strings"),
            null),
        arguments(
            stringArrayInputDefinition,
            "[\" stringy \", \" and\",  \"whitespace \"]",
            List.of("stringy", "and", "whitespace"),
            null),
        arguments(
            stringArrayInputDefinition, "I am not an array", null, stringArrayTypeErrorMessage),
        arguments(
            stringArrayInputDefinition,
            Arrays.asList("string", 2, 3),
            null,
            stringArrayTypeErrorMessage),
        arguments(
            stringArrayInputDefinition, Arrays.asList(1, 2, 3), null, stringArrayTypeErrorMessage),
        arguments(stringArrayInputDefinition, null, null, notNullOrEmptyErrorMessage),
        arguments(
            stringArrayInputDefinition,
            Collections.emptyList(),
            Collections.emptyList(), // cast is successful but validation fails
            notNullOrEmptyErrorMessage),
        arguments(
            stringArrayInputDefinition,
            Arrays.asList("", "array with empty string"),
            null,
            stringArrayTypeErrorMessage),
        arguments(
            stringArrayInputDefinition,
            Arrays.asList(null, "array with null"),
            null,
            stringArrayTypeErrorMessage),

        // FILE_ARRAY
        arguments(
            fileArrayVcfInputDefinition,
            List.of("path/to/file.vcf.gz", "path/to/file2.vcf.gz"),
            List.of("path/to/file.vcf.gz", "path/to/file2.vcf.gz"),
            null),
        arguments(
            fileArrayVcfInputDefinition,
            List.of("  path/to/file/with/extra/whitespace.vcf.gz  "),
            List.of("path/to/file/with/extra/whitespace.vcf.gz"),
            null),
        arguments(
            fileArrayVcfInputDefinition,
            "[\"path/to/file1.vcf.gz\", \"path/to/file2.vcf.gz\"]",
            List.of("path/to/file1.vcf.gz", "path/to/file2.vcf.gz"),
            null),
        arguments(
            fileArrayVcfInputDefinition,
            "[\" path/to/file1.vcf.gz \", \"path/to/file2.vcf.gz \"]",
            List.of("path/to/file1.vcf.gz", "path/to/file2.vcf.gz"),
            null),
        arguments(
            fileArrayVcfInputDefinition,
            "this/is/not/an/array.vcf.gz",
            null,
            vcfArrayTypeErrorMessage),
        arguments(
            fileArrayVcfInputDefinition,
            Arrays.asList("path/to/file.vcf.gz", "just/a/string"),
            Arrays.asList(
                "path/to/file.vcf.gz", "just/a/string"), // cast is successful but validation fails
            vcfArrayTypeErrorMessage),
        arguments(
            fileArrayVcfInputDefinition,
            Arrays.asList("path/to/file.vcf.gz", 2.5),
            null,
            vcfArrayTypeErrorMessage),
        arguments(
            fileArrayBedInputDefinition,
            Arrays.asList("file.vcf.gz", "file.bed"),
            Arrays.asList("file.vcf.gz", "file.bed"),
            "%s must be an array of paths to files ending in .bed"
                .formatted(commonInputName)), // cast is successful but validation fails
        arguments(fileArrayVcfInputDefinition, null, null, notNullOrEmptyErrorMessage),
        arguments(
            fileArrayVcfInputDefinition,
            Collections.emptyList(),
            Collections.emptyList(),
            notNullOrEmptyErrorMessage), // cast is successful but validation fails
        arguments(
            fileArrayVcfInputDefinition,
            Arrays.asList(null, "list/with/null.vcf.gz"),
            null,
            vcfArrayTypeErrorMessage));
  }

  @ParameterizedTest
  @MethodSource("castValidations")
  <T> void castInputValues(
      PipelineInputDefinition inputDefinition, Object inputValue, T expectedCastValue) {
    if (expectedCastValue != null) {
      if (expectedCastValue instanceof List expectedListCastValue) {
        List<String> listCastValue =
            inputDefinition.getType().cast("fieldName", inputValue, new TypeReference<>() {});

        assertEquals(expectedCastValue, listCastValue);
        // Ensure that base class matches up
        if (!expectedListCastValue.isEmpty()) {
          assertEquals(expectedListCastValue.get(0).getClass(), listCastValue.get(0).getClass());
        }
      } else {
        assertEquals(
            expectedCastValue,
            inputDefinition.getType().cast("fieldName", inputValue, new TypeReference<>() {}));
        // Ensure that class matches up
        assertEquals(
            expectedCastValue.getClass(),
            inputDefinition
                .getType()
                .cast("fieldName", inputValue, new TypeReference<>() {})
                .getClass());
      }
    } else {
      // cast should return null
      assertNull(inputDefinition.getType().cast("fieldName", inputValue, new TypeReference<>() {}));
    }
  }

  @ParameterizedTest
  @MethodSource("castValidations")
  <T> void validateInputValues(
      PipelineInputDefinition inputDefinition,
      Object inputValue,
      T expectedCastValue,
      String expectedErrorMessage) {
    if (expectedErrorMessage == null) {
      // should validate
      assertNull(inputDefinition.getType().validate(inputDefinition, inputValue));
    } else {
      // should not validate
      assertEquals(
          expectedErrorMessage, inputDefinition.getType().validate(inputDefinition, inputValue));
    }
  }
}
