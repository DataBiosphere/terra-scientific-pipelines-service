package bio.terra.pipelines.common.utils;

import static bio.terra.pipelines.common.utils.PipelineVariableTypesEnum.*;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.junit.jupiter.params.provider.Arguments.arguments;

import bio.terra.pipelines.db.entities.PipelineInputDefinition;
import bio.terra.pipelines.testutils.BaseTest;
import com.fasterxml.jackson.core.type.TypeReference;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.stream.Stream;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.Arguments;
import org.junit.jupiter.params.provider.MethodSource;

class PipelineVariableTypesEnumTest extends BaseTest {

  private static Stream<Arguments> castValidations() {
    String commonInputName = "inputName";
    Double commonMinValue = 0.0;
    Double commonMaxValue = 200.0;
    Double commonFloatMinValue = 0.0;
    Double commonFloatMaxValue = 3.0;
    // error messages
    String stringTypeErrorMessage = "%s must be a string".formatted(commonInputName);
    String integerTypeErrorMessage = "%s must be an integer".formatted(commonInputName);
    String typeErrorMessageRange =
        "%s must be between %s and %s".formatted(commonInputName, commonMinValue, commonMaxValue);
    String typeErrorMessageMin =
        "%s must be at least %s".formatted(commonInputName, commonMinValue);
    String typeErrorMessageMax = "%s must be at most %s".formatted(commonInputName, commonMaxValue);
    String floatTypeErrorMessage = "%s must be a float".formatted(commonInputName);
    String floatTypeErrorMessageRange =
        "%s must be between %s and %s"
            .formatted(commonInputName, commonFloatMinValue, commonFloatMaxValue);
    String floatTypeErrorMessageMin =
        "%s must be at least %s".formatted(commonInputName, commonFloatMinValue);
    String floatTypeErrorMessageMax =
        "%s must be at most %s".formatted(commonInputName, commonFloatMaxValue);
    String booleanTypeErrorMessage = "%s must be a boolean".formatted(commonInputName);
    String vcfTypeErrorMessage =
        "%s must be a path to a file ending in .vcf.gz".formatted(commonInputName);
    String stringArrayTypeErrorMessage =
        "%s must be an array of strings".formatted(commonInputName);
    String tsvTypeErrorMessage =
        "%s must be a path to a file ending in .tsv".formatted(commonInputName);
    String vcfArrayTypeErrorMessage =
        "%s must be an array of paths to files ending in .vcf.gz and containing only alphanumeric characters or the following symbols: -_.=:\\/" // contains colon
            .formatted(commonInputName);
    String notNullOrEmptyErrorMessage = "%s must not be null or empty".formatted(commonInputName);
    String stringPatternErrorMessage =
        "%s must only contain alphanumeric characters or the following symbols: -_.=\\/"
            .formatted(commonInputName);
    String stringArrayPatternErrorMessage =
        "%s must only contain strings with alphanumeric characters or the following symbols: -_.=\\/"
            .formatted(commonInputName);
    String filePatternErrorMessage =
        "%s must only contain alphanumeric characters or the following symbols: -_.=:\\/" // contains colon
            .formatted(commonInputName);

    // the only used information in the input definitions is name, type, fileSuffix, minValue, and
    // maxValue
    PipelineInputDefinition integerInputDefinition =
        new PipelineInputDefinition(
            1L,
            commonInputName,
            "integer_input",
            null,
            null,
            INTEGER,
            null,
            true,
            true,
            false,
            null,
            null,
            null);
    PipelineInputDefinition integerInputDefinitionRange =
        new PipelineInputDefinition(
            1L,
            commonInputName,
            "integer_input",
            null,
            null,
            INTEGER,
            null,
            true,
            true,
            false,
            null,
            commonMinValue,
            commonMaxValue);
    PipelineInputDefinition integerInputDefinitionOnlyMin =
        new PipelineInputDefinition(
            1L,
            commonInputName,
            "integer_input",
            null,
            null,
            INTEGER,
            null,
            true,
            true,
            false,
            null,
            commonMinValue,
            null);
    PipelineInputDefinition integerInputDefinitionOnlyMax =
        new PipelineInputDefinition(
            1L,
            commonInputName,
            "integer_input",
            null,
            null,
            INTEGER,
            null,
            true,
            true,
            false,
            null,
            null,
            commonMaxValue);
    PipelineInputDefinition floatInputDefinition =
        new PipelineInputDefinition(
            1L,
            commonInputName,
            "float_input",
            null,
            null,
            FLOAT,
            null,
            true,
            true,
            false,
            null,
            null,
            null);
    PipelineInputDefinition floatInputDefinitionRange =
        new PipelineInputDefinition(
            1L,
            commonInputName,
            "float_input",
            null,
            null,
            FLOAT,
            null,
            true,
            true,
            false,
            null,
            commonFloatMinValue,
            commonFloatMaxValue);
    PipelineInputDefinition floatInputDefinitionOnlyMin =
        new PipelineInputDefinition(
            1L,
            commonInputName,
            "float_input",
            null,
            null,
            FLOAT,
            null,
            true,
            true,
            false,
            null,
            commonFloatMinValue,
            null);
    PipelineInputDefinition floatInputDefinitionOnlyMax =
        new PipelineInputDefinition(
            1L,
            commonInputName,
            "float_input",
            null,
            null,
            FLOAT,
            null,
            true,
            true,
            false,
            null,
            null,
            commonFloatMaxValue);
    PipelineInputDefinition stringInputDefinition =
        new PipelineInputDefinition(
            1L,
            commonInputName,
            "string_input",
            null,
            null,
            STRING,
            null,
            true,
            true,
            false,
            null,
            null,
            null);
    PipelineInputDefinition booleanInputDefinition =
        new PipelineInputDefinition(
            1L,
            commonInputName,
            "boolean_input",
            null,
            null,
            BOOLEAN,
            null,
            true,
            true,
            false,
            null,
            null,
            null);
    PipelineInputDefinition fileVcfInputDefinition =
        new PipelineInputDefinition(
            1L,
            commonInputName,
            "file_vcf_input",
            null,
            null,
            FILE,
            ".vcf.gz",
            true,
            true,
            false,
            null,
            null,
            null);
    PipelineInputDefinition fileBedInputDefinition =
        new PipelineInputDefinition(
            1L,
            commonInputName,
            "file_bed_input",
            null,
            null,
            FILE,
            ".bed",
            true,
            true,
            false,
            null,
            null,
            null);
    PipelineInputDefinition fileNoSuffixInputDefinition =
        new PipelineInputDefinition(
            1L,
            commonInputName,
            "file_no_suffix_input",
            null,
            null,
            FILE,
            "",
            true,
            true,
            false,
            null,
            null,
            null);
    PipelineInputDefinition manifestTsvInputDefinition =
        new PipelineInputDefinition(
            1L,
            commonInputName,
            "manifest_input",
            null,
            null,
            MANIFEST,
            ".tsv",
            true,
            true,
            false,
            null,
            null,
            null);
    PipelineInputDefinition stringArrayInputDefinition =
        new PipelineInputDefinition(
            1L,
            commonInputName,
            "string_array_input",
            null,
            null,
            STRING_ARRAY,
            null,
            true,
            true,
            false,
            null,
            null,
            null);
    PipelineInputDefinition fileArrayVcfInputDefinition =
        new PipelineInputDefinition(
            1L,
            commonInputName,
            "file_array_vcf_input",
            null,
            null,
            FILE_ARRAY,
            ".vcf.gz",
            true,
            true,
            false,
            null,
            null,
            null);
    PipelineInputDefinition fileArrayBedInputDefinition =
        new PipelineInputDefinition(
            1L,
            commonInputName,
            "file_array_bed_input",
            null,
            null,
            FILE_ARRAY,
            ".bed",
            true,
            true,
            false,
            null,
            null,
            null);

    return Stream.of(
        // arguments: input definition, input value to cast, expected cast value if
        // successful, error message if cast or validation fails

        // INTEGER
        arguments(integerInputDefinition, 123, 123, null),
        arguments(integerInputDefinition, "123", 123, null),
        arguments(integerInputDefinition, "I am a string", null, integerTypeErrorMessage),
        arguments(integerInputDefinition, 2.3, null, integerTypeErrorMessage),
        arguments(integerInputDefinition, "2.3", null, integerTypeErrorMessage),
        arguments(integerInputDefinition, null, null, integerTypeErrorMessage),
        arguments(integerInputDefinition, "", null, integerTypeErrorMessage),

        // INTEGER with min and/or max set (will cast correctly but may fail validation if out of
        // range)
        arguments(integerInputDefinitionRange, 123, 123, null),
        arguments(integerInputDefinitionRange, 200, 200, null),
        arguments(integerInputDefinitionRange, 0, 0, null),
        arguments(integerInputDefinitionRange, -5, -5, typeErrorMessageRange),
        arguments(integerInputDefinitionRange, 250, 250, typeErrorMessageRange),
        arguments(integerInputDefinitionOnlyMax, -5, -5, null),
        arguments(integerInputDefinitionOnlyMax, 250, 250, typeErrorMessageMax),
        arguments(integerInputDefinitionOnlyMin, 250, 250, null),
        arguments(integerInputDefinitionOnlyMin, -5, -5, typeErrorMessageMin),

        // FLOAT
        arguments(floatInputDefinition, 2.3, 2.3, null),
        arguments(floatInputDefinition, "2.3", 2.3, null),
        arguments(floatInputDefinition, "-2.3", -2.3, null),
        arguments(floatInputDefinition, "  2.3  ", 2.3, null),
        arguments(floatInputDefinition, "I am a string", null, floatTypeErrorMessage),
        arguments(floatInputDefinition, true, null, floatTypeErrorMessage),
        arguments(floatInputDefinition, null, null, floatTypeErrorMessage),
        arguments(floatInputDefinition, "", null, floatTypeErrorMessage),

        // FLOAT with min and/or max set (will cast correctly but may fail validation if out of
        // range)
        arguments(floatInputDefinitionRange, 2, 2.0, null),
        arguments(floatInputDefinitionRange, 3, 3.0, null),
        arguments(floatInputDefinitionRange, 0, 0.0, null),
        arguments(floatInputDefinitionRange, -5.4, -5.4, floatTypeErrorMessageRange),
        arguments(floatInputDefinitionRange, 5.4, 5.4, floatTypeErrorMessageRange),
        arguments(floatInputDefinitionOnlyMin, 5.4, 5.4, null),
        arguments(floatInputDefinitionOnlyMin, -5.2, -5.2, floatTypeErrorMessageMin),
        arguments(floatInputDefinitionOnlyMax, -5.2, -5.2, null),
        arguments(floatInputDefinitionOnlyMax, 5.2, 5.2, floatTypeErrorMessageMax),

        // STRING
        arguments(stringInputDefinition, "IAmAString", "IAmAString", null),
        arguments(stringInputDefinition, "I_Am-A=String.", "I_Am-A=String.", null),
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
        arguments(
            fileVcfInputDefinition,
            "path/to/!invalid/file.vcf.gz",
            "path/to/!invalid/file.vcf.gz",
            filePatternErrorMessage), // cast is successful but validation fails
        arguments(
            fileVcfInputDefinition,
            "gs://bucket_name/path/to/file.vcf.gz",
            "gs://bucket_name/path/to/file.vcf.gz",
            null),
        arguments(
            fileVcfInputDefinition,
            "path\\to\\valid\\windows\\file.vcf.gz",
            "path\\to\\valid\\windows\\file.vcf.gz",
            null),

        // MANIFEST
        arguments(manifestTsvInputDefinition, "path/to/file.tsv", "path/to/file.tsv", null),
        arguments(manifestTsvInputDefinition, "   path/to/file.tsv   ", "path/to/file.tsv", null),
        arguments(
            manifestTsvInputDefinition,
            "path/to/file.csv",
            "path/to/file.csv",
            tsvTypeErrorMessage), // cast is successful but validation fails
        arguments(manifestTsvInputDefinition, 3, null, tsvTypeErrorMessage),
        arguments(manifestTsvInputDefinition, null, null, tsvTypeErrorMessage),
        arguments(manifestTsvInputDefinition, "", null, tsvTypeErrorMessage),
        arguments(
            manifestTsvInputDefinition,
            "path/to/!invalid/file.tsv",
            "path/to/!invalid/file.tsv",
            filePatternErrorMessage), // cast is successful but validation fails
        arguments(
            manifestTsvInputDefinition,
            "gs://bucket_name/path/to/file.tsv",
            "gs://bucket_name/path/to/file.tsv",
            null),
        arguments(
            manifestTsvInputDefinition,
            "path\\to\\valid\\windows\\file.tsv",
            "path\\to\\valid\\windows\\file.tsv",
            null),

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
        arguments(
            stringArrayInputDefinition,
            Arrays.asList("validString1", "validString2", "validString3", "invalid string!"),
            Arrays.asList(
                "validString1",
                "validString2",
                "validString3",
                "invalid string!"), // cast is successful but validation fails
            stringArrayPatternErrorMessage),

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
            "%s must be an array of paths to files ending in .bed and containing only alphanumeric characters or the following symbols: -_.=:\\/"
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
            vcfArrayTypeErrorMessage),
        arguments(
            fileArrayVcfInputDefinition,
            List.of("valid/path/to/file.vcf.gz", "invalid!/path/to/file2.vcf.gz"),
            List.of(
                "valid/path/to/file.vcf.gz",
                "invalid!/path/to/file2.vcf.gz"), // cast is successful but validation fails
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

  @Test
  void isFileLike() {
    // manifests and files are files
    assertTrue(MANIFEST.isFileLike());
    assertTrue(FILE.isFileLike());

    // other types are not files
    assertFalse(STRING.isFileLike());
    assertFalse(INTEGER.isFileLike());
    assertFalse(STRING_ARRAY.isFileLike());
    assertFalse(
        FILE_ARRAY.isFileLike()); // for isFile purposes we don't yet support FILE_ARRAY processing
    assertFalse(FLOAT.isFileLike());
    assertFalse(BOOLEAN.isFileLike());
  }
}
