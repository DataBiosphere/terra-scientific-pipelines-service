package bio.terra.pipelines.common.utils;

import static bio.terra.pipelines.common.utils.PipelineVariableTypesEnum.FILE;
import static bio.terra.pipelines.common.utils.PipelineVariableTypesEnum.FILE_ARRAY;
import static bio.terra.pipelines.common.utils.PipelineVariableTypesEnum.INTEGER;
import static bio.terra.pipelines.common.utils.PipelineVariableTypesEnum.STRING;
import static bio.terra.pipelines.common.utils.PipelineVariableTypesEnum.STRING_ARRAY;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.params.provider.Arguments.arguments;

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
    // error messages
    String stringTypeErrorMessage = "input_name must be a string";
    String integerTypeErrorMessage = "input_name must be an integer";
    String vcfTypeErrorMessage = "input_name must be a path to a file ending in .vcf.gz";
    String stringArrayTypeErrorMessage = "input_name must be an array of strings";
    String vcfArrayTypeErrorMessage =
        "input_name must be an array of paths to files ending in .vcf.gz";
    String notNullOrEmptyErrorMessage = "input_name must not be null or empty";

    return Stream.of(
        // arguments: type enum, file suffix, input value to cast, expected cast value if
        // successful, error
        // message if fails

        // INTEGER
        arguments(INTEGER, null, 123, 123, null),
        arguments(INTEGER, null, "123", 123, null),
        arguments(INTEGER, null, "I am a string", null, integerTypeErrorMessage),
        arguments(INTEGER, null, 2.3, null, integerTypeErrorMessage),
        arguments(INTEGER, null, "2.3", null, integerTypeErrorMessage),
        arguments(INTEGER, null, null, null, integerTypeErrorMessage),
        arguments(INTEGER, null, "", null, integerTypeErrorMessage),

        // STRING
        arguments(STRING, null, "I am a string", "I am a string", null),
        arguments(STRING, null, "\"I am a string\"", "\"I am a string\"", null),
        arguments(
            STRING, null, "$tr1nG.w1th.0th3r.ch@r@ct3r$", "$tr1nG.w1th.0th3r.ch@r@ct3r$", null),
        arguments(
            STRING,
            null,
            "    I am a string with extra whitespace    ",
            "I am a string with extra whitespace",
            null),
        arguments(
            STRING,
            null,
            List.of("this", "is", "not", "a", "string"),
            null,
            stringTypeErrorMessage),
        arguments(STRING, null, 123, null, stringTypeErrorMessage),
        arguments(STRING, null, null, null, stringTypeErrorMessage),
        arguments(STRING, null, "", null, stringTypeErrorMessage),

        // FILE
        arguments(FILE, ".vcf.gz", "path/to/file.vcf.gz", "path/to/file.vcf.gz", null),
        arguments(FILE, ".vcf.gz", "   path/to/file.vcf.gz   ", "path/to/file.vcf.gz", null),
        arguments(
            FILE,
            ".vcf.gz",
            "path/to/file.vcf",
            "path/to/file.vcf",
            vcfTypeErrorMessage), // cast is successful but validation fails
        arguments(FILE, ".vcf.gz", 3, null, vcfTypeErrorMessage),
        arguments(FILE, ".vcf.gz", null, null, vcfTypeErrorMessage),
        arguments(FILE, ".vcf.gz", "", null, vcfTypeErrorMessage),
        arguments(
            FILE,
            ".bed",
            "path/to/not/a/bed/file",
            "path/to/not/a/bed/file",
            "input_name must be a path to a file ending in .bed"), // cast is successful but
        // validation fails

        // STRING_ARRAY
        arguments(
            STRING_ARRAY,
            null,
            List.of("this", "is", "a", "list", "of", "strings"),
            List.of("this", "is", "a", "list", "of", "strings"),
            null),
        arguments(STRING_ARRAY, null, List.of("singleton list"), List.of("singleton list"), null),
        arguments(
            STRING_ARRAY,
            null,
            List.of(
                "this ", " is", " a ", "  list", "of  ", "  strings  ", "with extra whitespace"),
            List.of("this", "is", "a", "list", "of", "strings", "with extra whitespace"),
            null),
        arguments(
            STRING_ARRAY,
            null,
            "[\"this\", \"is\", \"a\", \"stringy\", \"list\",  \"of\", \"strings\"]",
            List.of("this", "is", "a", "stringy", "list", "of", "strings"),
            null),
        arguments(
            STRING_ARRAY,
            null,
            "[\" stringy \", \" and\",  \"whitespace \"]",
            List.of("stringy", "and", "whitespace"),
            null),
        arguments(STRING_ARRAY, null, "I am not an array", null, stringArrayTypeErrorMessage),
        arguments(
            STRING_ARRAY, null, Arrays.asList("string", 2, 3), null, stringArrayTypeErrorMessage),
        arguments(STRING_ARRAY, null, Arrays.asList(1, 2, 3), null, stringArrayTypeErrorMessage),
        arguments(STRING_ARRAY, null, null, null, notNullOrEmptyErrorMessage),
        arguments(
            STRING_ARRAY,
            null,
            Collections.emptyList(),
            Collections.emptyList(), // cast is successful but validation fails
            notNullOrEmptyErrorMessage),
        arguments(
            STRING_ARRAY,
            null,
            Arrays.asList("", "array with empty string"),
            null,
            stringArrayTypeErrorMessage),
        arguments(
            STRING_ARRAY,
            null,
            Arrays.asList(null, "array with null"),
            null,
            stringArrayTypeErrorMessage),

        // FILE_ARRAY
        arguments(
            FILE_ARRAY,
            ".vcf.gz",
            List.of("path/to/file.vcf.gz", "path/to/file2.vcf.gz"),
            List.of("path/to/file.vcf.gz", "path/to/file2.vcf.gz"),
            null),
        arguments(
            FILE_ARRAY,
            ".vcf.gz",
            List.of("  path/to/file/with/extra/whitespace.vcf.gz  "),
            List.of("path/to/file/with/extra/whitespace.vcf.gz"),
            null),
        arguments(
            FILE_ARRAY,
            ".vcf.gz",
            "[\"path/to/file1.vcf.gz\", \"path/to/file2.vcf.gz\"]",
            List.of("path/to/file1.vcf.gz", "path/to/file2.vcf.gz"),
            null),
        arguments(
            FILE_ARRAY,
            ".vcf.gz",
            "[\" path/to/file1.vcf.gz \", \"path/to/file2.vcf.gz \"]",
            List.of("path/to/file1.vcf.gz", "path/to/file2.vcf.gz"),
            null),
        arguments(
            FILE_ARRAY, ".vcf.gz", "this/is/not/an/array.vcf.gz", null, vcfArrayTypeErrorMessage),
        arguments(
            FILE_ARRAY,
            ".vcf.gz",
            Arrays.asList("path/to/file.vcf.gz", "just/a/string"),
            Arrays.asList(
                "path/to/file.vcf.gz", "just/a/string"), // cast is successful but validation fails
            vcfArrayTypeErrorMessage),
        arguments(
            FILE_ARRAY,
            ".vcf.gz",
            Arrays.asList("path/to/file.vcf.gz", 2.5),
            null,
            vcfArrayTypeErrorMessage),
        arguments(
            FILE_ARRAY,
            ".bed",
            Arrays.asList("file.vcf.gz", "file.bed"),
            Arrays.asList("file.vcf.gz", "file.bed"),
            "input_name must be an array of paths to files ending in .bed"), // cast is successful
        // but validation fails
        arguments(FILE_ARRAY, ".vcf.gz", null, null, notNullOrEmptyErrorMessage),
        arguments(
            FILE_ARRAY,
            ".vcf.gz",
            Collections.emptyList(),
            Collections.emptyList(),
            notNullOrEmptyErrorMessage), // cast is successful but validation fails
        arguments(
            FILE_ARRAY,
            ".vcf.gz",
            Arrays.asList(null, "list/with/null.vcf.gz"),
            null,
            vcfArrayTypeErrorMessage));
  }

  @ParameterizedTest
  @MethodSource("castValidations")
  <T> void castInputValues(
      PipelineVariableTypesEnum inputType,
      String fileSuffix,
      Object inputValue,
      T expectedCastValue) {
    if (expectedCastValue != null) {
      if (expectedCastValue instanceof List expectedListCastValue) {
        List<String> listCastValue =
            inputType.cast("fieldName", inputValue, new TypeReference<>() {});

        assertEquals(expectedCastValue, listCastValue);
        // Ensure that base class matches up
        if (!expectedListCastValue.isEmpty()) {
          assertEquals(expectedListCastValue.get(0).getClass(), listCastValue.get(0).getClass());
        }
      } else {
        assertEquals(
            expectedCastValue, inputType.cast("fieldName", inputValue, new TypeReference<>() {}));
        // Ensure that class matches up
        assertEquals(
            expectedCastValue.getClass(),
            inputType.cast("fieldName", inputValue, new TypeReference<>() {}).getClass());
      }
    } else {
      // cast should return null
      assertNull(inputType.cast("fieldName", inputValue, new TypeReference<>() {}));
    }
  }

  @ParameterizedTest
  @MethodSource("castValidations")
  <T> void validateInputValues(
      PipelineVariableTypesEnum inputType,
      String fileSuffix,
      Object inputValue,
      T expectedCastValue,
      String expectedErrorMessage) {
    if (expectedErrorMessage == null) {
      // should validate
      assertNull(inputType.validate("input_name", fileSuffix, inputValue));
    } else {
      // should not validate
      assertEquals(expectedErrorMessage, inputType.validate("input_name", fileSuffix, inputValue));
    }
  }
}
