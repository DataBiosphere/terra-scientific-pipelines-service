package bio.terra.pipelines.common.utils;

import static bio.terra.pipelines.common.utils.PipelineInputTypesEnum.INTEGER;
import static bio.terra.pipelines.common.utils.PipelineInputTypesEnum.STRING;
import static bio.terra.pipelines.common.utils.PipelineInputTypesEnum.STRING_ARRAY;
import static bio.terra.pipelines.common.utils.PipelineInputTypesEnum.VCF;
import static bio.terra.pipelines.common.utils.PipelineInputTypesEnum.VCF_ARRAY;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertTrue;
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

class PipelineInputTypesEnumTest extends BaseTest {

  private static Stream<Arguments> castValidations() {
    // error messages
    String stringTypeErrorMessage = "input_name must be a string";
    String integerTypeErrorMessage = "input_name must be an integer";
    String vcfTypeErrorMessage = "input_name must be a path to a VCF file ending in .vcf.gz";
    String stringArrayTypeErrorMessage = "input_name must be an array of strings";
    String vcfArrayTypeErrorMessage =
        "input_name must be an array of paths to VCF files ending in .vcf.gz";
    String notNullErrorMessage = "input_name must not be null";
    String notEmptyErrorMessage = "input_name must not be empty";
    String notEmptyListErrorMessage = "input_name must not be an empty list";

    return Stream.of(
        // arguments: type enum, input value to cast, expected cast value if successful, error
        // message if fails

        // INTEGER
        arguments(INTEGER, 123, 123, null),
        arguments(INTEGER, "123", 123, null),
        arguments(INTEGER, "I am a string", null, integerTypeErrorMessage),
        arguments(INTEGER, 2.3, null, integerTypeErrorMessage),
        arguments(INTEGER, "2.3", null, integerTypeErrorMessage),
        arguments(INTEGER, null, null, integerTypeErrorMessage),
        arguments(INTEGER, "", null, integerTypeErrorMessage),

        // STRING
        arguments(STRING, "I am a string", "I am a string", null),
        arguments(STRING, "\"I am a string\"", "\"I am a string\"", null),
        arguments(STRING, "$tr1nG.w1th.0th3r.ch@r@ct3r$", "$tr1nG.w1th.0th3r.ch@r@ct3r$", null),
        arguments(
            STRING,
            "    I am a string with extra whitespace    ",
            "I am a string with extra whitespace",
            null),
        arguments(
            STRING, List.of("this", "is", "not", "a", "string"), null, stringTypeErrorMessage),
        arguments(STRING, 123, null, stringTypeErrorMessage),
        arguments(STRING, null, null, stringTypeErrorMessage),
        arguments(STRING, "", null, stringTypeErrorMessage),

        // VCF
        arguments(VCF, "path/to/file.vcf.gz", "path/to/file.vcf.gz", null),
        arguments(VCF, "   path/to/file.vcf.gz   ", "path/to/file.vcf.gz", null),
        arguments(
            VCF,
            "path/to/file.vcf",
            "path/to/file.vcf",
            vcfTypeErrorMessage), // cast is successful but validation fails
        arguments(VCF, 3, null, vcfTypeErrorMessage),
        arguments(VCF, null, null, vcfTypeErrorMessage),
        arguments(VCF, "", null, vcfTypeErrorMessage),

        // STRING_ARRAY
        arguments(
            STRING_ARRAY,
            List.of("this", "is", "a", "list", "of", "strings"),
            List.of("this", "is", "a", "list", "of", "strings"),
            null),
        arguments(STRING_ARRAY, List.of("singleton list"), List.of("singleton list"), null),
        arguments(
            STRING_ARRAY,
            List.of(
                "this ", " is", " a ", "  list", "of  ", "  strings  ", "with extra whitespace"),
            List.of("this", "is", "a", "list", "of", "strings", "with extra whitespace"),
            null),
        arguments(
            STRING_ARRAY,
            "[\"this\", \"is\", \"a\", \"stringy\", \"list\",  \"of\", \"strings\"]",
            List.of("this", "is", "a", "stringy", "list", "of", "strings"),
            null),
        arguments(
            STRING_ARRAY,
            "[\" stringy \", \" and\",  \"whitespace \"]",
            List.of("stringy", "and", "whitespace"),
            null),
        arguments(STRING_ARRAY, "I am not an array", null, stringArrayTypeErrorMessage),
        arguments(STRING_ARRAY, Arrays.asList("string", 2, 3), null, stringArrayTypeErrorMessage),
        arguments(STRING_ARRAY, Arrays.asList(1, 2, 3), null, stringArrayTypeErrorMessage),
        arguments(STRING_ARRAY, null, null, stringArrayTypeErrorMessage),
        arguments(
            STRING_ARRAY,
            Collections.emptyList(),
            Collections.emptyList(), // cast is successful but validation fails
            stringArrayTypeErrorMessage),
        arguments(
            STRING_ARRAY,
            Arrays.asList("", "array with empty string"),
            null,
            stringArrayTypeErrorMessage),
        arguments(
            STRING_ARRAY,
            Arrays.asList(null, "array with null"),
            null,
            stringArrayTypeErrorMessage),

        // VCF_ARRAY
        arguments(
            VCF_ARRAY,
            List.of("path/to/file.vcf.gz", "path/to/file2.vcf.gz"),
            List.of("path/to/file.vcf.gz", "path/to/file2.vcf.gz"),
            null),
        arguments(
            VCF_ARRAY,
            List.of("  path/to/file/with/extra/whitespace.vcf.gz  "),
            List.of("path/to/file/with/extra/whitespace.vcf.gz"),
            null),
        arguments(
            VCF_ARRAY,
            "[\"path/to/file1.vcf.gz\", \"path/to/file2.vcf.gz\"]",
            List.of("path/to/file1.vcf.gz", "path/to/file2.vcf.gz"),
            null),
        arguments(
            VCF_ARRAY,
            "[\" path/to/file1.vcf.gz \", \"path/to/file2.vcf.gz \"]",
            List.of("path/to/file1.vcf.gz", "path/to/file2.vcf.gz"),
            null),
        arguments(VCF_ARRAY, "this/is/not/an/array.vcf.gz", null, vcfArrayTypeErrorMessage),
        arguments(
            VCF_ARRAY,
            Arrays.asList("path/to/file.vcf.gz", "just/a/string"),
            Arrays.asList(
                "path/to/file.vcf.gz", "just/a/string"), // cast is successful but validation fails
            vcfArrayTypeErrorMessage),
        arguments(
            VCF_ARRAY, Arrays.asList("path/to/file.vcf.gz", 2.5), null, vcfArrayTypeErrorMessage),
        arguments(VCF_ARRAY, null, null, vcfArrayTypeErrorMessage),
        arguments(
            VCF_ARRAY,
            Collections.emptyList(),
            Collections.emptyList(),
            vcfArrayTypeErrorMessage), // cast is successful but validation fails
        arguments(
            VCF_ARRAY,
            Arrays.asList(null, "list/with/null.vcf.gz"),
            null,
            vcfArrayTypeErrorMessage));
  }

  @ParameterizedTest
  @MethodSource("castValidations")
  <T> void castInputValues(
      PipelineInputTypesEnum inputType, Object inputValue, T expectedCastValue) {
    if (expectedCastValue != null) {
      if (expectedCastValue instanceof List expectedListCastValue) {
        List<String> listCastValue =
            inputType.cast("fieldName", inputValue, new TypeReference<>() {});
        assertTrue(expectedListCastValue.containsAll(listCastValue));
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
      PipelineInputTypesEnum inputType,
      Object inputValue,
      T expectedCastValue,
      String expectedErrorMessage) {
    if (expectedErrorMessage == null) {
      // should validate
      assertTrue(inputType.validate("input_name", inputValue).isEmpty());
    } else {
      // should not validate
      assertEquals(expectedErrorMessage, inputType.validate("input_name", inputValue).get());
    }
  }
}
