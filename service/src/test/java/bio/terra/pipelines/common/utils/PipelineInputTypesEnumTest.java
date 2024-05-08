package bio.terra.pipelines.common.utils;

import static bio.terra.pipelines.common.utils.PipelineInputTypesEnum.INTEGER;
import static bio.terra.pipelines.common.utils.PipelineInputTypesEnum.STRING;
import static bio.terra.pipelines.common.utils.PipelineInputTypesEnum.STRING_ARRAY;
import static bio.terra.pipelines.common.utils.PipelineInputTypesEnum.VCF;
import static bio.terra.pipelines.common.utils.PipelineInputTypesEnum.VCF_ARRAY;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.params.provider.Arguments.arguments;

import bio.terra.common.exception.ValidationException;
import bio.terra.pipelines.testutils.BaseTest;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Stream;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.Arguments;
import org.junit.jupiter.params.provider.MethodSource;

class PipelineInputTypesEnumTest extends BaseTest {

  @Test
  void testPipelineInputTypesEnum() {
    STRING.cast("fieldName", "value");
    INTEGER.cast("fieldName", 1);
    VCF.cast("fieldName", "value.vcf.gz");
    STRING_ARRAY.cast("fieldName", List.of("value1", "value2"));
    VCF_ARRAY.cast("fieldName", List.of("value1.vcf.gz", "value2.vcf.gz"));
  }

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
    String emptyArrayErrorMessage = "input_name must not be an empty list";

    return Stream.of(
        // arguments: type enum, input value to cast, expected cast value if successful, error
        // message if fails

        // INTEGER
        arguments(INTEGER, 123, 123, null),
        arguments(INTEGER, "123", 123, null),
        arguments(INTEGER, "I am a string", null, integerTypeErrorMessage),
        arguments(INTEGER, 2.3, null, integerTypeErrorMessage),
        arguments(INTEGER, "2.3", null, integerTypeErrorMessage),
        arguments(INTEGER, null, null, notNullErrorMessage),
        arguments(INTEGER, "", null, notEmptyErrorMessage),

        // STRING
        arguments(STRING, "I am a string", "I am a string", null),
        arguments(STRING, "\"I am a string\"", "\"I am a string\"", null),
        arguments(STRING, "$tr1nG.w1th.0th3r.ch@r@ct3r$", "$tr1nG.w1th.0th3r.ch@r@ct3r$", null),
        arguments(STRING, "    I am a string    ", "I am a string", null),
        arguments(
            STRING, List.of("this", "is", "not", "a", "string"), null, stringTypeErrorMessage),
        arguments(STRING, 123, null, stringTypeErrorMessage),
        arguments(STRING, null, null, notNullErrorMessage),
        arguments(STRING, "", null, notEmptyErrorMessage),

        // VCF
        arguments(VCF, "path/to/file.vcf.gz", "path/to/file.vcf.gz", null),
        arguments(VCF, "   path/to/file.vcf.gz   ", "path/to/file.vcf.gz", null),
        arguments(VCF, "path/to/file.vcf", null, vcfTypeErrorMessage),
        arguments(VCF, 3, null, stringTypeErrorMessage),
        arguments(VCF, null, null, notNullErrorMessage),
        arguments(VCF, "", null, notEmptyErrorMessage),

        // STRING_ARRAY
        arguments(
            STRING_ARRAY,
            List.of("this", "is", "a", "list", "of", "strings"),
            List.of("this", "is", "a", "list", "of", "strings"),
            null),
        arguments(
            STRING_ARRAY,
            List.of("this ", " is", " a ", "  list", "of  ", "  strings  "),
            List.of("this", "is", "a", "list", "of", "strings"),
            null),
        arguments(STRING_ARRAY, "I am not an array", null, stringArrayTypeErrorMessage),
        arguments(STRING_ARRAY, Arrays.asList("string", 2, 3), null, stringArrayTypeErrorMessage),
        arguments(STRING_ARRAY, Arrays.asList(1, 2, 3), null, stringArrayTypeErrorMessage),
        arguments(STRING_ARRAY, null, null, notNullErrorMessage),
        arguments(STRING_ARRAY, List.of(), null, stringArrayTypeErrorMessage),
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
        arguments(VCF_ARRAY, List.of("path/to/file.vcf.gz"), List.of("path/to/file.vcf.gz"), null),
        arguments(
            VCF_ARRAY, List.of("  path/to/file.vcf.gz  "), List.of("path/to/file.vcf.gz"), null),
        arguments(VCF_ARRAY, "this/is/not/an/array.vcf.gz", null, vcfArrayTypeErrorMessage),
        arguments(
            VCF_ARRAY,
            Arrays.asList("path/to/file.vcf.gz", "not a path"),
            null,
            vcfArrayTypeErrorMessage),
        arguments(
            VCF_ARRAY, Arrays.asList("path/to/file.vcf.gz", 2.5), null, vcfArrayTypeErrorMessage),
        arguments(VCF_ARRAY, null, null, notNullErrorMessage),
        arguments(VCF_ARRAY, List.of(), null, vcfArrayTypeErrorMessage),
        arguments(
            VCF_ARRAY,
            Arrays.asList(null, "list/with/null.vcf.gz"),
            null,
            vcfArrayTypeErrorMessage));
  }

  @ParameterizedTest
  @MethodSource("castValidations")
  void castInputValues(
      PipelineInputTypesEnum inputType,
      Object inputValue,
      Object expectedCastValue,
      String expectedErrorMessage) {
    if (!(expectedCastValue == null)) {
      assertEquals(expectedCastValue, inputType.cast("fieldName", inputValue));
    } else {
      ValidationException exception =
          assertThrows(ValidationException.class, () -> inputType.cast("input_name", inputValue));
      assertEquals(expectedErrorMessage, exception.getMessage());
    }
  }
}
