package bio.terra.pipelines.common.utils;

import static org.junit.jupiter.api.Assertions.*;
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
    PipelineInputTypesEnum.STRING.cast("fieldName", "value");
    PipelineInputTypesEnum.INTEGER.cast("fieldName", 1);
    PipelineInputTypesEnum.VCF.cast("fieldName", "value.vcf.gz");
    PipelineInputTypesEnum.STRING_ARRAY.cast("fieldName", List.of("value1", "value2"));
    PipelineInputTypesEnum.VCF_ARRAY.cast("fieldName", List.of("value1.vcf.gz", "value2.vcf.gz"));
  }

  private static Stream<Arguments> castValidations() {
    return Stream.of(
        // arguments: type enum, input value to cast, expected cast value
        arguments(PipelineInputTypesEnum.INTEGER, 123, 123),
        arguments(PipelineInputTypesEnum.INTEGER, "123", 123),
        arguments(PipelineInputTypesEnum.STRING, "I am a string", "I am a string"),
        arguments(PipelineInputTypesEnum.STRING, "\"I am a string\"", "\"I am a string\""),
        arguments(PipelineInputTypesEnum.VCF, "path/to/file.vcf.gz", "path/to/file.vcf.gz"),
        arguments(
            PipelineInputTypesEnum.STRING_ARRAY,
            List.of("this", "is", "a", "list", "of", "strings"),
            List.of("this", "is", "a", "list", "of", "strings")),
        arguments(
            PipelineInputTypesEnum.VCF_ARRAY,
            List.of("path/to/file.vcf.gz"),
            List.of("path/to/file.vcf.gz")));
  }

  @ParameterizedTest
  @MethodSource("castValidations")
  void castInputValueOk(
      PipelineInputTypesEnum inputType, Object inputValue, Object expectedCastValue) {
    assertEquals(expectedCastValue, inputType.cast("fieldName", inputValue));
  }

  private static Stream<Arguments> castValidationFailures() {
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
        // arguments: type enum, input value to cast, expected error message
        // basic type checks that should fail validation (produce an error message)
        arguments(
            PipelineInputTypesEnum.STRING,
            List.of("this", "is", "not", "a", "string"),
            stringTypeErrorMessage),
        arguments(PipelineInputTypesEnum.STRING, 123, stringTypeErrorMessage),
        arguments(PipelineInputTypesEnum.INTEGER, "I am a string", integerTypeErrorMessage),
        arguments(PipelineInputTypesEnum.INTEGER, 2.3, integerTypeErrorMessage),
        arguments(PipelineInputTypesEnum.INTEGER, "2.3", integerTypeErrorMessage),
        arguments(PipelineInputTypesEnum.VCF, "path/to/file.vcf", vcfTypeErrorMessage),
        arguments(PipelineInputTypesEnum.VCF, 3, stringTypeErrorMessage),
        arguments(
            PipelineInputTypesEnum.STRING_ARRAY, "I am not an array", stringArrayTypeErrorMessage),
        arguments(
            PipelineInputTypesEnum.STRING_ARRAY,
            Arrays.asList(1, 2, 3),
            stringArrayTypeErrorMessage),
        arguments(
            PipelineInputTypesEnum.VCF_ARRAY,
            "this/is/not/an/array.vcf.gz",
            vcfArrayTypeErrorMessage),
        arguments(
            PipelineInputTypesEnum.VCF_ARRAY,
            Arrays.asList("path/to/file.vcf.gz", "not a path"),
            vcfArrayTypeErrorMessage),
        // null and empty checks
        arguments(PipelineInputTypesEnum.STRING, null, notNullErrorMessage),
        arguments(PipelineInputTypesEnum.STRING, "", notEmptyErrorMessage),
        arguments(PipelineInputTypesEnum.INTEGER, null, notNullErrorMessage),
        arguments(PipelineInputTypesEnum.INTEGER, "", notEmptyErrorMessage),
        arguments(PipelineInputTypesEnum.STRING_ARRAY, null, notNullErrorMessage),
        arguments(PipelineInputTypesEnum.STRING_ARRAY, List.of(), emptyArrayErrorMessage),
        arguments(
            PipelineInputTypesEnum.STRING_ARRAY,
            Arrays.asList("", "array with empty string"),
            stringArrayTypeErrorMessage),
        arguments(
            PipelineInputTypesEnum.STRING_ARRAY,
            Arrays.asList(null, "array with null"),
            stringArrayTypeErrorMessage),
        arguments(PipelineInputTypesEnum.VCF_ARRAY, null, notNullErrorMessage),
        arguments(PipelineInputTypesEnum.VCF_ARRAY, List.of(), emptyArrayErrorMessage),
        arguments(
            PipelineInputTypesEnum.VCF_ARRAY,
            Arrays.asList(null, "list/with/null.vcf.gz"),
            vcfArrayTypeErrorMessage));
  }

  @ParameterizedTest
  @MethodSource("castValidationFailures")
  void castInputValueFail(
      PipelineInputTypesEnum inputType, Object inputValue, String expectedErrorMessage) {
    ValidationException exception =
        assertThrows(ValidationException.class, () -> inputType.cast("input_name", inputValue));
    assertEquals(expectedErrorMessage, exception.getMessage());
  }
}
