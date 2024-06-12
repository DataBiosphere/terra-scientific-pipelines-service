package bio.terra.pipelines.db.entities;

import bio.terra.common.exception.InternalServerErrorException;
import com.fasterxml.jackson.core.JsonProcessingException;
import com.fasterxml.jackson.core.type.TypeReference;
import com.fasterxml.jackson.databind.ObjectMapper;
import jakarta.persistence.AttributeConverter;
import java.util.Map;

// inspired by https://www.baeldung.com/jpa-persisting-enums-in-jpa
// @Converter(autoApply = true)
public class PipelineRunOutputConverter implements AttributeConverter<Map<String, String>, String> {

  ObjectMapper objectMapper = new ObjectMapper();

  @Override
  public String convertToDatabaseColumn(Map<String, String> pipelineOutputMap) {
    try {
      return objectMapper.writeValueAsString(pipelineOutputMap);
    } catch (JsonProcessingException e) {
      throw new InternalServerErrorException("Error converting pipeline outputs to string", e);
    }
  }

  @Override
  public Map<String, String> convertToEntityAttribute(String pipelineOutputString) {
    try {
      return objectMapper.readValue(pipelineOutputString, new TypeReference<>() {});
    } catch (JsonProcessingException e) {
      throw new InternalServerErrorException("Error extracting pipeline outputs to map", e);
    }
  }
}
