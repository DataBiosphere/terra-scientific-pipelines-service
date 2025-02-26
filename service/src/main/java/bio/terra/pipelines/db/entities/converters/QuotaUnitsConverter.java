package bio.terra.pipelines.db.entities.converters;

import bio.terra.pipelines.common.utils.QuotaUnitsEnum;
import jakarta.persistence.AttributeConverter;
import jakarta.persistence.Converter;

// inspired by https://www.baeldung.com/jpa-persisting-enums-in-jpa
@Converter(autoApply = true)
public class QuotaUnitsConverter implements AttributeConverter<QuotaUnitsEnum, String> {
  @Override
  public String convertToDatabaseColumn(QuotaUnitsEnum quotaUnitsEnum) {
    return quotaUnitsEnum.getValue();
  }

  @Override
  public QuotaUnitsEnum convertToEntityAttribute(String quotaUnitsString) {
    return QuotaUnitsEnum.valueOf(quotaUnitsString.toUpperCase());
  }
}
