package bio.terra.pipelines.db.entities.converters;

import bio.terra.pipelines.common.utils.QuotaAllocationSourceEnum;
import jakarta.persistence.AttributeConverter;
import jakarta.persistence.Converter;

// inspired by https://www.baeldung.com/jpa-persisting-enums-in-jpa
@Converter(autoApply = true)
public class QuotaAllocationSourceConverter
    implements AttributeConverter<QuotaAllocationSourceEnum, String> {
  @Override
  public String convertToDatabaseColumn(QuotaAllocationSourceEnum quotaAllocationSourceEnum) {
    return quotaAllocationSourceEnum.getValue();
  }

  @Override
  public QuotaAllocationSourceEnum convertToEntityAttribute(String quotaSourceString) {
    return QuotaAllocationSourceEnum.valueOf(quotaSourceString.toUpperCase());
  }
}
