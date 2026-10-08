package bio.terra.pipelines.db.entities.converters;

import bio.terra.pipelines.common.utils.QuotaAllocationStatusEnum;
import jakarta.persistence.AttributeConverter;
import jakarta.persistence.Converter;

// inspired by https://www.baeldung.com/jpa-persisting-enums-in-jpa
@Converter(autoApply = true)
public class QuotaAllocationStatusConverter
    implements AttributeConverter<QuotaAllocationStatusEnum, String> {
  @Override
  public String convertToDatabaseColumn(QuotaAllocationStatusEnum quotaAllocationStatusEnum) {
    return quotaAllocationStatusEnum.toString();
  }

  @Override
  public QuotaAllocationStatusEnum convertToEntityAttribute(String quotaAllocationStatusString) {
    return QuotaAllocationStatusEnum.valueOf(quotaAllocationStatusString);
  }
}
