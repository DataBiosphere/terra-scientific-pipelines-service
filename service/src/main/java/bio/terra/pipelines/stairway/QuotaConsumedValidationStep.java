package bio.terra.pipelines.stairway;

import bio.terra.common.exception.InternalServerErrorException;
import bio.terra.pipelines.common.utils.FlightUtils;
import bio.terra.pipelines.common.utils.PipelinesEnum;
import bio.terra.pipelines.dependencies.rawls.RawlsService;
import bio.terra.pipelines.dependencies.rawls.RawlsServiceApiException;
import bio.terra.pipelines.dependencies.sam.SamService;
import bio.terra.pipelines.dependencies.stairway.JobMapKeys;
import bio.terra.pipelines.service.QuotasService;
import bio.terra.pipelines.stairway.imputation.ImputationJobMapKeys;
import bio.terra.rawls.model.Entity;
import bio.terra.stairway.*;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * This step calls Rawls to fetch outputs from a data table row for a given quota consumed job. It
 * specifically fetches the quota consumed value from the data table row using the quota_consumed
 * key.  If successful, it stores the quota consumed value in the working map.
 *
 * <p>This step expects quota consumed to be provided in the working map
 */
public class QuotaConsumedValidationStep implements Step {

    private final RawlsService rawlsService;
    private final SamService samService;
    private final QuotasService quotasService;
    private final Logger logger = LoggerFactory.getLogger(FetchQuotaConsumedFromDataTableStep.class);

    public QuotaConsumedValidationStep(RawlsService rawlsService, SamService samService, QuotasService quotasService) {
        this.rawlsService = rawlsService;
        this.samService = samService;
        this.quotasService = quotasService;
    }

    @Override
    @SuppressWarnings(
            "java:S2259") // suppress warning for possible NPE when calling pipelineName.getValue(),
    //  since we do validate that pipelineName is not null in `validateRequiredEntries`
    public StepResult doStep(FlightContext flightContext) {
        String jobId = flightContext.getFlightId();

        // validate and extract parameters from input map
        var inputParameters = flightContext.getInputParameters();
        FlightUtils.validateRequiredEntries(
                inputParameters,
                JobMapKeys.PIPELINE_NAME,
                JobMapKeys.USER_ID);

        PipelinesEnum pipelineName = inputParameters.get(JobMapKeys.PIPELINE_NAME, PipelinesEnum.class);
        String userId = inputParameters.get(JobMapKeys.USER_ID, String.class);

        // validate and extract parameters from working map
        FlightMap workingMap = flightContext.getWorkingMap();
        FlightUtils.validateRequiredEntries(
                workingMap,
                ImputationJobMapKeys.QUOTA_CONSUMED);

        int quotaConsumed = workingMap.get(ImputationJobMapKeys.QUOTA_CONSUMED, Integer.class);

        Entity entity = new Entity();

        // extract quota_consumed from entity
        try {
            quotaConsumed = (int) entity.getAttributes().get("quota_consumed");
            if (quotaConsumed <= 0) {
                return new StepResult(
                        StepStatus.STEP_RESULT_FAILURE_FATAL,
                        new InternalServerErrorException("Quota consumed is unexpectedly not greater than 0"));
            }
        } catch (NullPointerException e) {
            return new StepResult(
                    StepStatus.STEP_RESULT_FAILURE_FATAL,
                    new InternalServerErrorException("Quota consumed is unexpectedly null"));
        }

        // store the quota consumed value in the working map to use in a subsequent step
        workingMap.put(ImputationJobMapKeys.QUOTA_CONSUMED, quotaConsumed);

        return StepResult.getStepResultSuccess();
    }

    @Override
    public StepResult undoStep(FlightContext flightContext) {
        // nothing to undo
        return StepResult.getStepResultSuccess();
    }
}
