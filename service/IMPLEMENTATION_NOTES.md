# Implementation Notes

## Input validation

We validate user-provided inputs to the pipeline (supplied via the "pipelineInputs" parameter in the createJob payload).

We support required and optional inputs. 
Required inputs:
- must not be null
- must not be empty, if they are strings or arrays
- must be the correct type
Optional inputs:
- will have a default value, if not provided by the user
- must not be null, if provided
- must be the correct type, if provided
Extra inputs that are not defined in the pipeline's input schema are considered an error. 
All user-provided file inputs must be either local file paths or GS URIs (beginning with "gs://").
Leading and trailing whitespace is trimmed from string inputs (including VCF inputs).

Validation is called by ApiPipelinesController.createJob() and performed in PipelinesService.validateInputs().
validateInputs() calls both validateRequiredInputs() and validateInputTypes(). 
These methods return a list of error messages, if any, that are accumulated into a ValidationException that is 
thrown inside validateInputs(). Because ValidationException ultimately extends ErrorReportException, 
which is handled by the global exception handler, these error messages are returned directly to the caller 
in the API response.

## PipelineRun Status updates, Metrics, and Stairway hooks
We use our pipelineRuns table as the source of truth for our pipeline runs, and we pull in the pipelineRun 
status from Stairway when the pipelineRun completes. In the case of a successful pipelineRun, we update its 
status in the final step of the stairway flight. In the case of a failed pipelineRun, we update its status 
using a Stairway hook ([StairwaySetPipelineRunStatusHook](src/main/java/bio/terra/pipelines/common/utils/StairwaySetPipelineRunStatusHook.java)).

We use a Stairway hook so that the pipelineRun will be marked as failed regardless of when in the flight the 
failure occurs, and whether it's a roll-back-able error or a dismal failure.

Similarly, for metrics reporting, we use a Stairway hook ([StairwayFailedMetricsCounterHook](src/main/java/bio/terra/pipelines/common/utils/StairwayFailedMetricsCounterHook.java)) 
to increment the failed metrics counter when a pipelineRun fails.

Because Stairway hooks are applied at the Stairway instance level and not per-flight, we conditionally run the
logic in each of these hooks only if the corresponding flight map key is present and set to true in the flight map.
