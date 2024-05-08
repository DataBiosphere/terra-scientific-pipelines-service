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
Extra inputs that are not defined in the pipeline's input schema are ignored and logged at the WARN level. 
Leading and trailing whitespace is trimmed from string inputs (including VCF inputs).

Validation is called by ApiPipelinesController.createJob() and performed in PipelinesService.validateInputs().
validateInputs() calls both validateRequiredInputs() and validateInputTypes(). 
These methods return a list of error messages, if any, that are accumulated into a ValidationException that is 
thrown inside validateInputs(). Because ValidationException ultimately extends ErrorReportException, 
which is handled by the global exception handler, these error messages are returned directly to the caller 
in the API response.
