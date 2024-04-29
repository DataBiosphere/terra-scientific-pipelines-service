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

Validation is performed in PipelinesService.createJob(), which calls validateInputs(), which calls both 
validateRequiredInputs() and validateInputTypes(). These methods return a list of error messages, if any, which are 
accumulated into a ValidationException that is thrown inside validateInputs(). Because ValidationException ultimately 
extends ErrorReportException, which is handled by the global exception handler, these error messages are returned 
directly to the caller in the response. 

Inside validateInputTypes(), we call the cast() method of the appropriate input type on the input value. The cast method 
can throw either a ValidationException or an IllegalArgumentException. ValidationExceptions contain custom 
error messages that we return directly to the caller. IllegalArgumentExceptions, which are thrown by
the objectMapper, are caught and translated into a generic "wrong type" error message. 
