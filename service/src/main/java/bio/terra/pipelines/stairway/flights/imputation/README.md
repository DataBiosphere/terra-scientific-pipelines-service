# Flight class versioning guidelines

In Stairway, a Flight is a sequence of Steps that are executed in order to perform a specific task. Over time, as 
requirements change or new features are added, it may be necessary to update or modify existing Flight classes. 
To manage these changes effectively, it is important to implement a versioning strategy for Flights.

Breaking changes to a Flight can cause the entire application to fail upon release, if any flights are in progress at the time.

Breaking changes to Flights include but may not be limited to:
- new required input parameters
- new Steps
- removing Steps
- changing the order of Steps
A good rule of thumb is that if you are changing the Flight class itself, it is likely a breaking change. 
Changes to Steps are unlikely to cause the application to fail but could in theory cause a running flight to fail.

To avoid this, we recommend the following versioning strategy:
When making breaking changes to a Flight, create a new version of the Flight class with a new subpackage name with a date format following `vYYYYMMDD`. 
Reference the new Flight in the code that will use it going forward (currently, this is [PipelineRunsService](service/src/main/java/bio/terra/pipelines/service/PipelineRunsService.java)). 
Do not delete the previous Flight version.

Create a new test directory for this version and make updates to the Flight's test as needed.

The old Flight and its test should remain in the codebase to support any in-progress executions, but they can be 
removed in a future PR/release once the old Flight is no longer in use.
