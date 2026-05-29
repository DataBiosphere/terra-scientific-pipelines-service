# Wdl Based Pipelines Overview

## Developer Information
### [Flight](https://github.com/DataBiosphere/terra-scientific-pipelines-service/tree/main/service/src/main/java/bio/terra/pipelines/stairway/flights/wdlbasedpipelinerun)
All wdl based pipelines use the same flight. This flight expects a variety of inputs so that it can be configured properly. Some of the more important ones.

* Pipeline Inputs
    * From pipeline_input_definitions table
* Pipeline outputs
    * From pipeline_output_definitions table
* User inputs
    * From pipeline_inputs table
* Quota consumed wdl
    * Generated using the [ToolConfigService](https://github.com/DataBiosphere/terra-scientific-pipelines-service/blob/main/service/src/main/java/bio/terra/pipelines/service/ToolConfigService.java) while creating the [flight](https://github.com/DataBiosphere/terra-scientific-pipelines-service/blob/main/service/src/main/java/bio/terra/pipelines/service/PipelineRunsService.java) (in the startPipelineRun function)
* QC wdl
    * Same as Quota consumed wdl
* Pipeline Wdl
    * Same as quota consumed wdl
* Pipeline version
    * From the prepare payload (stored in the pipeline_runs table)

This flight will first write all the inputs the wdl(s) need to the data tables in the runner workspace. It will then launch a quota consumed wdl submission to make sure the user has enough quota to run their pipeline.  Then it launches a QC wdl submission to make sure the users inputs pass all of our QC checks.  It will then launch the pipeline wdl.  Once that completes it will write the outputs to the pipeline_outputs table and mark the pipeline run as a success.  Any differences between what wdls to run and what inputs/outputs to expect are all driven by the pipeline version supplied to the flight.
