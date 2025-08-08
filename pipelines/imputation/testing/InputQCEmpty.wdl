version 1.0

workflow InputQC {
    String pipeline_version = "0.0.1"

    input {
        File multi_sample_vcf
        Array[String] contigs
    }

    call ReturnBoolAndString

    output {
        Boolean passed_qc = ReturnBoolAndString.passed_qc
        String errorMessages = ReturnBoolAndString.errorMessages
    }
}


task ReturnBoolAndString {
    String ubuntu_docker = "ubuntu:20.04"

    command {
        touch empty_file
    }

    runtime {
        docker: ubuntu_docker
        disk: "10 GB"
        memory: "1000 MiB"
        cpu: 1
        maxRetries: 2
    }
    output {
        Boolean passed_qc = true
        String errorMessages = ""
    }
}
