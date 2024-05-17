version 1.0

workflow VerifyGvcf {

    input {
        File test_gvcf
        File truth_gvcf
    }

    call CompareVcfs {
        input:
            file1 = test_gvcf,
            file2 = truth_gvcf,
            patternForLinesToExcludeFromComparison = "^##"
    }
}

task CompareVcfs {
    input {
        File file1
        File file2
        String patternForLinesToExcludeFromComparison
    }
    Int disk_size_gb = ceil(3 * size(file1, "GiB")) + ceil(3 * size(file2, "GiB")) + 50
    command {
        set -eo pipefail

        comm --nocheck-order <(gunzip -c -f ~{file1} | grep -v '~{patternForLinesToExcludeFromComparison}')  <(gunzip -c -f ~{file2} | grep -v '~{patternForLinesToExcludeFromComparison}') > comm_output.txt

        if [ -s comm_output.txt ]; then
        # file is not empty so that means comm found differences
        cat comm_output.txt
        exit 1
        fi
    }

    runtime {
        docker: "gcr.io/gcp-runtimes/ubuntu_16_0_4:latest"
        disks: "local-disk ${disk_size_gb} SSD"
        memory: "8 GiB"
    }
}
