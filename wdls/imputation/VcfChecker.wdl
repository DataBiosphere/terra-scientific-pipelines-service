version 1.0

workflow VerifyGvcf {

    input {
        File test_gvcf
        File test_gvcf_index
        File truth_gvcf
        File truth_gvcf_index

        Boolean? done
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
        String patternForLinesToExcludeFromComparison = ""
    }
    Int disk_size_gb = ceil(size(file1, "GiB")) + ceil(size(file2, "GiB")) + 20
    command {
        set -eo pipefail

        if [ -z ~{patternForLinesToExcludeFromComparison} ]; then
        diff <(gunzip -c -f ~{file1}) <(gunzip -c -f ~{file2})
        else
        echo "It's defined!"
        diff <(gunzip -c -f ~{file1} | grep -v '~{patternForLinesToExcludeFromComparison}') <(gunzip -c -f ~{file2} | grep -v '~{patternForLinesToExcludeFromComparison}')
        fi
    }

    runtime {
        docker: "gcr.io/gcp-runtimes/ubuntu_16_0_4:latest"
        disks: "local-disk ${disk_size_gb} SSD"
        memory: "32 GiB"
    }
}
