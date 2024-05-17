version 1.0

workflow VerifyGvcf {

    input {
        File test_gvcf
        File truth_gvcf

        Boolean? done
    }

    call CompareVcfs {
        input:
            file1 = test_gvcf,
            file2 = truth_gvcf,
            patternForLinesToExcludeFromComparison = "^##"
    }

    output {
        File vcf_differences = CompareVcfs.comm_output
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

        gunzip -c -f ~{file1} | grep -v '~{patternForLinesToExcludeFromComparison}' > file_1.vcf
        gunzip -c -f ~{file2} | grep -v '~{patternForLinesToExcludeFromComparison}' > file_2.vcf

        comm file_1.vcf file_2.vcf > output.txt

        if [ -s output.txt ]; then
        # The output is not-empty so that means there was a diff
        exit 1
        fi
    }

    runtime {
        docker: "gcr.io/gcp-runtimes/ubuntu_16_0_4:latest"
        disks: "local-disk ${disk_size_gb} SSD"
        memory: "32 GiB"
    }

    output {
        File comm_output = "output.txt"
    }
}
