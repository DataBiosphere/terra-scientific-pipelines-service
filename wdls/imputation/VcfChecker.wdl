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
        String patternForLinesToExcludeFromComparison
    }
    Int disk_size_gb = ceil(3 * size(file1, "GiB")) + ceil(3 * size(file2, "GiB")) + 50
    command {
        set -eo pipefail

        gunzip -c -f ~{file1} | grep -v '~{patternForLinesToExcludeFromComparison}' > file_1.vcf
        gunzip -c -f ~{file2} | grep -v '~{patternForLinesToExcludeFromComparison}' > file_2.vcf

        diff file_1.vcf file_2.vcf
    }

    runtime {
        docker: "gcr.io/gcp-runtimes/ubuntu_16_0_4:latest"
        disks: "local-disk ${disk_size_gb} SSD"
        memory: "32 GiB"
    }
}
