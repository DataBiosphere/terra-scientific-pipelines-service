version 1.0

workflow CompareVcfs {

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
    Int disk_size_gb = ceil(10 * size(file1, "GiB")) + ceil(10 * size(file2, "GiB")) + 50
    command {
        set -eo pipefail

        gunzip -c -f ~{file1} | grep -v '~{patternForLinesToExcludeFromComparison}' > file_1.vcf
        gunzip -c -f ~{file2} | grep -v '~{patternForLinesToExcludeFromComparison}' > file_2.vcf

        diff -q file_1.vcf file_2.vcf
    }

    runtime {
        docker: "gcr.io/gcp-runtimes/ubuntu_16_0_4:latest"
        disks: "local-disk ${disk_size_gb} SSD"
        memory: "64 GiB"
    }
}
