version 1.0

workflow CompareVcfs {

    input {
        File test_gvcf
        File truth_gvcf

        Int memory_gb = 64
    }

    call CompareVcfs {
        input:
            file1 = test_gvcf,
            file2 = truth_gvcf,
            patternForLinesToExcludeFromComparison = "##" # ignore headers but do compare sample names
    }
}

task CompareVcfsAllowingQualityDifferences {
    input {
        File file1
        File file2

        String patternForLinesToExcludeFromComparison = ""

        Int memory_gb = 64
    }

    Int disk_size_gb = ceil(10 * size(file1, "GiB")) + ceil(10 * size(file2, "GiB")) + 50

    command {
        set -eo pipefail

        if [ -z ~{patternForLinesToExcludeFromComparison} ]; then
          diff <(gunzip -c -f ~{file1}) <(gunzip -c -f ~{file2})
        else
          echo "patternForLinesToExcludeFromComparison is defined: ~{patternForLinesToExcludeFromComparison}"
          diff <(gunzip -c -f ~{file1} | grep -v '~{patternForLinesToExcludeFromComparison}') <(gunzip -c -f ~{file2} | grep -v '~{patternForLinesToExcludeFromComparison}')
        fi
    }

    runtime {
        docker: "gcr.io/gcp-runtimes/ubuntu_16_0_4@sha256:025124e2f1cf4d29149958f17270596bffe13fc6acca6252977c572dd5ba01bf"
        disks: "local-disk ${disk_size_gb} SSD"
        memory: "${memory_gb} GiB"
        preemptible: 3
    }
}
