version 1.0

workflow CompareBcfs {

    input {
        File test_bcf
        File truth_bcf

        Int memory_gb = 64
        Int preemptible_tries = 0
    }

    call CompareBcfsTask {
        input:
            file1 = test_bcf,
            file2 = truth_bcf,
            patternForLinesToExcludeFromComparison = "##", # ignore headers but do compare sample names,
            memory_gb = memory_gb,
            preemptible_tries = preemptible_tries
    }
}

task CompareBcfsTask {
    input {
        File file1
        File file2

        String patternForLinesToExcludeFromComparison = ""

        Int memory_gb = 64
        Int preemptible_tries = 0
    }

    Int disk_size_gb = ceil(10 * size(file1, "GiB")) + ceil(10 * size(file2, "GiB")) + 50

    command {
        set -eo pipefail

        if [ -z "~{patternForLinesToExcludeFromComparison}" ]; then
            diff <(bcftools view -O v ~{file1}) <(bcftools view -O v ~{file2})
        else
            echo "patternForLinesToExcludeFromComparison is defined: '~{patternForLinesToExcludeFromComparison}'"
            diff <(bcftools view -O v ~{file1} | grep -v '~{patternForLinesToExcludeFromComparison}') <(bcftools view -O v ~{file2} | grep -v '~{patternForLinesToExcludeFromComparison}')
        fi
    }

    runtime {
        docker: "us.gcr.io/broad-gotc-prod/bcftools-vcftools:sps_sv_docker_images"
        disks: "local-disk ${disk_size_gb} SSD"
        memory: "${memory_gb} GiB"
        preemptible: preemptible_tries
    }
}
