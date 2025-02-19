version 1.0

workflow CompareVcfs {

    input {
        File test_gvcf
        File truth_gvcf

        Int memory_gb = 64
    }

    call CompareVcfsAllowingQualityDifferences {
        input:
            file1 = test_gvcf,
            file2 = truth_gvcf
    }
}

task CompareVcfsAllowingQualityDifferences {
    input {
        File file1
        File file2

        Int memory_gb = 64
    }

    Int disk_size_gb = ceil(10 * size(file1, "GiB")) + ceil(10 * size(file2, "GiB")) + 50

    command {
        exit_code=0

        cmp <(gunzip -c -f ~{file1} | grep -v '^##') <(gunzip -c -f ~{file2} | grep -v '^##')
        if [ $? -ne 0 ]; then
        exit_code=1
        echo "Error: VCF ~{file1} differs in content from ~{file2}" >&2
        cmp <(gunzip -c -f ~{file1} | grep -v '^##' | cut -f 1-5,7-) <(gunzip -c -f ~{file2} | grep -v '^##' | cut -f 1-5,7-)
        if [ $? -eq 0 ]; then
        echo "However they ONLY differ in the quality column" >&2
        fi
        fi

        exit $exit_code
    }

    runtime {
        docker: "gcr.io/gcp-runtimes/ubuntu_16_0_4@sha256:025124e2f1cf4d29149958f17270596bffe13fc6acca6252977c572dd5ba01bf"
        disks: "local-disk ${disk_size_gb} SSD"
        memory: "${memory_gb} GiB"
        preemptible: 3
    }
}
