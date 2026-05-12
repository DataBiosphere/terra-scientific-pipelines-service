version 1.0

# Given a vcf, create a multisample vcf containing that sample duplicated vcf_duplicate_count number of times
workflow DuplicateVcfAndMerge {
    input{
        File vcf_file
        File vcf_file_index
        Int vcf_duplicate_count
        String base_name
        Int merge_memory_mb
        Int merge_disk_space
        Int merge_num_threads
    }

    call generateFofn {
        input:
            file_to_duplicate = vcf_file,
            duplicate_number = vcf_duplicate_count,
            base_name = base_name
    }

    call mergeFofn {
        input:
            vcf_files = read_lines(generateFofn.fofn),
            vcf_index = vcf_file_index,
            base_name = base_name,
            memory_mb = merge_memory_mb,
            disk_space = merge_disk_space,
            num_threads = merge_num_threads
    }

    output {
        File merged_vcf = mergeFofn.merged_vcf
        File merged_vcf_index = mergeFofn.merged_vcf_idx
    }
}

task generateFofn {
    input {
        String file_to_duplicate
        String base_name
        Int duplicate_number
    }

    command <<<
        set -e
        python3 <<CODE
        with open("~{base_name}.fofn", "w") as f:
        f.write(("~{file_to_duplicate}" + "\n") * ~{duplicate_number})
        CODE
    >>>

    runtime {
        docker: "python:3.12"
        memory: "2000 MiB"
        cpu: 1
        disks: "local-disk 20 HDD"
    }

    output {
        File fofn = "~{base_name}.fofn"
    }
}

task mergeFofn {
    input {
        Array[File] vcf_files
        File vcf_index
        String base_name
        Int memory_mb
        Int disk_space
        Int num_threads
    }


    command <<<
        bcftools merge -m all --force-samples --threads ~{num_threads} ~{sep=" " vcf_files} -Oz -o "~{base_name}.merged.vcf.gz"
        bcftools index -t "~{base_name}.merged.vcf.gz"
    >>>

    output {
        File merged_vcf = "${base_name}.merged.vcf.gz"
        File merged_vcf_idx = "${base_name}.merged.vcf.gz.tbi"
    }

    runtime {
        docker: "biocontainers/bcftools:v1.9-1-deb_cv1"
        memory: "${memory_mb} MiB"
        cpu: "${num_threads}"
        disks: "local-disk ${disk_space} HDD"
    }
}
