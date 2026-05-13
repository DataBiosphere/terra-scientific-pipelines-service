version 1.0

# This script is under review. It is not actively tested or maintained at this time.
workflow RecombineVariantAndHomRefVcfs {
    input {
        File variant_vcf
        File variant_vcf_index
        File hom_ref_vcf
        File hom_ref_vcf_index
        String output_basename
    }

    call RecombineVariantAndHomRefVcfs {
        input:
            variant_vcf = variant_vcf,
            variant_vcf_index = variant_vcf_index,
            hom_ref_vcf = hom_ref_vcf,
            hom_ref_vcf_index = hom_ref_vcf_index,
            output_basename = output_basename
    }

    output {
        File recombined_vcf = RecombineVariantAndHomRefVcfs.recombined_vcf
        File recombined_vcf_index = RecombineVariantAndHomRefVcfs.recombined_vcf_index
    }
}

task RecombineVariantAndHomRefVcfs {
    input {
        File variant_vcf
        File variant_vcf_index
        File hom_ref_vcf
        File hom_ref_vcf_index
        String output_basename = "merged.all_variants"
        Int memory_mb = 4000
        String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.6.1.0"
    }

    Int disk_size = ceil(2*(size(variant_vcf, "GiB") + size(hom_ref_vcf, "GiB")) + 20)
    Int command_mem = memory_mb - 1500
    Int max_heap = memory_mb - 1000

    command {
        set -e -o pipefail

        cat <<'EOF' > expand_sites_only_vcf.py
        import sys
        import argparse

        def add_columns_to_tsv(fill_value="0|0:0", num_columns=10):
            """
            Add specified number of columns to each row in a TSV while preserving header lines starting with '#'
            Reads from stdin and writes to stdout
            """
            # Create the additional columns string
            additional_columns = '\t'.join([fill_value] * num_columns)

            for line in sys.stdin:
                line = line.rstrip('\n')

                if line.startswith('#'):
                    # Pass through header lines unchanged
                    print(line)
                else:
                    print(line + '\t' + 'GT:DS' + '\t' + additional_columns)

        parser = argparse.ArgumentParser(description='number of samples to add default values for')
        parser.add_argument('-n', '--num-columns', type=int, default=10,
            help='number of samples to add default values for (default: 10)')

        args = parser.parse_args()

        add_columns_to_tsv(num_columns=args.num_columns)
        EOF

        ln -sf ~{variant_vcf} input.imputed_variants.vcf.gz
        ln -sf ~{variant_vcf_index} input.imputated_variants.vcf.gz.tbi
        ln -sf ~{hom_ref_vcf} input.hom_ref_sites_only.vcf.gz
        ln -sf ~{hom_ref_vcf_index} input.hom_ref_sites_only.vcf.gz.tbi

        echo "grabbing header from input imputed variants vcf"
        bcftools view -h input.imputed_variants.vcf.gz > header.txt

        echo "finding sample count"
        sample_count=$(bcftools query -l input.imputed_variants.vcf.gz | wc -l)
        echo "found sample count: $sample_count"

        echo "reheadering hom ref sites only vcf to have same header as input imputed variants vcf"
        bcftools reheader -h header.txt -o reheadered_sites_only.vcf.gz input.hom_ref_sites_only.vcf.gz

        # requires python 3
        echo "running python script to add GT:DS columns with default values to hom ref sites only vcf"
        gunzip -c reheadered_sites_only.vcf.gz | python3 expand_sites_only_vcf.py -n $sample_count | bgzip -c > reheadered_sites_only_expanded.vcf.gz

        echo "indexing reheadered and expanded hom ref sites only vcf"
        tabix reheadered_sites_only_expanded.vcf.gz

        # requires gatk jar
        echo "merging input imputed variants vcf and reheadered and expanded hom ref sites only vcf"
        gatk --java-options "-Xms~{command_mem}m -Xmx~{max_heap}m" \
        MergeVcfs \
        -I input.imputed_variants.vcf.gz \
        -I reheadered_sites_only_expanded.vcf.gz \
        -O ~{output_basename}.vcf.gz
    }

    output {
        File recombined_vcf = "~{output_basename}.vcf.gz"
        File recombined_vcf_index = "~{output_basename}.vcf.gz.tbi"
    }

    runtime {
        docker: gatk_docker
        preemptible: 0
        retries: 1
        memory: "${memory_mb} MiB"
        cpu: 2
        disks: "local-disk ${disk_size} HDD"
    }
}
