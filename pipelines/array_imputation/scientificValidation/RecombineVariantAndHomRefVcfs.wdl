version 1.0

# This script is under review. It is not actively tested or maintained at this time.
workflow RecombineVariantAndHomRefVcfs {
    input {
        File variant_vcf
        File variant_vcf_index
        File hom_ref_vcf
        File hom_ref_vcf_index
        File ref_dict
        Array[String] contigs
        String output_basename

    }


    call RecombineVariantAndHomRefVcfs {
        input:
            variant_vcf = variant_vcf,
            variant_vcf_index = variant_vcf_index,
            hom_ref_vcf = hom_ref_vcf,
            hom_ref_vcf_index = hom_ref_vcf_index
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
    }

    command {
        set -e -o pipefail

        cat <<'EOF' > script.py
#!/usr/bin/env python3

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

def main():
    parser = argparse.ArgumentParser(description='number of samples to add default values for')
    parser.add_argument('-n', '--num-columns', type=int, default=10,
        help='number of samples to add default values for (default: 10)')

    args = parser.parse_args()

    add_columns_to_tsv(num_columns=args.num_columns)

if __name__ == "__main__":
    main()

        EOF

        ln -sf ~{variant_vcf} input.imputed_variants.vcf.gz
        ln -sf ~{variant_vcf_index} input.imputated_variants.vcf.gz.tbi
        ln -sf ~{hom_ref_vcf} input.hom_ref_sites_only.vcf.gz
        ln -sf ~{hom_ref_vcf_index} input.hom_ref_sites_only.vcf.gz.tbi

        bcftools view -h input.imputed_variants.vcf.gz > header.txt

        sample_count=$(bcftools query -l input.imputed_variants.vcf.gz | wc -l)

        bcftools reheader -h header.txt -o reheadered_sites_only.vcf.gz input.hom_ref_sites_only.vcf.gz

        # requires python 3
        gunzip -c reheadered_sites_only.vcf.gz | python3 expand_sites_only_vcf.py -n $sample_count | bgzip -c > reheadered_sites_only_expanded.vcf.gz

        tabix reheadered_sites_only_expanded.vcf.gz

        # requires gatk jar
        java -jar gatk.jar MergeVcfs -I input.imputed_variants.vcf.gz -I reheadered_sites_only_expanded.vcf.gz  -O ~{output_basename}.vcf.gz
        bcftools index -t ~{output_basename}.vcf.gz
    }

    output {
        File recombined_vcf = "~{output_basename}.vcf.gz"
        File recombined_vcf_index = "~{output_basename}.vcf.gz.tbi"
    }

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/ubuntu:20.04"
    }
}
