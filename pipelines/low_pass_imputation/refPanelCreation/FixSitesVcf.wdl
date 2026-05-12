version 1.0

workflow FixSitesVcf {
    input {
        File input_vcf
        String output_filename = "output"
    }

    call CutVcfColumns {
        input:
            input_vcf = input_vcf,
            output_filename = output_filename
    }

    output {
        File output_vcf = CutVcfColumns.output_vcf
        File output_vcf_index = CutVcfColumns.output_vcf_index

        File output_table = CutVcfColumns.output_table
        File output_table_index = CutVcfColumns.output_table_index
    }
}

task CutVcfColumns {
    input {
        File input_vcf
        String output_filename
        Int disk_size = ceil(2.5 * size(input_vcf, "GiB") + 20)
    }

    command <<<
        set -euo pipefail

        # Decompress, cut columns 1-8, and recompress in one pipeline
        gunzip -c ~{input_vcf} | cut -f 1-8 | bgzip -c > cut.vcf.gz
        bcftools index -t cut.vcf.gz

        # Trim alleles at sites with more than 30 ALTs, taking top most frequent alleles
        python3 << CODE
        import pysam

        with pysam.VariantFile('cut.vcf.gz') as in_vcf:
            with pysam.VariantFile('with-annotations.vcf.gz', 'w', header=in_vcf.header) as out_vcf:
            for rec in in_vcf:
                if len(rec.alleles) > 31:    # 31 = 1 REF + 30 ALTs
                    sorted_alts = sorted(rec.alleles[1:], key=lambda x: dict(zip(rec.alleles[1:], rec.info['AF']))[x], reverse=True)
                    rec.alleles = tuple([rec.alleles[0]] + sorted_alts[:30])
                out_vcf.write(rec)
        CODE

        # Drop allele-specific INFO fields to avoid future confusion
        bcftools annotate -x 'INFO/AC,INFO/AF,INFO/AN,INFO/HC' with-annotations.vcf.gz -o ~{output_filename}.vcf.gz

        # Index the VCF using bcftools instead of tabix
        bcftools index -t ~{output_filename}.vcf.gz

        # Make sites table for bcftools call to use
        bcftools query -f'%CHROM\t%POS\t%REF,%ALT\n' ~{output_filename}.vcf.gz | bgzip -c > ~{output_filename}.tsv.gz && tabix -s1 -b2 -e2 ~{output_filename}.tsv.gz
    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/vcfeval_pysam:v1.0"
        memory: "4 GiB"
        cpu: 2
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 3
    }

    output {
        File output_vcf = "~{output_filename}.vcf.gz"
        File output_vcf_index = "~{output_filename}.vcf.gz.tbi"

        File output_table = "~{output_filename}.tsv.gz"
        File output_table_index = "~{output_filename}.tsv.gz.tbi"
    }
}
