version 1.0

workflow GatkConcordanceValidation {
    input {
        Array[String] chromosomes = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10",
                                    "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18",
                                    "chr19", "chr20", "chr21", "chr22"]

        File eval_vcf
        File eval_vcf_index
        File truth_vcf
        File truth_vcf_index
        File af_annotation_vcf
        File af_annotation_vcf_index
        File sample_to_ancestry_af_annotation

        String output_basename

        File? sample_mapping
        Int? n_bins
        Float? right_edge_first_bin

        Boolean run_full_genome = true
        Boolean run_chromosomes = true

        Int preemptible = 0
    }

    if (run_chromosomes){
        scatter(chr in chromosomes) {
            call PearsonCorrelationByAF as PearsonByAF_chr {
                input:
                    evalVcf = eval_vcf,
                    evalVcfIndex = eval_vcf_index,
                    af_resource = af_annotation_vcf,
                    af_resource_index = af_annotation_vcf_index,
                    sample_to_ancestry_af_annotation = sample_to_ancestry_af_annotation,
                    truthVcf = truth_vcf,
                    truthVcfIndex = truth_vcf_index,
                    intervals = chr,
                    output_basename = output_basename + "." + chr,
                    n_bins = n_bins,
                    sample_mapping = sample_mapping,
                    right_edge_first_bin = right_edge_first_bin,
                    preemptible = preemptible
            }

            call AddConstantColumn as AddConstantColumn_chr {
                input:
                    input_tsv = PearsonByAF_chr.correlations,
                    constant_value = chr,
                    column_name = "CHROMOSOME",
                    output_filename = output_basename + "." + chr + ".correlations"
            }
        }
        call ConcatenateTsvs as ConcatenateChrSpecificCorrelationTsvs {
            input:
                input_tsvs = AddConstantColumn_chr.output_file,
                output_filename = output_basename + ".chr_specific_correlations",
                preemptible = preemptible
        }
    }

    if (run_full_genome) {
        call PearsonCorrelationByAF as PearsonByAF_WholeGenome {
            input:
                evalVcf = eval_vcf,
                evalVcfIndex = eval_vcf_index,
                af_resource = af_annotation_vcf,
                af_resource_index = af_annotation_vcf_index,
                sample_to_ancestry_af_annotation = sample_to_ancestry_af_annotation,
                truthVcf = truth_vcf,
                truthVcfIndex = truth_vcf_index,
                output_basename = output_basename,
                n_bins = n_bins,
                sample_mapping = sample_mapping,
                right_edge_first_bin = right_edge_first_bin,
                preemptible = preemptible
        }

        call AddConstantColumn as AddConstantColumn_whole_genome {
            input:
                input_tsv = PearsonByAF_WholeGenome.correlations,
                constant_value = "WholeGenome",
                column_name = "CHROMOSOME",
                output_filename = output_basename + ".whole_genome_correlations"
        }
    }



output {
    File? correlations_chr_concatenated = ConcatenateChrSpecificCorrelationTsvs.output_file
    File? correlations_whole_genome = AddConstantColumn_whole_genome.output_file
  }
}


task PearsonCorrelationByAF {
    input {
        File evalVcf
        File evalVcfIndex
        File truthVcf
        File truthVcfIndex
        String output_basename
        File af_resource
        File af_resource_index
        File sample_to_ancestry_af_annotation
        File? sample_mapping
        File? sites
        String? intervals
        String? dosage_field
        Int? n_bins
        Float? right_edge_first_bin
        Float? min_af_for_accuracy_metrics
        Int? n_calibration_bins
        Int mem_gb = 6
        Int preemptible = 3
    }

    Int disk_size = ceil(size([evalVcf, truthVcf, af_resource, sample_to_ancestry_af_annotation], "GB") +  20)

    command <<<
        set -xeuo pipefail

        gatk --java-options "-Xmx~{mem_gb - 2}G" EvaluateGenotypingPerformance \
        -eval ~{evalVcf} \
        -truth ~{truthVcf} \
        --arguments_file ~{sample_to_ancestry_af_annotation} \
        --resource ~{af_resource} \
        ~{"--ids " + sites} \
        ~{"-L " + intervals} \
        ~{"--dosage-field " + dosage_field} \
        ~{"--sample-map " + sample_mapping} \
        -O ~{output_basename}.correlations.tsv \
        -OA ~{output_basename}.accuracy.tsv \
        ~{"-nbins " + n_bins} \
        ~{"-first-bin-right-edge " + right_edge_first_bin} \
        ~{"--min-af-for-accuracy-metrics " + min_af_for_accuracy_metrics} \
        --allow-differing-ploidies \
        --output-gp-calibration ~{output_basename}.gp_calibration.tsv \
        ~{"--n-calibration-bins " + n_calibration_bins} \
        --output-accuracy-af ~{output_basename}.accuracy_af.tsv
    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/ckachulis/gatk-array-correlation@sha256:c5eb54fdc4a9dabe4a6dda25af1ad1fe1f10f93c91bd0653ec2a49e4253c1f2e"
        disks: "local-disk " + disk_size +" SSD"
        memory: mem_gb + " GB"
        preemptible: preemptible
    }

    output {
        File correlations = "~{output_basename}.correlations.tsv"
        File accuracy = "~{output_basename}.accuracy.tsv"
        File accuracy_af = "~{output_basename}.accuracy_af.tsv"
        File gp_calibration = "~{output_basename}.gp_calibration.tsv"
    }
}

task AddConstantColumn {
    input {
        File input_tsv
        String constant_value
        String column_name
        String output_filename = "output"

        Int disk_size = ceil(2 * size(input_tsv, "GB") + 10)
        Int mem_gb = 2
        Int preemptible = 3
    }

    String output_tsv = output_filename + ".tsv"

    command <<<
        set -xeuo pipefail

        python <<CODE
        import pandas as pd

        # Read the TSV file
        df = pd.read_csv("~{input_tsv}", sep='\t', comment='#')

        # Add the new column with the constant value
        df["~{column_name}"] = "~{constant_value}"

        # Write the result to a new TSV file
        df.to_csv("~{output_tsv}", sep='\t', index=False)
        CODE
    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/python-data-slim"
        memory: mem_gb + " GB"
        cpu: 1
        disks: "local-disk " + disk_size + " HDD"
        preemptible: preemptible
    }

    output {
        File output_file = output_tsv
    }
}

task ConcatenateTsvs {
    input {
        Array[File] input_tsvs
        String output_filename = "concatenated"

        Int disk_size = ceil(2 * size(input_tsvs, "GB") + 10)
        Int mem_gb = 2
        Int preemptible = 3
    }

    String output_tsv = output_filename + ".tsv"

    command <<<
        set -xeuo pipefail

        python <<CODE
        import pandas as pd

        # Create and concatenate dataframes in one step with a list comprehension
        input_files = ['~{sep="','" input_tsvs}']
        combined_df = pd.concat([pd.read_csv(file, sep='\t', comment='#') for file in input_files])

        # Write the concatenated dataframe to a new TSV file
        combined_df.to_csv("~{output_tsv}", sep='\t', index=False)
        CODE
    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/python-data-slim"
        memory: mem_gb + " GB"
        cpu: 1
        disks: "local-disk " + disk_size + " HDD"
        preemptible: preemptible
    }

    output {
        File output_file = output_tsv
    }
}
