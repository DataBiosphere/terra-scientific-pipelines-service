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
        Int batch_size = 50
        Int? n_bins
        Float? right_edge_first_bin

        Boolean run_full_genome = true
        Boolean run_chromosomes = true

        Int preemptible = 0
    }

    call GetSamplesFromVcf as EvalVcfSamples {
        input:
            vcf = eval_vcf
    }

    Float batch_size_float = batch_size
    Int num_sample_batches = ceil(EvalVcfSamples.sample_count / batch_size_float)

    scatter (batch_index in range(num_sample_batches)) {
        Int batch_start_field = (batch_index * batch_size) + 10
        Int batch_end_field = if (((batch_index + 1) * batch_size) >= EvalVcfSamples.sample_count) then EvalVcfSamples.sample_count + 9 else (((batch_index + 1) * batch_size) + 9)

        if (num_sample_batches > 1) {
            call SelectSamplesWithCut as SelectEvalVcfBatch {
                input:
                    vcf = eval_vcf,
                    cut_start_field = batch_start_field,
                    cut_end_field = batch_end_field,
                    basename = output_basename + ".sample_batch_" + batch_index
            }
        }

        File batch_eval_vcf = select_first([SelectEvalVcfBatch.output_vcf, eval_vcf])
        File batch_eval_vcf_index = select_first([SelectEvalVcfBatch.output_vcf_index, eval_vcf_index])

        if (run_chromosomes) {
            scatter(chr in chromosomes) {
                call PearsonCorrelationByAF as PearsonByAF_chr_batch {
                    input:
                        evalVcf = batch_eval_vcf,
                        evalVcfIndex = batch_eval_vcf_index,
                        af_resource = af_annotation_vcf,
                        af_resource_index = af_annotation_vcf_index,
                        sample_to_ancestry_af_annotation = sample_to_ancestry_af_annotation,
                        truthVcf = truth_vcf,
                        truthVcfIndex = truth_vcf_index,
                        intervals = chr,
                        output_basename = output_basename + "." + chr + ".sample_batch_" + batch_index,
                        n_bins = n_bins,
                        sample_mapping = sample_mapping,
                        right_edge_first_bin = right_edge_first_bin,
                        preemptible = preemptible
                }

                call AddConstantColumn as AddConstantColumn_chr_batch {
                    input:
                        input_tsv = PearsonByAF_chr_batch.correlations,
                        constant_value = chr,
                        column_name = "CHROMOSOME",
                        output_filename = output_basename + "." + chr + ".sample_batch_" + batch_index + ".correlations"
                }
            }
            call ConcatenateTsvs as ConcatenateChrSpecificCorrelationTsvs_batch {
                input:
                    input_tsvs = AddConstantColumn_chr_batch.output_file,
                    output_filename = output_basename + ".chr_specific_correlations.sample_batch_" + batch_index,
                    preemptible = preemptible
            }
        }

        if (run_full_genome) {
            call PearsonCorrelationByAF as PearsonByAF_WholeGenome_batch {
                input:
                    evalVcf = batch_eval_vcf,
                    evalVcfIndex = batch_eval_vcf_index,
                    af_resource = af_annotation_vcf,
                    af_resource_index = af_annotation_vcf_index,
                    sample_to_ancestry_af_annotation = sample_to_ancestry_af_annotation,
                    truthVcf = truth_vcf,
                    truthVcfIndex = truth_vcf_index,
                    output_basename = output_basename + ".sample_batch_" + batch_index,
                    n_bins = n_bins,
                    sample_mapping = sample_mapping,
                    right_edge_first_bin = right_edge_first_bin,
                    preemptible = preemptible
            }

            call AddConstantColumn as AddConstantColumn_whole_genome_batch {
                input:
                    input_tsv = PearsonByAF_WholeGenome_batch.correlations,
                    constant_value = "WholeGenome",
                    column_name = "CHROMOSOME",
                    output_filename = output_basename + ".sample_batch_" + batch_index + ".whole_genome_correlations"
            }
        }
    }

    if (run_chromosomes) {
        call ConcatenateTsvs as ConcatenateChrSpecificCorrelationTsvs_sampled {
            input:
                input_tsvs = select_all(ConcatenateChrSpecificCorrelationTsvs_batch.output_file),
                output_filename = output_basename + ".chr_specific_correlations",
                preemptible = preemptible
        }
    }

    if (run_full_genome) {
        call ConcatenateTsvs as ConcatenateWholeGenomeCorrelations_sampled {
            input:
                input_tsvs = select_all(AddConstantColumn_whole_genome_batch.output_file),
                output_filename = output_basename + ".whole_genome_correlations",
                preemptible = preemptible
        }
    }

    output {
        File? correlations_chr_concatenated = ConcatenateChrSpecificCorrelationTsvs_sampled.output_file
        File? correlations_whole_genome = ConcatenateWholeGenomeCorrelations_sampled.output_file
    }
}

task GetSamplesFromVcf {
    input {
        File vcf
    }

    Int disk_size = ceil(size(vcf, "GB") + 10)

    command <<<
        set -euo pipefail

        bcftools query -l "~{vcf}" > eval_vcf_samples.txt
        wc -l < eval_vcf_samples.txt > eval_vcf_samples_count.txt
    >>>

    runtime {
        docker: "us.gcr.io/broad-gotc-prod/bcftools-vcftools:2.0.0-1.24-0.1.17-1784569943"
        memory: "2 GB"
        cpu: 1
        disks: "local-disk " + disk_size + " SSD"
        preemptible: 0
        maxRetries: 1
    }

    output {
        File sample_names_file = "eval_vcf_samples.txt"
        Int sample_count = read_int("eval_vcf_samples_count.txt")
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
        Int preemptible = 0
    }

    Int disk_size = ceil((2 * size(evalVcf, "GB")) + size([truthVcf, af_resource, sample_to_ancestry_af_annotation], "GB") + 20)

    command <<<
        set -xeuo pipefail

        gatk --java-options "-Xmx~{mem_gb - 2}G" EvaluateGenotypingPerformance \
        -eval "~{evalVcf}" \
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
        maxRetries: 1
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
        maxRetries: 1
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
        maxRetries: 1
    }

    output {
        File output_file = output_tsv
    }
}

task SelectSamplesWithCut {
    input {
        File vcf

        Int cut_start_field
        Int cut_end_field
        String basename

        Int disk_size_gb = ceil(1.5 * size(vcf, "GiB")) + 10
        String bcftools_docker = "us.gcr.io/broad-gotc-prod/imputation-bcf-vcf:1.0.7-1.10.2-0.1.16-1669908889"
        Int cpu = 2
        Int memory_mb = 6000
    }

    command <<<
        set -euo pipefail

        mkfifo fifo_bgzip
        mkfifo fifo_cut

        bcftools view -h --no-version ~{vcf} | awk '!/^#CHROM/' > header.vcf
        n_lines=$(wc -l header.vcf | cut -d' ' -f1)

        cat header.vcf
        echo $n_lines

        bgzip -d ~{vcf} -o fifo_bgzip &
        tail +$((n_lines)) fifo_bgzip | cut -f 1-9,~{cut_start_field}-~{cut_end_field} > fifo_cut &

        cat header.vcf fifo_cut | bgzip -o ~{basename}.vcf.gz
        bcftools index -t ~{basename}.vcf.gz
    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/bcftools_bgzip:beagle_imputation_v1.0.0"
        disks: "local-disk " + disk_size_gb + " SSD"
        memory: memory_mb + " MiB"
        preemptible: 0
        maxRetries: 1
        noAddress: true
        cpu: cpu
    }

    output {
        File output_vcf = "~{basename}.vcf.gz"
        File output_vcf_index = "~{basename}.vcf.gz.tbi"
    }
}
