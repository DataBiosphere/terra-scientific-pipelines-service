version 1.0

workflow GLIMPSE2Concordance {
    input {
        Array[File] panel_vcfs
        Array[File] panel_vcf_idxs
        Array[File] imputed_vcfs
        Array[File] imputed_vcf_idxs
        File trh_bed
        File trh_bed_idx
        Array[String] regions
        String output_prefix
    }

    Array[String] trh_bins = ["outTRH", "inTRH"]
    Array[String] length_bins = ["SV_DEL", "DEL", "SNP", "INS", "SV_INS"]

    scatter (idx in range(length(regions))) {
        call AnnotateImputed { input:
                                   imputed_vcf = imputed_vcfs[idx],
                                   imputed_vcf_idx = imputed_vcf_idxs[idx],
                                   trh_bed = trh_bed,
                                   trh_bed_idx = trh_bed_idx,
                                   region = regions[idx],
                                   output_prefix = output_prefix + "." + regions[idx]
        }

        scatter (trh_bin in trh_bins) {
            scatter (length_bin in length_bins) {
                call FilterAndConcordance { input:
                                                annotated_bcf = AnnotateImputed.annotated_vcf,
                                                annotated_bcf_index = AnnotateImputed.annotated_vcf_idx,
                                                panel_vcf = panel_vcfs[idx],
                                                panel_vcf_idx = panel_vcf_idxs[idx],
                                                trh_bin = trh_bin,
                                                length_bin = length_bin,
                                                region = regions[idx],
                                                output_prefix = output_prefix + "." + regions[idx] + "." + trh_bin + "." + length_bin
                }
            }
        }

        # Generate Per-Chromosome Plots
        Array[File] chrom_rsquare_grp = flatten(flatten(FilterAndConcordance.rsquare_grp_files))
        Array[File] chrom_error_spl = flatten(flatten(FilterAndConcordance.error_spl_files))

        call PlotResults as PlotResultsPerChrom { input:
                                                      rsquare_grp_files = chrom_rsquare_grp,
                                                      error_spl_files = chrom_error_spl,
                                                      output_prefix = output_prefix + "." + regions[idx],
                                                      panel_name = basename(panel_vcfs[idx]),
                                                      imputed_name = basename(imputed_vcfs[idx])
        }
    }

    # Generate Aggregate Plots
    Array[File] all_rsquare_grp_files = flatten(flatten(flatten(FilterAndConcordance.rsquare_grp_files)))
    Array[File] all_rsquare_spl_files = flatten(flatten(flatten(FilterAndConcordance.rsquare_spl_files)))
    Array[File] all_error_grp_files = flatten(flatten(flatten(FilterAndConcordance.error_grp_files)))
    Array[File] all_error_spl_files = flatten(flatten(flatten(FilterAndConcordance.error_spl_files)))
    Array[File] all_error_cal_files = flatten(flatten(flatten(FilterAndConcordance.error_cal_files)))

    call PlotResults as PlotResultsAggregate { input:
                                                   rsquare_grp_files = all_rsquare_grp_files,
                                                   error_spl_files = all_error_spl_files,
                                                   output_prefix = output_prefix + ".aggregate",
                                                   panel_name = "Aggregate Panel",
                                                   imputed_name = "Aggregate Imputed"
    }

    output {
        Array[File] concordance_aggregate_plots_pdf = [PlotResultsAggregate.r2_plot_inTRH_pdf, PlotResultsAggregate.r2_plot_outTRH_pdf, PlotResultsAggregate.nrd_plot_inTRH_pdf, PlotResultsAggregate.nrd_plot_outTRH_pdf]
        Array[File] concordance_aggregate_plots_png = [PlotResultsAggregate.r2_plot_inTRH_png, PlotResultsAggregate.r2_plot_outTRH_png, PlotResultsAggregate.nrd_plot_inTRH_png, PlotResultsAggregate.nrd_plot_outTRH_png]

        Array[File] concordance_per_chrom_plots_pdf = flatten([PlotResultsPerChrom.r2_plot_inTRH_pdf, PlotResultsPerChrom.r2_plot_outTRH_pdf, PlotResultsPerChrom.nrd_plot_inTRH_pdf, PlotResultsPerChrom.nrd_plot_outTRH_pdf])
        Array[File] concordance_per_chrom_plots_png = flatten([PlotResultsPerChrom.r2_plot_inTRH_png, PlotResultsPerChrom.r2_plot_outTRH_png, PlotResultsPerChrom.nrd_plot_inTRH_png, PlotResultsPerChrom.nrd_plot_outTRH_png])

        Array[Array[File]] concordance_results = [all_rsquare_grp_files, all_rsquare_spl_files, all_error_grp_files, all_error_spl_files, all_error_cal_files]
    }
}

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    String? disk_type
    Int? preemptible_tries
    Int? max_retries
    String? docker
}

task AnnotateImputed {
    input {
        File imputed_vcf
        File imputed_vcf_idx
        File trh_bed
        File trh_bed_idx
        String region
        String output_prefix

        RuntimeAttr? runtime_attr_override
    }

    Int disk_gb = 10 + 3 * ceil(size(imputed_vcf, "GiB") + size(trh_bed, "GiB"))

    command <<<
        set -euox pipefail

        bcftools annotate ~{imputed_vcf} \
        --threads $(nproc) \
        -r ~{region} \
        -a ~{trh_bed} -c CHROM,FROM,TO -m +TRH \
        --write-index=csi -Ob -o ~{output_prefix}.imputed.annotated.bcf
    >>>

    output {
        File annotated_vcf = "~{output_prefix}.imputed.annotated.bcf"
        File annotated_vcf_idx = "~{output_prefix}.imputed.annotated.bcf.csi"
    }

    RuntimeAttr default_attr = object {
                                   cpu_cores:          2,
                                   mem_gb:             8,
                                   disk_gb:            disk_gb,
                                   boot_disk_gb:       10,
                                   disk_type:          "SSD",
                                   preemptible_tries:  2,
                                   max_retries:        1,
                                   docker:             "us.gcr.io/broad-dsp-lrma/lr-gcloud-samtools:0.1.23"
                               }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " " + select_first([runtime_attr.disk_type, default_attr.disk_type])
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task FilterAndConcordance {
    input {
        File annotated_bcf
        File annotated_bcf_index
        File panel_vcf
        File panel_vcf_idx
        String trh_bin
        String length_bin
        String region
        String output_prefix

        RuntimeAttr? runtime_attr_override
    }

    Int disk_gb = 10 + 2 * ceil(size(annotated_bcf, "GiB") + size(panel_vcf, "GiB"))

    command <<<
        set -euox pipefail

        # Determine TRH Expression
        if [ "~{trh_bin}" == "outTRH" ]; then
        TRH_EXP="TRH!=1"
        else
        TRH_EXP="TRH==1"
        fi

        # Determine Length Expression
        if [ "~{length_bin}" == "SV_DEL" ]; then LEN_EXP="(STRLEN(ALT)-STRLEN(REF) <= -50)"; fi
        if [ "~{length_bin}" == "DEL" ]; then LEN_EXP="((-50 < STRLEN(ALT)-STRLEN(REF)) && (STRLEN(ALT)-STRLEN(REF) <= -1))"; fi
        if [ "~{length_bin}" == "SNP" ]; then LEN_EXP="((STRLEN(REF) == 1) && (STRLEN(ALT) == 1))"; fi
        if [ "~{length_bin}" == "INS" ]; then LEN_EXP="((STRLEN(REF) != 1) && (0 <= STRLEN(ALT)-STRLEN(REF)) && (STRLEN(ALT)-STRLEN(REF) < 50))"; fi
        if [ "~{length_bin}" == "SV_INS" ]; then LEN_EXP="(50 <= STRLEN(ALT)-STRLEN(REF))"; fi

        echo "Filtering with: $TRH_EXP & $LEN_EXP"

        # 1. Unfiltered and GP>0.9 Evaluations
        bcftools view ~{annotated_bcf} \
        -i "$TRH_EXP & $LEN_EXP" \
        --threads $(nproc) \
        --write-index=csi -Ob -o ~{output_prefix}.bcf

        echo "~{region} ~{panel_vcf} ~{panel_vcf} ~{output_prefix}.bcf" > ~{output_prefix}.concordance-input.txt

        wget https://github.com/odelaneau/GLIMPSE/releases/download/v2.0.1/GLIMPSE2_concordance_static
        chmod +x GLIMPSE2_concordance_static

        ./GLIMPSE2_concordance_static \
        --min-tar-gp 0.0 0.9 \
        --gt-val \
        --use-alt-af \
        --out-r2-per-site \
        --bins 0.00001 0.00002 0.00005 0.0001 0.0002 0.0005 0.001 0.002 0.005 0.01 0.02 0.05 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.95 0.99 0.995 0.999 1.0 \
        --input ~{output_prefix}.concordance-input.txt \
        --threads $(nproc) \
        --output ~{output_prefix}.concordance-result

        # 2. INFO>0.5 Evaluations
        bcftools view ~{output_prefix}.bcf \
        -i 'INFO/INFO>0.5' \
        --threads $(nproc) \
        --write-index=csi -Ob -o ~{output_prefix}.INFO05.bcf

        echo "~{region} ~{panel_vcf} ~{panel_vcf} ~{output_prefix}.INFO05.bcf" > ~{output_prefix}.INFO05.concordance-input.txt

        ./GLIMPSE2_concordance_static \
        --min-tar-gp 0.0 \
        --gt-val \
        --use-alt-af \
        --out-r2-per-site \
        --bins 0.00001 0.00002 0.00005 0.0001 0.0002 0.0005 0.001 0.002 0.005 0.01 0.02 0.05 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.95 0.99 0.995 0.999 1.0 \
        --input ~{output_prefix}.INFO05.concordance-input.txt \
        --threads $(nproc) \
        --output ~{output_prefix}.INFO05.concordance-result
    >>>

    output {
        Array[File] rsquare_grp_files = glob("*.rsquare.grp.txt.gz")
        Array[File] rsquare_spl_files = glob("*.rsquare.spl.txt.gz")
        Array[File] error_grp_files = glob("*.error.grp.txt.gz")
        Array[File] error_spl_files = glob("*.error.spl.txt.gz")
        Array[File] error_cal_files = glob("*.error.cal.txt.gz")
    }

    RuntimeAttr default_attr = object {
                                   cpu_cores:          4,
                                   mem_gb:             32,
                                   disk_gb:            disk_gb,
                                   boot_disk_gb:       10,
                                   disk_type:          "SSD",
                                   preemptible_tries:  2,
                                   max_retries:        1,
                                   docker:             "us.gcr.io/broad-dsp-lrma/lr-gcloud-samtools:0.1.23"
                               }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " " + select_first([runtime_attr.disk_type, default_attr.disk_type])
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task PlotResults {
    input {
        Array[File] rsquare_grp_files
        Array[File] error_spl_files
        String output_prefix
        String panel_name
        String imputed_name

        RuntimeAttr? runtime_attr_override
    }

    Int disk_gb = 10 + ceil(size(rsquare_grp_files, "GiB") + size(error_spl_files, "GiB"))

    command <<<
        set -euox pipefail

        # Write file paths to list files to prevent "Argument list too long" errors
        cat ~{write_lines(rsquare_grp_files)} > rsquare_files.list
        cat ~{write_lines(error_spl_files)} > error_files.list

        cat << 'EOF' > plot_script.py
        import sys
        import os
        import numpy as np
        import pandas as pd
        import matplotlib.pyplot as plt
        import seaborn as sns

        # Read files from the written lists
        with open(sys.argv[1], 'r') as f:
        rsquare_files = [line.strip() for line in f if line.strip()]

        with open(sys.argv[2], 'r') as f:
        error_files = [line.strip() for line in f if line.strip()]

        panel_name = sys.argv[3]
        imputed_name = sys.argv[4]
        output_prefix = sys.argv[5]

        # 1. Load Dosage R2 Data
        r2_vs_af_df_values = []
        for filepath in rsquare_files:
        if not filepath.strip(): continue
        filename = os.path.basename(filepath)

        parts = filename.split('_GPfilt_')
        prefix_parts = parts[0].split('.')
        is_info05 = 'INFO05' in prefix_parts

        if is_info05:
        length_bin = prefix_parts[-3]
        trh_bin = prefix_parts[-4]
        label_type = 'INFO05'
        else:
        length_bin = prefix_parts[-2]
        trh_bin = prefix_parts[-3]
        label_type = 'NORMAL'

        min_tar_gp = float(parts[1].split('.rsquare')[0])

        df = pd.read_csv(filepath, sep=' ', comment='#', names=['AF_BIN_INDEX', 'AF_BIN_COUNT', 'AF_BIN_MEAN', 'R2_GT', 'R2_DS'])
        for _, row in df.iterrows():
        if row['AF_BIN_COUNT'] > 0:
        r2_vs_af_df_values.append([
        trh_bin, length_bin, min_tar_gp, label_type,
        int(row['AF_BIN_INDEX']), row['AF_BIN_COUNT'], row['AF_BIN_MEAN'], row['R2_DS']
        ])

        r2_vs_af_df = pd.DataFrame(r2_vs_af_df_values, columns=['TRH_BIN', 'LENGTH_BIN', 'MIN_TAR_GP', 'LABEL_TYPE', 'AF_BIN_INDEX', 'AF_BIN_COUNT', 'AF_BIN_MEAN', 'R2_DS'])

        # Statistically aggregate across regions using a weighted average
        if not r2_vs_af_df.empty:
        r2_vs_af_df['WEIGHTED_R2'] = r2_vs_af_df['R2_DS'] * r2_vs_af_df['AF_BIN_COUNT']
        r2_vs_af_df['WEIGHTED_AF'] = r2_vs_af_df['AF_BIN_MEAN'] * r2_vs_af_df['AF_BIN_COUNT']

        r2_vs_af_df = r2_vs_af_df.groupby(['TRH_BIN', 'LENGTH_BIN', 'MIN_TAR_GP', 'LABEL_TYPE', 'AF_BIN_INDEX']).agg({
        'AF_BIN_COUNT': 'sum',
        'WEIGHTED_R2': 'sum',
        'WEIGHTED_AF': 'sum'
        }).reset_index()

        # Recalculate true means from the aggregated weights
        r2_vs_af_df['R2_DS'] = r2_vs_af_df['WEIGHTED_R2'] / r2_vs_af_df['AF_BIN_COUNT']
        r2_vs_af_df['AF_BIN_MEAN'] = r2_vs_af_df['WEIGHTED_AF'] / r2_vs_af_df['AF_BIN_COUNT']

        # 2. Load Error/Concordance Rate Data
        sample_df_values = []
        for filepath in error_files:
        if not filepath.strip(): continue
        filename = os.path.basename(filepath)

        parts = filename.split('_GPfilt_')
        prefix_parts = parts[0].split('.')
        is_info05 = 'INFO05' in prefix_parts

        if is_info05:
        length_bin = prefix_parts[-3]
        trh_bin = prefix_parts[-4]
        label_type = 'INFO05'
        else:
        length_bin = prefix_parts[-2]
        trh_bin = prefix_parts[-3]
        label_type = 'NORMAL'

        min_tar_gp = float(parts[1].split('.error')[0])

        cols = 'GCsV id sample_name #val_gt_RR #val_gt_RA #val_gt_AA #filtered_gp RR_hom_matches RA_het_matches AA_hom_matches RR_hom_mismatches RA_het_mismatches AA_hom_mismatches RR_hom_mismatches_rate_percent RA_het_mismatches_rate_percent AA_hom_mimatches non_reference_discordanc_rate_percent best_gt_rsquared imputed_ds_rsquared'.split(' ')
        df = pd.read_csv(filepath, sep=' ', comment='#', names=cols)

        for _, row in df.iterrows():
        if row['GCsV'] != 'GCsV': continue

        # We need the number of non-reference variants to weight the error rate correctly
        n_nonref = float(row['#val_gt_RA']) + float(row['#val_gt_AA'])
        nrd = float(row['non_reference_discordanc_rate_percent'])

        sample_df_values.append([
        trh_bin, length_bin, min_tar_gp, label_type, row['sample_name'],
        n_nonref, nrd
        ])

        sample_df = pd.DataFrame(sample_df_values, columns=['TRH_BIN', 'LENGTH_BIN', 'MIN_TAR_GP', 'LABEL_TYPE', 'sample_name', 'N_NONREF', 'NRD'])

        # Statistically aggregate across regions using a weighted average
        if not sample_df.empty:
        sample_df['WEIGHTED_NRD'] = sample_df['NRD'] * sample_df['N_NONREF']

        sample_df = sample_df.groupby(['TRH_BIN', 'LENGTH_BIN', 'MIN_TAR_GP', 'LABEL_TYPE', 'sample_name']).agg({
        'N_NONREF': 'sum',
        'WEIGHTED_NRD': 'sum'
        }).reset_index()

        # Recalculate true discordance rate, handling division by zero for sparse samples
        sample_df['non_reference_discordanc_rate_percent'] = np.where(
        sample_df['N_NONREF'] > 0,
        sample_df['WEIGHTED_NRD'] / sample_df['N_NONREF'],
        np.nan
        )
        sample_df = sample_df.dropna(subset=['non_reference_discordanc_rate_percent'])

        # 3. Calculate metrics for title
        num_samples = sample_df['sample_name'].nunique() if not sample_df.empty else 0
        title_metadata = f"Panel: {panel_name}\nTarget: {imputed_name}\nEvaluated samples: {num_samples}"

        # 4. Generate Dosage R2 Plots
        for trh_bin in ['outTRH', 'inTRH']:
        fig, ax = plt.subplots(1, 5, figsize=(10, 2))
        for i, length_bin in enumerate(['SV_DEL', 'DEL', 'SNP', 'INS', 'SV_INS']):
        ax2 = ax[i].twinx()
        for label_type, min_tar_gp in [('NORMAL', 0.0), ('NORMAL', 0.9), ('INFO05', 0.0)]:
        if label_type == 'INFO05':
        label = 'INFO > 0.5'
        ls_val = 'dashed'
        else:
        label = 'unfiltered' if min_tar_gp == 0.0 else f'GP > {min_tar_gp}'
        ls_val = 'solid' if min_tar_gp == 0.0 else 'dotted'

        if not r2_vs_af_df.empty:
        x = (r2_vs_af_df['TRH_BIN'] == trh_bin) & (r2_vs_af_df['LENGTH_BIN'] == length_bin) & (r2_vs_af_df['MIN_TAR_GP'] == min_tar_gp) & (r2_vs_af_df['LABEL_TYPE'] == label_type)

        if not r2_vs_af_df[x].empty:
        ax[i].plot(r2_vs_af_df[x]['AF_BIN_MEAN'], r2_vs_af_df[x]['R2_DS'], label=label,
        ls=ls_val, color='C0')
        ax2.plot(r2_vs_af_df[x]['AF_BIN_MEAN'], r2_vs_af_df[x]['AF_BIN_COUNT'], label=label,
        ls=ls_val, color='C1')

        ax[i].set_xscale('log')
        ax[i].set_xlim([1E-4, 1])
        ax[i].set_ylim([0, 1])

        ax2.set_yscale('log')
        ax2.set_ylim([10**2, 10**9])
        ax2.set_yticks([10**j for j in range(2, 10)])

        if i == 0:
        ax[i].set_ylabel('$r^2_{dosage}$', fontsize=14, color='C0')
        ax2.set_yticklabels([])
        elif i == 4:
        ax2.set_ylabel('number of variants', fontsize=14, color='C1', rotation=270, va='bottom')
        ax[i].set_yticklabels([])
        else:
        ax[i].set_yticklabels([])
        ax2.set_yticklabels([])

        if i == 2:
        trh_tag = {'outTRH': 'non-TR/homopolymer', 'inTRH': 'TR/homopolymer'}[trh_bin]
        ax[i].set_title(f"{title_metadata}\n{trh_tag}\n", fontsize=10)
        ax[i].set_xlabel(f'panel allele frequency\n\n{length_bin}\nALT length - REF length (bp)', fontsize=12)
        ax[i].legend(loc='lower right', fontsize=8)
        else:
        length_bin_label = {'SV_DEL': '(-inf, -50]', 'DEL': '(-50, -1]', 'SNP': 'SNP', 'INS': '[0, 50)', 'SV_INS': '[50, inf)'}[length_bin]
        ax[i].set_xlabel(f'\n\n{length_bin_label}', fontsize=12)

        plt.savefig(f'{output_prefix}.{trh_bin}.r2.png', bbox_inches='tight')
        plt.savefig(f'{output_prefix}.{trh_bin}.r2.pdf', bbox_inches='tight')
        plt.close()

        # 5. Generate Error Plots
        for trh_bin in ['outTRH', 'inTRH']:
        fig, ax = plt.subplots(1, 1, figsize=(6, 3))
        trh_tag = {'outTRH': 'non-TR/homopolymer', 'inTRH': 'TR/homopolymer'}[trh_bin]
        ax.set_title(f"{title_metadata}\n{trh_tag}", fontsize=10)

        plt_df_values = []
        for length_bin in ['SV_DEL', 'DEL', 'SNP', 'INS', 'SV_INS']:
        length_bin_label = {'SV_DEL': '(-inf, -50]', 'DEL': '(-50, -1]', 'SNP': 'SNP', 'INS': '[0, 50)', 'SV_INS': '[50, inf)'}[length_bin]
        for label_type, min_tar_gp in [('NORMAL', 0.0), ('NORMAL', 0.9), ('INFO05', 0.0)]:
        if label_type == 'INFO05':
        min_tar_gp_label = 'INFO > 0.5'
        else:
        min_tar_gp_label = 'unfiltered' if min_tar_gp == 0.0 else f'GP > {min_tar_gp}'

        if not sample_df.empty:
        x = (sample_df['TRH_BIN'] == trh_bin) & (sample_df['LENGTH_BIN'] == length_bin) & (sample_df['MIN_TAR_GP'] == min_tar_gp) & (sample_df['LABEL_TYPE'] == label_type)
        bin_df = sample_df[x]

        for s in range(bin_df.shape[0]):
        plt_df_values.append([length_bin_label, min_tar_gp_label, 1 - 0.01 * bin_df['non_reference_discordanc_rate_percent'].values[s]])

        if plt_df_values:
        plt_df = pd.DataFrame(plt_df_values, columns=['LENGTH_BIN_TEXT', 'MIN_TAR_GP_TEXT', 'non_reference_discordanc_rate_percent'])
        hue_order = ['unfiltered', 'GP > 0.9', 'INFO > 0.5']
        sns.boxplot(data=plt_df, x='LENGTH_BIN_TEXT', y='non_reference_discordanc_rate_percent', hue='MIN_TAR_GP_TEXT', hue_order=hue_order, ax=ax)

        # Apply dashed borders to INFO>0.5 boxes (the last third of drawn patches)
        for j, box in enumerate(ax.patches):
        if j >= 10:
        box.set_linestyle('dashed')

        ax.set_xlabel('ALT length - REF length (bp)', fontsize=8)
        ax.set_ylabel('non-reference concordance rate', fontsize=8)
        ax.set_ylim([0, 1.01])

        handles, labels = ax.get_legend_handles_labels()
        if len(handles) > 2:
        handles[2].set_linestyle('dashed')
        ax.legend(handles=handles, labels=labels, loc='lower center', fontsize=8)

        plt.tight_layout()
        plt.savefig(f'{output_prefix}.{trh_bin}.nrd.png', bbox_inches='tight')
        plt.savefig(f'{output_prefix}.{trh_bin}.nrd.pdf', bbox_inches='tight')
        plt.close()
        EOF

        python3 plot_script.py rsquare_files.list error_files.list "~{panel_name}" "~{imputed_name}" "~{output_prefix}"
    >>>

    output {
        File r2_plot_inTRH_png = "~{output_prefix}.inTRH.r2.png"
        File r2_plot_outTRH_png = "~{output_prefix}.outTRH.r2.png"
        File nrd_plot_inTRH_png = "~{output_prefix}.inTRH.nrd.png"
        File nrd_plot_outTRH_png = "~{output_prefix}.outTRH.nrd.png"

        File r2_plot_inTRH_pdf = "~{output_prefix}.inTRH.r2.pdf"
        File r2_plot_outTRH_pdf = "~{output_prefix}.outTRH.r2.pdf"
        File nrd_plot_inTRH_pdf = "~{output_prefix}.inTRH.nrd.pdf"
        File nrd_plot_outTRH_pdf = "~{output_prefix}.outTRH.nrd.pdf"
    }

    #########################
    RuntimeAttr default_attr = object {
                                   cpu_cores:          2,
                                   mem_gb:             8,
                                   disk_gb:            disk_gb,
                                   boot_disk_gb:       10,
                                   disk_type:          "SSD",
                                   preemptible_tries:  2,
                                   max_retries:        1,
                                   docker:             "jupyter/scipy-notebook:latest"
                               }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " " + select_first([runtime_attr.disk_type, default_attr.disk_type])
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}
