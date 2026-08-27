# WDLs and notebooks used for scientific validation of SV imputation

## GLIMPSE2Concordance

### Purpose
This WDL evaluates concordance between imputed variant call sets and a reference
panel using GLIMPSE2's concordance tool. Variants are stratified by tandem
repeat / homopolymer (TRH) status and by variant length bin (SV_DEL, DEL, SNP,
INS, SV_INS). Concordance is computed both on unfiltered genotypes, genotypes
filtered by posterior probability (GP > 0.9), and variants filtered by imputation
INFO score (INFO > 0.5). Per-chromosome and aggregate plots are produced for
dosage r² and non-reference discordance rate (NRD).

#### Inputs
* `panel_vcfs` / `panel_vcf_idxs` — Reference panel VCFs (one per chromosome/region)
* `imputed_vcfs` / `imputed_vcf_idxs` — Imputed VCFs (one per chromosome/region)
* `trh_bed` / `trh_bed_idx` — BED file defining tandem repeat / homopolymer regions
* `regions` — List of chromosomal regions matching the VCF arrays
* `output_prefix` — Prefix for all output file names
* `remap_file` *(optional)* — Two-column file used to rename samples before concordance

#### Outputs
* `concordance_aggregate_plots_pdf` / `concordance_aggregate_plots_png`
  * Aggregate r² and NRD plots for inTRH and outTRH variants (4 files each)
* `concordance_per_chrom_plots_pdf` / `concordance_per_chrom_plots_png`
  * Per-chromosome r² and NRD plots
* `concordance_results`
  * Raw GLIMPSE2 concordance result files:
    `*.rsquare.grp.txt.gz`, `*.rsquare.spl.txt.gz`, `*.error.grp.txt.gz`,
    `*.error.spl.txt.gz`, `*.error.cal.txt.gz`

---

## VcfdistEvaluation

### Purpose
This WDL evaluates imputation accuracy against truth call sets using
[vcfdist](https://github.com/TimD1/vcfdist), computing precision, recall, and
F1 score for SNPs, INDELs, and SVs. Samples are evaluated independently within
confident regions, and results can be stratified across multiple BED-defined
variant subsets. A summary TSV aggregating metrics across all samples and
stratifications is produced.

#### Inputs
* `eval_sample_names` — Sample names as they appear in the evaluation VCF
* `truth_sample_names` — Sample names as they appear in the truth VCFs (e.g. `"syndip"` for dipcall truth sets)
* `eval_vcf` / `eval_vcf_idx` — Evaluation VCF containing all samples
* `truth_vcfs` / `truth_vcf_idxs` *(optional)* — Per-sample truth VCFs
* `confident_regions_bed_files` — Per-sample BED files defining confident regions
* `region` — Chromosomal region to evaluate
* `reference_fasta` / `reference_fasta_fai` — Reference genome FASTA and index
* `do_naively_phase` — Convert unphased (`/`) genotypes to phased (`|`) before evaluation
* `vcfdist_bed_files` — BED files defining variant stratifications
* `labels_per_stratification` — Human-readable labels for each stratification BED
* `vcfdist_extra_args` *(optional)* — Additional arguments passed to vcfdist
* `summarize_evaluations_docker` — Docker image used for the summarization step

#### Outputs
* `vcfdist_summary`
  * Per-stratification, per-sample vcfdist output files:
    summary VCF, precision-recall TSVs, query/truth TSVs, phasing TSVs
* `evaluation_summary_tsv`
  * Aggregated TSV with mean precision, recall, and F1 score for SNPs, INDELs,
    and SVs across all evaluated samples, broken down by stratification label

---

## MendelianConsistency

### Purpose
This WDL evaluates Mendelian consistency of imputed variant call sets across
family trios. Variants are annotated with tandem repeat / homopolymer (TRH)
status and allele frequency from a panel sites-only VCF, then Mendelian error
rates are computed per trio and per locus. Metrics are stratified by TRH status,
variant length bin (SV_DEL, DEL, SNP, INS, SV_INS), and allele frequency bin.
Error rates are computed on unfiltered genotypes, genotypes filtered by posterior
probability (GP > 0.9), and variants filtered by imputation INFO score
(INFO > 0.5). Per-chromosome and aggregate boxplots are produced for both
per-trio and per-locus error rates.

#### Inputs
* `panel_sites_only_vcfs` / `panel_sites_only_vcf_idxs` — Reference panel sites-only VCFs split to biallelic (one per chromosome/region)
* `imputed_vcfs` / `imputed_vcf_idxs` — Imputed VCFs split to biallelic (one per chromosome/region)
* `trh_bed` / `trh_bed_idx` — BED file defining tandem repeat / homopolymer regions
* `pedigree` — Tab-delimited pedigree file (FAM format) defining family trios
* `output_prefix` — Prefix for all output file names

#### Outputs
* `mendelian_aggregate_plots_pdf` / `mendelian_aggregate_plots_png`
  * Aggregate per-trio and per-locus Mendelian error rate boxplots for inTRH
    and outTRH variants (4 files each)
* `mendelian_per_chrom_plots_pdf` / `mendelian_per_chrom_plots_png`
  * Per-chromosome per-trio and per-locus Mendelian error rate boxplots

---

## plot_concordance_differences.ipynb

### Purpose
This Jupyter notebook visually compares GLIMPSE2 concordance results between
two workflow runs (Run A and Run B). It downloads rsquare and error output files
from GCS (or reads them locally), aggregates metrics using weighted averages
across chromosomes, and produces absolute-difference plots showing how dosage r²
and NRD change between the two runs. Plots are stratified by TRH bin and variant
length bin and are saved as PNG and PDF.

#### Inputs *(configured at the top of the notebook)*
* `PATH_RSQUARE_A` / `PATH_RSQUARE_B` — GCS prefix or local path containing rsquare files for each run
* `PATH_ERROR_A` / `PATH_ERROR_B` — GCS prefix or local path containing error files for each run
* `LABEL_RUN_A` / `LABEL_RUN_B` — Display labels for the two runs
* `OUTPUT_PREFIX` — Prefix for saved plot files

#### Outputs
* `<output_prefix>.<trh_bin>.diff.png` / `.pdf`
  * Absolute difference plots for each TRH bin (inTRH, outTRH)

---

## plot_vcfdist_evaluation_comparison.ipynb

### Purpose
This Jupyter notebook visually compares vcfdist evaluation summaries between
two workflow runs (Run A and Run B). It reads the `evaluation_summary.tsv`
produced by `VcfdistEvaluation` from local paths or GCS, and generates
bar-chart comparison plots for precision, recall, and F1 score across variant
types (SNP, INDEL, SV) and stratification labels. An optional fixed y-axis
range can be set so that plots from separate notebook runs remain directly
comparable.

#### Inputs *(configured at the top of the notebook)*
* `PATH_TSV_A` / `PATH_TSV_B` — Local or GCS path to each run's `evaluation_summary.tsv`
* `LABEL_RUN_A` / `LABEL_RUN_B` — Display labels for the two runs
* `OUTPUT_PREFIX` — Prefix for saved plot files
* `FIXED_DIFF_Y_RANGE` — Fixed y-axis range for difference plots (set to `None` for autoscale)

#### Outputs
* `<output_prefix>_comparison.<variant_type>.png` / `.pdf`
  * Comparison and difference plots per variant type (SNP, INDEL, SV)

