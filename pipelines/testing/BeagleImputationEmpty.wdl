version 1.0

workflow ImputationBeagleEmpty {

    String pipeline_version = "0.0.1"

    input {
        Int chunkLength = 25000000
        Int chunkOverlaps = 5000000 # this is the padding that will be added to the beginning and end of each chunk to reduce edge effects

        File multi_sample_vcf
        File multi_sample_vcf_index

        Boolean perform_extra_qc_steps = false # these are optional additional extra QC steps from Amit's group that should only be
        # run for large sample sets, especially a diverse set of samples (it's further limiting called at sites to 95% and by HWE)
        Float optional_qc_max_missing = 0.05
        Float optional_qc_hwe = 0.000001
        File ref_dict # for reheadering / adding contig lengths in the header of the ouptut VCF, and calculating contig lengths
        Array[String] contigs
        String reference_panel_path # path to the bucket where the reference panel files are stored for all contigs
        String genetic_maps_path # path to the bucket where genetic maps are stored for all contigs
        String output_callset_name # the output callset name
        Boolean split_output_to_single_sample = false

        Int chunks_fail_threshold = 1 # require fewer than this many chunks to fail in order to pass

        # file extensions used to find reference panel files
        String interval_list_suffix = ".interval_list"
        String bref3_suffix = ".bref3"
    }

    output {
        File imputed_multi_sample_vcf = multi_sample_vcf
        File imputed_multi_sample_vcf_index = multi_sample_vcf_index
        File chunks_info = multi_sample_vcf
        File failed_chunks = multi_sample_vcf
        File n_failed_chunks = multi_sample_vcf
    }
}
