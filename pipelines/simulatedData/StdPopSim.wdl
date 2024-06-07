version 1.0

workflow StdPopSim {

    input {
        Int numEuropeanSamples
        Int numAfricanSamples
        String basename
        Int disk_size_gb
        Array[String] contigs
    }

    scatter(contig in contigs){
        call makeVcfFromStdPopSim {
            input:
                numEuropeanSamples = numEuropeanSamples,
                numAfricanSamples = numAfricanSamples,
                basename = basename,
                disk_size_gb = disk_size_gb,
                contig = contig

        }
    }

    output{
        Array[File] vcf = makeVcfFromStdPopSim.vcf_output
        Array[File] vcf_indices = makeVcfFromStdPopSim.vcf_output_index
    }
}

task makeVcfFromStdPopSim {
    input {
        Int numEuropeanSamples
        Int numAfricanSamples
        String basename
        String contig
        Int disk_size_gb
    }

    command {
        set -eo pipefail
        
        stdpopsim HomSap -s 1046 -g HapMapII_GRCh38 -c ~{contig} -o sim_tree.ts -d OutOfAfrica_2T12 EUR:~{numEuropeanSamples} AFR:~{numAfricanSamples}
        tskit vcf --contig-id ~{contig} sim_tree.ts | bgzip -c > ~{basename}.vcf.gz
        tabix ~{basename}.vcf.gz
    }

    runtime {
        docker: " jsotoimputation.azurecr.io/stdpopsim:latest"
        disks: "local-disk ${disk_size_gb} SSD"
        memory: "16 GiB"
    }

    output {
        File output_vcf = "~{basename}.vcf.gz"
        File output_vcf_index = "~{basename}.vcf.gz.tbi"
    }
}
