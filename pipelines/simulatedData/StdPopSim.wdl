version 1.0

workflow StdPopSim {

    input {
        Int numEuropeanSamples
        Int numAfricanSamples
        String basename
        Int disk_size_gb_make_tree
        Int disk_size_gb_make_vcf
        Array[String] contigs
    }

    scatter(contig in contigs){
        call makeTreeFromStdPopSim {
            input:
                numEuropeanSamples = numEuropeanSamples,
                numAfricanSamples = numAfricanSamples,
                disk_size_gb = disk_size_gb_make_tree,
                contig = contig

        }

        call makeVcfFromStdPopSimTree {
            input:
                treeFile = makeTreeFromStdPopSim.output_tree,
                basename = basename,
                disk_size_gb = disk_size_gb_make_vcf,
                contig = contig

        }
    }

    output{
        Array[File] vcf = makeVcfFromStdPopSimTree.output_vcf
        Array[File] vcf_indices = makeVcfFromStdPopSimTree.output_vcf_index
    }
}

task makeTreeFromStdPopSim {
    input {
        Int numEuropeanSamples
        Int numAfricanSamples
        String contig
        Int disk_size_gb
    }

    command {
        set -eo pipefail
        stdpopsim HomSap -s 1046 -g HapMapII_GRCh38 -c ~{contig} -o sim_tree.ts -d OutOfAfrica_2T12 EUR:~{numEuropeanSamples} AFR:~{numAfricanSamples}
    }

    runtime {
        docker: "jsotobroad/stdpopsim:latest"
        disks: "local-disk ${disk_size_gb} SSD"
        memory: "64 GiB"
    }

    output {
        File output_tree = "sim_tree.ts"
    }
}

task makeVcfFromStdPopSimTree {
    input {
        File treeFile
        String basename
        String contig
        Int disk_size_gb
    }

    command {
        set -eo pipefail

        tskit vcf --contig-id ~{contig} ~{treeFile} | bgzip -c > ~{basename}.~{contig}.vcf.gz
        tabix ~{basename}.~{contig}.vcf.gz
    }

    runtime {
        docker: "jsotobroad/stdpopsim:latest"
        disks: "local-disk ${disk_size_gb} SSD"
        memory: "64 GiB"
    }

    output {
        File output_vcf = "~{basename}.~{contig}.vcf.gz"
        File output_vcf_index = "~{basename}.~{contig}.vcf.gz.tbi"
    }
}
