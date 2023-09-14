include { GUNZIP as GUNZIP_FASTA } from '../modules/gunzip.nf'
include { XML_TO_TSV             } from '../modules/xml_to_tsv.nf'
include { DIGEST_GENOME          } from '../modules/hicup_digest_genome.nf'
include { BOWTIE2_BUILD_INDEX    } from '../modules/bowtie2_build_index.nf'

workflow PREPARE_GENOME {
    take:
    prepare_genome_for_tools
    dynamic_params

    main:

    // Uncompress genome fasta file if required
    if (params.fasta.endsWith('.gz')) {
        ch_fasta = GUNZIP_FASTA (
             file( dynamic_params.genomeFasta ) // needs to be wrapped in file for GUNZIP to recognize as input
        )

    } else {
        ch_fasta = file( dynamic_params.genomeFasta )
    }

    if ("bowtie" in prepare_genome_for_tools) {
        ch_bowtie2_index = BOWTIE_BUILD_INDEX (
            ch_fasta,
            dynamic_params.genomeSizeType
        )

    } else {
        def bwt_base = file( dynamic_params.bowtieIndex ).getSimpleName()
        def bwt_dir = file( dynamic_params.bowtieIndex ).getParent()
        ch_bowtie_index = [ bwt_base, bwt_dir ]
    }

    if ("chromSizes" in prepare_genome_for_tools) {
        ch_genome_sizes = XML_TO_TSV (
            file( dynamic_params.genomeSizes )
        )

    } else {
        ch_genome_sizes = file( dynamic_params.genomeSizes )
    }

    emit:
    index      = ch_bowtie_index
    sizes      = ch_genome_sizes

}
