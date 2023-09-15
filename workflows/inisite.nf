/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
WorkflowInisite.initialise(params, log)

if (params.genome && params.genomes && !params.igenomes_ignore) {
    igenomes_bowtie     = WorkflowInisite.getGenomeAttribute(params, 'bowtie')
    igenomes_fasta      = WorkflowInisite.getGenomeAttribute(params, 'fasta')
    igenomes_chromSizes = WorkflowInisite.getGenomeAttribute(params, 'chromSizes')

} else {
    igenomes_bowtie = ''
    igenomes_fasta = ''
    igenomes_chromSizes = ''
}

// Check input path parameters to see if they exist
checkPathParamList = [
    params.input,
    params.fasta,
    params.chromSizes,
    igenomes_fasta,
    igenomes_chromSizes
]

for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// setting up for prepare genome subworkflow
def prepare_genome_for_tools = []

// if we do not have --genome
if (!params.genome) {
    // bowtie index
    if (params.fasta) {
        prepare_genome_for_tools << "bowtie"

    } else {
        log.error "Neither --genome nor --fasta are specified but needed for bowtie2 index."
        System.exit(1)
    }

    // HICUP genome digest
    if (!(params.fasta && params.chromSizes)) {
        log.error "--fasta or --chromSizes not given and --genome not specified"
        System.exit(1)
    }

    if (params.chromSizes.endsWith("xml")) {
        prepare_genome_for_tools << "chromSizes"
    }
// if --genome is specified we check if everything is there
} else {
    if (!igenomes_bowtie) {
        log.info "Bowtie index not found in igenomes config file. Computing from igenomes_fasta"
        prepare_genome_for_tools << "bowtie"
    }

    if (igenomes_chromSizes.endsWith("xml")) {
        prepare_genome_for_tools << "chromSizes"
    }
}

dynamic_params = [:]
dynamic_params.genomeFasta      = params.genome ? igenomes_fasta : params.fasta
dynamic_params.genomeSizes      = params.genome ? igenomes_chromSizes : params.chromSizes
dynamic_params.bowtieIndex      = igenomes_bowtie ? igenomes_bowtie : "computed from fasta"
// dynamic_params.genomeSizeType   = WorkflowInisite.getGenomeSizeType( dynamic_params.genomeSizes )
dynamic_params.genomeName       = params.genome ? params.genome : file(dynamic_params.genomeFasta).getSimpleName()

WorkflowInisite.paramsSummaryLog( params, dynamic_params, log )

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { INPUT_CHECK                   } from '../subworkflows/input_check.nf'
include { PREPARE_GENOME                } from '../subworkflows/prepare_genome.nf'
include { TRIM_GALORE                   } from '../modules/trim_galore.nf'
include { BOWTIE_ALIGN                  } from '../modules/bowtie_align.nf'
include { DEEPTOOLS_GENERATE_BIGWIG     } from '../modules/deeptools_generate_bigwigs.nf'
include { MACS2_CALL_PEAKS              } from '../modules/macs2_call_peaks.nf'
include { BEDTOOLS_INTERSECT_REPLICATES } from '../modules/bedtools_intersect_replicates.nf'
include { CLUSTERSCAN_CLUSTER_IS        } from '../modules/clusterscan_cluster_initiation_sites.nf'
include { BEDTOOLS_FILTER_IS            } from '../modules/bedtools_filter_initiation_sites.nf'
include { OVERLAP_SAMPLES               } from '../modules/overlap_samples.nf'
include { MULTIQC                       } from '../modules/multiqc.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow INISITE {
    ch_input = file( params.input )

    // check input sample sheet
    // adapted from nf-core/rnaseq
    INPUT_CHECK ( ch_input )
        .reads
        .groupTuple(by: [0])
        .branch {
            meta, fastq ->
                single  : fastq.size() == 1
                    return [ meta, fastq.flatten() ]
                multiple: fastq.size() > 1
                    return [ meta, fastq.flatten() ]
        }
        .set { ch_fastq }

    // prepare genome files
    if (!prepare_genome_for_tools.isEmpty()) {
        ch_genome = PREPARE_GENOME (
            prepare_genome_for_tools,
            dynamic_params
        )

    } else {
        ch_genome = [:]
        ch_genome.index     = file( dynamic_params.bowtie2Index ).getParent()
        ch_genome.sizes     = file( dynamic_params.genomeSizes )
    }

    // read QC
    TRIM_GALORE ( ch_cat_fastq )
        .reads
        .set { ch_trim_fastq }

    BOWTIE_ALIGN (
        ch_trim_fastq,
        ch_genome.bowtieIndex
    )

    DEEPTOOLS_GENERATE_BIGWIG (
        BOWTIE_ALIGN.alignments,
        params.binSize
    )

    // adapted from nf-core chipseq
    BOWTIE_ALIGN.out.alignments
        .combine ( BOWTIE_ALIGN.alignments )
        .map {
            meta1, bam1, bai1, meta2, bam2, bai2 ->
                meta1.control == meta2.id ? [ meta1, [ bam1, bam2 ]] : null
        }
        .set { ch_callpeaks_input }

    MACS2_CALL_PEAKS (
        ch_callpeaks_input,
        params.extensionSize,
        params.qValueCutoff
        macsGenomeSize
    )

    MACS2_CALL_PEAKS.out.bed
        .map {
            meta, bed ->
                [ meta.replicate, meta, bed ]
        }
        .groupTuple(by: [0])
        .branch {
            meta, beds ->
                single      :   beds.size() == 1
                    return [ meta, beds.flatten() ]

                replicates  :   beds.size() > 1
                    return [ meta, beds.flatten() ]
        }
        .set { ch_called_peaks }

    BEDTOOLS_INTERSECT_REPLICATES ( ch_called_peaks.replicates )

    BEDTOOLS_INTERSECT_REPLICATES.commonPeaks
        .mix ( ch_called_peaks.single )
        .set { ch_cluster_is }

    CLUSTERSCAN_CLUSTER_IS ( ch_cluster_is )

    BEDTOOLS_FILTER_IS ()

    OVERLAP_SAMPLES ()

    MULTIQC (
        TRIM_GALORE.out.fastqc.collect(),
        TRIM_GALORE.out.reports.collect(),
        ch_hicup.multiqc.collect()
    )
}

workflow.onComplete {
	println ( workflow.success ? "COMPLETED!" : "FAILED" )
}
