process BEDTOOLS_INTERSECT_REPLICATES {
    tag "meta.id"

    input:
    tuple val(meta), path(replicatePeaks)

    output:
    tuple val(meta), path("${meta.id}.common.bed"),     emit: commonPeaks

    shell:
    """
    """
}
