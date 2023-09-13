process BEDTOOLS_FILTER_IS {
    tag "meta.id"

    input:
    tuple val(meta), path(commonPeaks)
    tuple val(meta), path(clusters)

    output:
    tuple val(meta), path("${meta.id}_IS.bed"), emit: filteredbed

    shell:
    """
    bedtools intersect \
        -u -wa \
        -a !{commonPeaks} \
        -b !{clusters} > \
        !{meta.id}_IS.bed
    """
}
