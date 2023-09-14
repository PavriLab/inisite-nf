process DEEPTOOLS_GENERATE_BIGWIG {
    tag "meta.id"

    input:
    tuple val(meta), path(bam), path(index)
    val(binsize)

    output:
    tuple val(meta), path("${meta.id}.bw"),     emit: bigwig

    shell:
    """
    bamCoverage \
        -b !{bam} \
        -o bigwigs/!{meta.id}.bw \
        -of bigwig \
        -bs !{binsize} \
        -p 8 \
        --ignoreDuplicates \
        --normalizeUsing CPM \
        --exactScaling
    """
}
