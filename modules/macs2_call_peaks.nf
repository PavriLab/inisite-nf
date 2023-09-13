process MACS2_CALL_PEAKS {
    tag "meta.id"

    input:
    tuple val(meta), path(treatment)
    path(control)
    val(extensionSize)
    val(qValueCutoff)

    output:
    tuple val(meta), path("${meta.id}_MACS.bed"),           emit: bed
    tuple val(meta), path("${meta.id}_peaks.narrowPeak"),   emit: narrowpeaks

    script:
    ext_args = control ? "-c ${control}" : ""

    """
    macs2 callpeak \
        !{ext_args} \
        -t !{treatment} \
        -f AUTO \
        -g !{meta.genomesize} \
        -n !{meta.id} \
        --nomodel \
        --extsize !{extensionSize} \
        -q !{qValueCutoff}

    grep -v "^#" !{meta.id}_peaks.xls | \
    grep -v "fold_enrichment" | \
    grep -v "^$" | \
    awk \'BEGIN{FS="\\t"; OFS="\\t"} {print $1, $2, $3, $10, $9, "+"}\' > \
    !{meta.id}_MACS.bed
    """
}
