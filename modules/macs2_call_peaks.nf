process MACS2_CALL_PEAKS {
    tag "meta.id"

    input:
    tuple val(meta), path(treatment), path(control)
    val(extensionSize)
    val(qValueCutoff)
    val(macsGenomeSize)

    output:
    tuple val(meta), path("${meta.id}_MACS.bed"),           emit: bed
    tuple val(meta), path("${meta.id}_peaks.narrowPeak"),   emit: narrowpeaks

    script:
    def control = control ? "-c ${control}" : ""

    """
    macs2 callpeak \
        !{control} \
        -t !{treatment} \
        -f AUTO \
        -g !{macsGenomeSize} \
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
