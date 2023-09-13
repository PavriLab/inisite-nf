process BOWTIE_BUILD_INDEX {

    tag "${bwt2_base}"

    input:
    file(genomeFasta)

    output:
    tuple val(bwt2_base), path("bowtie2Index"), emit: index

    script:
    bwt2_base = genomeFasta.getSimpleName()

    """
    mkdir bowtieIndex

    bowtie-build \
        ${genomeFasta} \
        bowtieIndex/${bwt_base} \
        --threads !{task.cpus}
    """
}
