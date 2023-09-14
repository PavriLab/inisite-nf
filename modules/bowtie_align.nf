process BOWTIE_ALIGN {
    tag "meta.id"

    input:
    tuple val(meta), path(reads)
    path(index)

    output:
    tuple val(meta), path("${meta.id}.bam"), path("${meta.id}.bam.bai"),    emit: alignments
    tuple val(meta), path("*.{flagstat,idxstats,stats}"),                   emit: stats

    shell:
    '''
    bowtie \
        -S \
        -p !{task.cpus} \
        --trim5 0 \
        --trim3 0 \
        -v 2 \
        --best \
        --strata \
        --tryhard \
        -m 1 \
        --phred33-quals \
        --chunkmbs 256 \
        !{index}/!{bwt_base} \
        !{reads} | \
    samtools sort \
        -@ !{task.cpus} \
        -O BAM \
        -o !{meta.id}.bam \

    samtools index !{meta.id}.bam

    samtools flagstat !{meta.id}.bam > !{meta.id}.bam.flagstat
    samtools idxstats !{meta.id}.bam > !{meta.id}.bam.idxstats
    samtools stats !{meta.id}.bam > !{meta.id}.bam.stats
    '''
}
