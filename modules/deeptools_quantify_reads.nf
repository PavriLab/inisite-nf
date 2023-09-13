process DEEPTOOLS_QUANTIFY_READS {
    tag "meta.id"

    input:

    output:

    shell:
    """
    """
}

process quantifyReads {

  tag { params.filePrefix }

  publishDir  path: "${params.outputDir}",
              mode: "copy",
              overwrite: "true",
              pattern: "*.count.tsv"

  input:
  set file(bamA), val(labelA), file(bamB), val(labelB) from quantificationChannel
  set file(bamAindex), file(bamBindex) from indexChannel
  file(mergedSites) from resultsMergeSites

  output:
  set val(labelA), val(labelB), file(mergedSites), file("${params.filePrefix}.count.tsv") into resultsQuantifyReads

  shell:
  '''
  multiBamSummary BED-file -b !{bamA} !{bamB} \
                           --BED !{mergedSites} \
                           -l !{labelA} !{labelB} \
                           -o !{params.filePrefix}.counts.npz \
                           --outRawCounts !{params.filePrefix}.count.tsv \
                           -p 16 \
                           --ignoreDuplicates
  '''
}
