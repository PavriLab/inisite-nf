#!/usr/bin/env nextflow

/*
* MIT License
*
* Copyright (c) 2020 Daniel Malzl, Tobias Neumann
*
* Permission is hereby granted, free of charge, to any person obtaining a copy
* of this software and associated documentation files (the "Software"), to deal
* in the Software without restriction, including without limitation the rights
* to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the Software is
* furnished to do so, subject to the following conditions:
*
* The above copyright notice and this permission notice shall be included in all
* copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
* OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
* SOFTWARE.
*/
def helpMessage() {
  log.info"""
  ================================================================
   inisite-nf
  ================================================================
  DESCRIPTION

  Reproducible calling of initiation sites from nascent strand sequencing data

  Usage:
  nextflow run dmalzl/inisite-nf --treatment treatment.bam --genome mm [OPTIONS]

  Options:
    --treatment       BAM file containing aligned nacsent strand reads
    --treatment2      BAM file containing aligned nascent strand reads for
                      an additional experiment
    --control         BAM files containing aligned control reads
    --control2        BAM files containing aligned control reads for
                      an additional experiment. If not given --control
                      is reused.


    --extensionSize   size to which the reads should be extended
    --qValueCutoff    MACS2  q-value cutoff
    --genome          genome size to use by MACS (see MACS documentation
                      for allowed values)

    --filePrefix      prefix to use for output files
    --outputDir       name of the directory to save results to
    
    Authors:
    Daniel Malzl (daniel.malzl@imp.ac.at)
    Tobias Neumann (tobias.neumann@imp.ac.at)
  """.stripIndent()
}

params.help = false

if (params.help) {
  helpMessage()
  exit 0
}

if (!params.treatment) {
  exit 1, "--treatment was not specified!"
} else if (!file(params.treatment).exists()) {
  exit 1, "--treatment was specified but ${params.treatment} does not exist!"
}

if (!params.control) {
  log.info "--control was not specified. Proceeding without!"
  params.control = null
} else if (!file(params.control).exists()) {
  exit 1, "--control was specified but ${params.control} does not exist!"
}

if (!params.genome) {
  exit 1, "--genome was not specfied!"
}

if (!file(params.outputDir).exists()) {
  file(params.outputDir).mkdir()
}

if (params.treatment2) {
  if (!file(params.treatment2).exists()) {
    exit 1, "--treatment2 was specified but ${params.treatment2} does not exist!"
  }

  params.control2 = params.control
  if (params.control2) {
    if (!file(params.control2).exists()) {
      exit 1, "--control2 was specified but ${params.control2} does not exist!"
    }
  }

  log.info ""
  log.info " parameters"
  log.info " ======================"
  log.info " treatment                : ${params.treatment}"
  log.info " treatment2               : ${params.treatment2}"
  log.info " control                  : ${params.control}"
  log.info " control2                 : ${params.control2}"
  log.info " extensionSize            : ${params.extensionSize}"
  log.info " qValueCutoff             : ${params.qValueCutoff}"
  log.info " genome                   : ${params.genome}"
  log.info " filePrefix               : ${params.filePrefix}"
  log.info " outputDir                : ${params.outputDir}"
  log.info " ======================"
  log.info ""

} else {
  log.info ""
  log.info " parameters"
  log.info " ======================"
  log.info " treatment                : ${params.treatment}"
  log.info " control                  : ${params.control}"
  log.info " extensionSize            : ${params.extensionSize}"
  log.info " qValueCutoff             : ${params.qValueCutoff}"
  log.info " genome                   : ${params.genome}"
  log.info " filePrefix               : ${params.filePrefix}"
  log.info " outputDir                : ${params.outputDir}"
  log.info " ======================"
  log.info ""

}

if (!params.treatment2 && !params.control) {
  fileList = [[file(params.treatment)]]
} else if (!params.treatment2 && params.control) {
  fileList = [[file(params.treatment), file(params.control)]]
} else if (params.treatment2 && !params.control) {
  fileList = [[file(params.treatment)],
              [file(params.treatment2)]]
} else {
  fileList = [[file(params.treatment), file(params.control)],
              [file(params.treatment2), file(params.control2)]]
}

paramChannel = Channel
                  .fromList([[params.extensionSize,
                              params.qValueCutoff,
                              params.genome,
                              params.outputDir]])

inputChannel = Channel
                  .fromList(fileList)
                  .combine(paramChannel)

if (params.control) {
  process callPeaksWithControl {

    tag { treatment.getSimpleName() }

    publishDir  path: "${params.outputDir}",
                mode: "copy",
                overwrite: "true",
                pattern: "*_MACS.bed"

    input:
    set file(treatment), file(control), val(extensionSize), val(qValueCutoff), val(genome), file(outputDir) from inputChannel

    output:
    file("${treatment.getSimpleName()}_peaks.narrowPeak") into resultsCallPeaks
    file("${treatment.getSimpleName()}_MACS.bed") into calledPeaks

    shell:
    '''
    macs2 callpeak -t !{treatment} -c !{control} -f AUTO -g !{genome} -n !{treatment.getSimpleName()} --nomodel --extsize !{extensionSize} -q !{qValueCutoff}

    grep -v "^#" !{treatment.getSimpleName()}_peaks.xls | grep -v "fold_enrichment" | grep -v "^$" | \\
   	awk \'BEGIN{FS="\\t"; OFS="\\t"} {print $1, $2, $3, $10, $9, "+"}\' > !{treatment.getSimpleName()}_MACS.bed
    '''
  }
} else {
  process callPeaksWithoutControl {

    tag { treatment.getSimpleName() }

    publishDir  path: "${params.outputDir}",
                mode: "copy",
                overwrite: "true",
                pattern: "*_MACS.bed"

    input:
    set file(treatment), val(extensionSize), val(qValueCutoff), val(genome), file(outputDir) from inputChannel

    output:
    file("${treatment.getSimpleName()}_peaks.narrowPeak") into resultsCallPeaks
    file("${treatment.getSimpleName()}_MACS.bed") into calledPeaks

    shell:
    '''
    macs2 callpeak -t !{treatment} -f AUTO -g !{genome} -n !{treatment.getSimpleName()} --nomodel --extsize !{extensionSize} -q !{qValueCutoff}

    grep -v "^#" !{treatment.getSimpleName()}_peaks.xls | grep -v "fold_enrichment" | grep -v "^$" | \\
   	awk \'BEGIN{FS="\\t"; OFS="\\t"} {print $1, $2, $3, $10, $9, "+"}\' > !{treatment.getSimpleName()}_MACS.bed
    '''
  }
}

if (params.treatment2) {
  process intersectTreatments {

    tag { params.filePrefix }

    publishDir  path: "${params.outputDir}",
                mode: "copy",
                overwrite: "true",
                pattern: "*.common.bed"

    input:
    set file(peakFile1), file(peakFile2) from resultsCallPeaks.collect()

    output:
    set val(params.filePrefix), file("${params.filePrefix}.common.bed") into resultsIntersectTreatments

    shell:
    '''
    bedtools intersect -wa -a !{peakFile1} -b !{peakFile2} > AB.intersect.bed
    bedtools intersect -wa -b !{peakFile1} -a !{peakFile2} > BA.intersect.bed
    cat *.intersect.bed | sort -k1,1 -k2,2n > !{params.filePrefix}.intersect.sort.bed
    bedtools merge -i !{params.filePrefix}.intersect.sort.bed -c 4 -o collapse,count,count_distinct > !{params.filePrefix}.common.bed
    '''
  }

  process clusterInitiationSitesFromIntersect {

    tag { filePrefix }

    publishDir  path: "${params.outputDir}",
                mode: "copy",
                saveAs: { filename -> "${filePrefix}_IZ.bed"},
                pattern: "*_clusters.bed"

    input:
    set val(filePrefix), file(commonPeaks) from resultsIntersectTreatments

    output:
    set val(filePrefix), file(commonPeaks), file("${filePrefix}_clusters.bed") into resultsCluster

    shell:
    '''
    clusterinitsites.py -p !{commonPeaks} -o !{filePrefix}
    '''
  }
} else {
  process clusterInitiationSitesFromPeaks {

    tag { params.filePrefix }

    publishDir  path: "${params.outputDir}",
                mode: "copy",
                saveAs: { filename -> "${params.filePrefix}_IZ.bed"},
                pattern: "*_clusters.bed"

    input:
    file(commonPeaks) from resultsCallPeaks

    output:
    set val(params.filePrefix), file(commonPeaks), file("${params.filePrefix}_clusters.bed") into resultsCluster

    shell:
    '''
    clusterinitsites.py -p !{commonPeaks} -o !{params.filePrefix}
    '''
  }
}

process filterInitiationSites {

  tag { filePrefix }

  publishDir  path: "${params.outputDir}",
              mode: "copy",
              overwrite: "true",
              pattern: "*_IS.bed"

  input:
  set val(filePrefix), file(commonPeaks), file(iniZones) from resultsCluster

  output:
  file("${filePrefix}_IS.bed") into resultsFilter

  shell:
  '''
  bedtools intersect -u -wa -a !{commonPeaks} -b !{iniZones} > !{filePrefix}_IS.bed
  '''
}

workflow.onComplete {
	println ( workflow.success ? "COMPLETED!" : "FAILED" )
}
