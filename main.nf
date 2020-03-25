#!/usr/bin/env nextflow

/*
* MIT License
*
* Copyright (c) 2019 Tobias Neumann, Daniel Malzl
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
  nextflow run dmalzl/inisite-nf

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

if (!params.outputDir) {
  log.info "--outputDir was not specified. Using default directory ./results"
  params.outputDir = "results"
  file(params.outputDir).mkdir()
} else if (!file(params.outputDir).exists()) {
  file(params.outputDir).mkdir()
}

if (!params.filePrefix) {
  params.filePrefix = "NSpeaks"
}

if (!params.extensionSize) {
  log.info "--extensionSize was not specified. Using default of 275"
  params.extensionSize = 275
}

if (!params.qValueCutoff) {
  log.info "--qValueCutoff was not specified. Using default of 0.05"
  params.qValueCutoff = 0.05
}

if (params.treatment2) {

  if (!file(params.treatment2).exists()) {
    exit 1, "--treatment2 was specified but ${params.treatment2} does not exist!"
  }

  if (!params.control2) {
    params.control2 = params.control
  } else if (!file(params.control2).exists()) {
    exit 1, "--control2 was specified but ${params.control2} does not exist!"
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

inputChannel = Channel
                  .from(params.treatment)

if (params.control) {
  process callPeaksWithControl {

    tag { name }

    input:
    from inputChannel

    output:
    into resultCallPeaks

    shell:
    """
    macs2 callpeak -t !{treatment} -c !{control} -f AUTO -g !{genome} -n !{prefix} --outdir !{outdir} --nomodel --extsize !{extensionSize} -q !{qValueCutoff}
    """
  }
} else {
  process callPeaksWithoutControl {

    tag { name }

    input:
    from inputChannel

    output:
    into resultCallPeaks

    shell:
    """
    macs2 callpeak -t !{treatment} -f AUTO -g !{genome} -n !{prefix} --outdir !{outdir} --nomodel --extsize !{extensionSize} -q !{qValueCutoff}
    """
  }
}

process computeInterPeakDistance {

  tag { name }

  input:
  from resultsCallPeaks

  output:
  into resultsInterPeakDistance

  shell:

}

process clusterIntiationSite {

  tag { name }

  input:
  from resultsInterPeakDistance

  output:
  into resultsCluster

  shell:

}

process filterInititiationSites {

  tag { name }

  input:
  from resultsCluster

  output:
  into resultsFilter

  shell:

}
