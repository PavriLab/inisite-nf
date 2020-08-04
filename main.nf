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
  nextflow run pavirlab/inisite-nf --treatment treatment.bam --genome mm9 [OPTIONS]

  Options:
    --treatment       fastq file containing aligned nacsent strand reads
    --treatment2      fastq file containing aligned nascent strand reads for
                      an additional experiment (Default: null)
    --control         fastq file containing aligned control reads (Default: null)
    --control2        fastq file containing aligned control reads for
                      an additional experiment. If not given --control
                      is reused. (Default: null)


    --extensionSize   size to which the reads should be extended (Default: 275bp)
    --qValueCutoff    MACS2  q-value cutoff (Default: 0.05)

    --filePrefix      prefix to use for output files (Default: NSpeaks)
    --outputDir       name of the directory to save results to (Default: results)

    References:
     --genome         Name of reference (hg38, mm10, ...)
     --fasta          Alternatively, path to genome fasta file which will be digested
     --genomeSize     effective genome size to use with MACS
     --bowtieIndex    Optional: Path to bowtie index

  Profiles:
     standard         local execution
     singularity      local execution with singularity
     cbe              CBE cluster execution with singularity

  Docker:
  pavrilab/inisite-nf:latest

  Authors:
  Tobias Neumann (tobias.neumann@imp.ac.at)
  Daniel Malzl (daniel.malzl@imp.ac.at)
  """.stripIndent()
}

params.help = false
igenomes_fasta = params.genomes[ params.genome ].fasta ?: false
igenomes_bowtie = params.genomes[ params.genome ].bowtie ?: false
igenomes_genomeSize = params.genomes[ params.genome ].genomeSize ?: false


if (params.help) {
  helpMessage()
  exit 0
}

if (params.genomeSize) {
    effectiveGenomeSize = params.genomeSize

} else if (igenomes_genomeSize) {
    effectiveGenomeSize = igenomes_genomeSize

} else {
  exit 1, "no effectve genome size given to use with MACS, exiting"

}

if (!params.bowtieIndex) {
  if (!params.fasta && !igenomes_bowtie) {
    exit 1, "Fasta needed for BowtieIndex not specified!"

  } else if (params.fasta && !igenomes_fasta) {
    Channel
        .fromPath(params.fasta, checkIfExists: true)
        .ifEmpty { exit 1, "Fasta needed for BowtieIndex but not found at ${params.fasta}" }
    fastaFile = params.fasta

  } else {
    Channel
        .fromPath(igenomes_fasta, checkIfExists: true)
        .ifEmpty { exit 1, "Fasta needed for BowtieIndex but not given and not found in igenomes or igenomes not present."}
    fastaFile = igenomes_fasta

  }

} else {
  log.info "bowtieIndex is specified explicitly. Fasta file will not be used if given!"
  fastaFile = "not used due to set --bowtieIndex"
}

if (!params.bowtieIndex) {
  if (params.fasta) {
    lastPath = fastaFile.lastIndexOf(File.separator)
    bwt_base = fastaFile.substring(lastPath+1) - ~/(\.fa)?(\.fasta)?(\.fas)?$/

    fastaForBowtie = Channel
                        .fromPath(fastaFile)
    bowtieIndexFile = 'computed from fasta'
    makeBowtieIndex = true

  }  else if (igenomes_bowtie) {
    lastPath = igenomes_bowtie.lastIndexOf(File.separator)
    bwt_dir  = igenomes_bowtie.substring(0,lastPath+1)
    bwt_base = igenomes_bowtie.substring(lastPath+1)

    bowtieIndex = Channel
                      .fromPath(bwt_dir , checkIfExists: true)
                      .ifEmpty { exit 1, "Genome index: Provided index not found: ${igenomes_bowtie}" }
    bowtieIndexFile = igenomes_bowtie
    makeBowtieIndex = false

  } else {
    exit 1, "BowtieIndex was not specified! Use either --fasta or make sure an igenomes database is properly configured."
  }

} else {
  bowtieIndexFile = params.bowtieIndex
  lastPath = params.bowtieIndex.lastIndexOf(File.separator)
  bwt_dir  = params.bowtieIndex.substring(0,lastPath+1)
  bwt_base = params.bowtieIndex.substring(lastPath+1)
  bowtieIndex = Channel
                    .fromPath(bwt_dir , checkIfExists: true)
                    .ifEmpty { exit 1, "Genome index: Provided index not found: ${params.bowtieIndex}" }
  makeBowtieIndex = false

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

  if (params.control && !params.control2) {
    exit 1, "both or none controls have to be given. You can use the same control for both arguments but be sure they have distinctive names!"
  
  } else if (!file(params.control2).exists()) {
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
  log.info " Fasta                    : ${fastaFile}"
  log.info " bowtieIndex              : ${bowtieIndexFile}"
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
  log.info " Fasta                    : ${fastaFile}"
  log.info " bowtieIndex              : ${bowtieIndexFile}"
  log.info " filePrefix               : ${params.filePrefix}"
  log.info " outputDir                : ${params.outputDir}"
  log.info " ======================"
  log.info ""

}

if (makeBowtieIndex) {
  process buildBowtieIndex {

    tag "${bwt_base}"

    input:
    file(fasta) from fastaForBowtie

    output:
    file("bowtieIndex") into bowtieIndex

    shell:
    """
    mkdir bowtieIndex

    bowtie-build ${fasta} bowtieIndex/${bwt_base} --threads !{task.cpus}
    """

  }
}

if (!params.treatment2 && !params.control) {
  files2Align = [[1, file(params.treatment)]]

} else if (!params.treatment2 && params.control) {
  files2Align = [[1, file(params.treatment)], [1, file(params.control)]]

} else if (params.treatment2 && !params.control) {
  files2Align = [[1, file(params.treatment)], [2, file(params.treatment2)]]

} else {
  files2Align = [[1, file(params.treatment)], [1, file(params.control)],
                 [2, file(params.treatment2)], [2, file(params.control2)]]
}

preprocessChannel = Channel
                        .fromList(files2Align)

process trimReads {
    tag { name }

    input:
    set val(num), file(fastqFile) from preprocessChannel

    output:
    set val(num), val(name), file("${name}_trimmed.fq") into alignChannel
    file("*_fastqc.{zip,html}") into fastqcResults
    file("*trimming_report.txt") into trimgaloreResults

    shell:
    lastPath = fastqFile.toString().lastIndexOf(File.separator)
    readBase = fastqFile.toString().substring(lastPath+1)
    name = fastqFile.toString() - ~/(\.fq)?(\.fq\.gz)?(\.fastq)?(\.fastq\.gz)?$/
    '''
    trim_galore --quality 20 \
                --fastqc \
                --length 18 \
                --illumina \
                --dont_gzip \
                --basename !{name} \
                --cores !{task.cpus} \
                !{fastqFile} \

    mv !{readBase}_trimming_report.txt !{name}_trimmed.fq_trimming_report.txt
    sed -i 's/Command line parameters:.*\$/Command line parameters: !{name}_trimmed/g' !{name}_trimmed.fq_trimming_report.txt
    '''

}

process alignReads {
    tag { name }

    publishDir  path: "${params.outputDir}/bams",
                mode: "copy",
                overwrite: "true",
                pattern: "*.bam"

    input:
    set val(num), val(name), file(trimmed) from alignChannel
    file(index) from bowtieIndex.collect()

    output:
    set val(num), val(name), file("${name}.bam") into alignOutputChannel
    file "*.{flagstat,idxstats,stats}" into bowtieMultiqcChannel

    shell:
    '''
    bowtie -S \
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
           !{trimmed} \
           !{name}.sam

    samtools sort -@ !{task.cpus} \
                  -O BAM \
                  -o !{name}.bam \
                  !{name}.sam

    samtools index !{name}.bam
    samtools flagstat !{name}.bam > !{name}.bam.flagstat
    samtools idxstats !{name}.bam > !{name}.bam.idxstats
    samtools stats !{name}.bam > !{name}.bam.stats
    '''

}

macsParamChannel = Channel
                      .fromList([[params.extensionSize,
                                  params.qValueCutoff,
                                  effectiveGenomeSize,
                                  params.outputDir]])

macsInputChannel = alignOutputChannel
                      .groupTuple()
                      .combine(macsParamChannel)
		      .println()

if (params.control) {
  process callPeaksWithControl {

    tag { treatname }

    publishDir  path: "${params.outputDir}/peaks",
                mode: "copy",
                overwrite: "true",
                pattern: "*_MACS.bed"

    input:
    set val(num), val(treatname), file(treatment), val(controlname), file(control), val(extensionSize), val(qValueCutoff), val(genomeSize), file(outputDir) from macsInputChannel

    output:
    file("${treatment.getSimpleName()}_peaks.narrowPeak") into resultsCallPeaks
    file("${treatment.getSimpleName()}_MACS.bed") into calledPeaks

    shell:
    '''
    macs2 callpeak -t !{treatment} -c !{control} -f AUTO -g !{genomeSize} -n !{treatment.getSimpleName()} --nomodel --extsize !{extensionSize} -q !{qValueCutoff}

    grep -v "^#" !{treatment.getSimpleName()}_peaks.xls | grep -v "fold_enrichment" | grep -v "^$" | \\
   	awk \'BEGIN{FS="\\t"; OFS="\\t"} {print $1, $2, $3, $10, $9, "+"}\' > !{treatment.getSimpleName()}_MACS.bed
    '''

  }

} else {
  process callPeaksWithoutControl {

    tag { treatment.getSimpleName() }

    publishDir  path: "${params.outputDir}/peaks",
                mode: "copy",
                overwrite: "true",
                pattern: "*_MACS.bed"

    input:
    set val(num), val(name), file(treatment), val(extensionSize), val(qValueCutoff), val(genome), file(outputDir) from macsInputChannel

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

    publishDir  path: "${params.outputDir}/peaks",
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

    publishDir  path: "${params.outputDir}/peaks",
                mode: "copy",
                saveAs: { filename -> "${filePrefix}_IZ.bed"},
                pattern: "*_clusters.bed"

    input:
    set val(filePrefix), file(commonPeaks) from resultsIntersectTreatments

    output:
    set val(filePrefix), file(commonPeaks), file("${filePrefix}_clusters.bed") into resultsCluster

    shell:
    clusterScanScript = "${NXF_HOME}/assets/pavrilab/hicer-nf/bin/clusterscan.py"

    '''
    clusterinitsites.py -p !{commonPeaks} --clusterScan !{clusterScanScript} -o !{filePrefix}
    '''

  }

} else {
  process clusterInitiationSitesFromPeaks {

    tag { params.filePrefix }

    publishDir  path: "${params.outputDir}/peaks",
                mode: "copy",
                saveAs: { filename -> "${params.filePrefix}_IZ.bed"},
                pattern: "*_clusters.bed"

    input:
    file(commonPeaks) from resultsCallPeaks

    output:
    set val(params.filePrefix), file(commonPeaks), file("${params.filePrefix}_clusters.bed") into resultsCluster

    shell:
    clusterScanScript = "${NXF_HOME}/assets/pavrilab/hicer-nf/bin/clusterscan.py"

    '''
    clusterinitsites.py -p !{commonPeaks} --clusterScan !{clusterScanScript} -o !{params.filePrefix}
    '''

  }

}

process filterInitiationSites {

  tag { filePrefix }

  publishDir  path: "${params.outputDir}/peaks",
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

process multiqc {
    tag { 'all' }

    publishDir path: "${params.outputDir}",
               mode: 'copy',
               overwrite: 'true'

    input:
    file (fastqc: 'fastqc/*') from fastqcResults.collect()
    file (trim: 'trim/*') from trimgaloreResults.collect()
    file (bowtie: 'bowtie/*') from bowtieMultiqcChannel.collect()

    output:
    file "*multiqc_report.html" into multiqc_report

    script:
    """
    export LC_ALL=C.UTF-8
    export LANG=C.UTF-8
    multiqc -f -x *.run .
    """
}

workflow.onComplete {
	println ( workflow.success ? "COMPLETED!" : "FAILED" )

}
