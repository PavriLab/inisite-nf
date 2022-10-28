# inisite-nf

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A519.10.0-brightgreen.svg)](https://www.nextflow.io/) [![DOI](https://zenodo.org/badge/249460546.svg)](https://zenodo.org/badge/latestdoi/249460546)


## Introduction
**inisite-nf** is a bioinformatics analysis pipeline used for mapping initiation sites from nascent strand sequencing data (NS-seq).

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner.

## Pipeline summary
The basic principle of the mapping workflow is abstracted from [Cayrou et al, Genome Research 2015](http://genome.cshlp.org/content/25/12/1873). In brief, we use the MACS peak caller to map NS-seq peaks genome-wide. There are two modes the input can be processed, dependent on if a second NS-seq data set is given or not.

**a.  Single NS-seq data**
1.  Adapter and quality trimming with trim_galore ([trim_galore](https://github.com/FelixKrueger/TrimGalore))
2.  Alignment of reads with bowtie ([bowtie](http://bowtie-bio.sourceforge.net/manual.shtml))
3.  Calling narrow peaks with or without control using MACS2 ([MACS2](https://github.com/taoliu/MACS))
4.  Spatial clustering of MACS peaks with ClusterScan ([ClusterScan](https://github.com/pyrevo/ClusterScan))
5.  Filtering MACS peaks by overlap with identified peak clusters with BEDTools ([BEDTools](https://bedtools.readthedocs.io/en/latest/))

**b.  Dual NS-seq data**
1.  Adapter and quality trimming with trim_galore ([trim_galore](https://github.com/FelixKrueger/TrimGalore))
2.  Alignment of reads with bowtie ([bowtie](http://bowtie-bio.sourceforge.net/manual.shtml))
3.  Calling narrow peaks for both datasets with or without control using MACS2 ([MACS2](https://github.com/taoliu/MACS))
4.  Identifying and retaining peaks common to both datasets with BEDTools ([BEDTools](https://bedtools.readthedocs.io/en/latest/))
5.  Spatial clustering of common peaks with ClusterScan ([ClusterScan](https://github.com/pyrevo/ClusterScan))
6.  Filtering common peaks by overlap with identified peak clusters with BEDTools ([BEDTools](https://bedtools.readthedocs.io/en/latest/))

## Quick Start
i. Install [`nextflow`](https://nf-co.re/usage/installation)

ii. Clone repository with 
```bash
nextflow pull pavrilab/inisite-nf
```

iii. Start running your own analysis!

**a. Single**
```bash
nextflow run pavrilab/iniseq-nf --treatment ns_seq1.fastq [--control control1.fastq] --genome mm9
```

**b. Dual**
```bash
nextflow run pavrilab/iniseq-nf --treatment ns_seq1.fastq --treatment2 ns_seq2.fastq [--control control1.fastq --control2 control2.fastq] --genome mm9
```

Note that we only support human, mouse, Drosophila and C. elegans by default (see `conf/igenomes.conf`). If you want to use a different genome either specify `--fasta`/`--bowtieIndex` and `--genomeSize` or add the respective genome to the igenomes.conf file.

## Main arguments
#### `-profile`
Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments. For example `-profile cbe` invokes the execution of processes using the [`slurm`](https://slurm.schedmd.com/documentation.html) workload manager. If no profile is given the pipeline will be executed locally.

#### `--genome`

The genome argument has two properties. Firstly, it is used to document the name of the reference genome of the organism the Hi-C data originates from and secondly it can be used to retrieve prespecified reference data from a local igenomes database, where the pipeline automatically takes the files it requires for processing the Hi-C data (i.e. a bowtie2 index, a genome fasta and a chromSizes file).

There are 31 different species supported in the iGenomes references. To run the pipeline, you must specify which to use with the `--genome` flag.

**Note that here we only support genomes that have precompiled effective genome sizes in MACS. All other genomes have to be added by the user**

You can find the keys to specify the genomes in the [iGenomes config file](../conf/igenomes.config). Common genomes that are supported are:

* Human
  * `--genome hg19`
* Mouse
  * `--genome mm9`
* _Drosophila_
  * `--genome dm6`

Note that you can use the same configuration setup to save sets of reference files for your own use, even if they are not part of the iGenomes resource. See the [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for instructions on where to save such a file.

The syntax for this reference configuration is as follows:

```nextflow
params {
  genomes {
    'GRCh37' {
      fasta      = "/path/to/genome/fasta/file" // Used if no bowtie index
      bowtie     = "/path/to/bowtie/index/basename
      genomeSize = effective size of the genome in bp to use for MACS (see MACS documentation)
    }
    // Any number of additional genomes, key is used with --genome
  }
}
```

#### `--treatment`
Raw NS-seq reads in FASTQ format

#### `--treatment2`
Optional second set of raw NS-seq reads in FASTQ format

#### `--control`
raw reads of background signal (e.g. sheared genomic DNA). If not given MACS calls peak without input.

#### `--control2`
raw control reads for `--treatment2`. If `--control` was specified this also must be specified. If both control arguments are not given MACS calls peaks without input for both.

## Generic arguments
#### `--genome`
Genome from which the NS-seq reads originate (same as `-g` option in MACS). Has to be specified.

#### `--fasta`
This parameter is used to specify the genome fasta file and is only required if no igenomes database is available or bowtie index is not available locally. The file is used for bowtie index computation (if not specified manually)

#### `--genomeSize`
effective size of the genome to use with MACS if not specified in igenomes

#### `--bowtieIndex`
index to use with bowtie. Is generated from fasta file if not specified in igenomes

#### `--extensionSize`
Integer value specifying the length to which each read is extended by MACS before peak calling

#### `--qValueCutoff`
Float value between 0 and 1 specifying the q-value cutoff used by MACS to identify significant pileups

#### `--filePrefix`
Prefix for the result files name

#### `--outputDir`
Folder to which results will be written (is created if not existing)

## Results
The pipeline's results comprise 3 to 5 files depending on the number of treatment files given (i.e. single or dual mode). If the pipeline is run in single mode the `--outputDir` will contain 3 files:

1.  aligned reads in BAM format
2.  A `*_MACS.bed` file containing the original peak called by macs2, named according to the treatment BAM basename
3.  A `*_IS.bed` file containing the cluster-filtered MACS peaks
4.  A `*_IZ.bed` file containing the clusters called by ClusterScan and used to filter the MACS peaks

If the pipeline is run in dual mode the `--outputDir` will contain 5 files:

1.  aligned reads in BAM format
2.  Two `*_MACS.bed` file containing the original peaks called by macs2 for each treatment BAM, named according to the treatment BAM basename
3.  A `*.common.bed` containing all MACS peaks found in both treatment BAMs
4.  A `*_IS.bed` file containing the cluster-filtered common MACS peaks
5.  A `*_IZ.bed` file containing the clusters called by ClusterScan and used to filter the common MACS peaks

## Credits
The pipeline was developed by [Daniel Malzl](mailto:daniel.malzl@gmx.at) for use at the [IMP](https://www.imp.ac.at/), Vienna.

Many thanks to others who have helped out along the way too, including (but not limited to): [@t-neumann](https://github.com/t-neumann), [@pditommaso](https://github.com/pditommaso).

## Citations
### Pipeline tools
* [Nextflow](https://www.ncbi.nlm.nih.gov/pubmed/28398311/)
  > Di Tommaso P, Chatzou M, Floden EW, Barja PP, Palumbo E, Notredame C. Nextflow enables reproducible computational workflows. Nat Biotechnol. 2017 Apr 11;35(4):316-319. doi: 10.1038/nbt.3820. PubMed PMID: 28398311.

* [BEDTools](https://www.ncbi.nlm.nih.gov/pubmed/20110278/)
  > Quinlan AR, Hall IM. BEDTools: a flexible suite of utilities for comparing genomic features. Bioinformatics. 2010 Mar 15;26(6):841-2. doi: 10.1093/bioinformatics/btq033. Epub 2010 Jan 28. PubMed PMID: 20110278. PubMed Central PMCID: PMC2832824.
  
* [Bowtie](https://pubmed.ncbi.nlm.nih.gov/19261174/)
  > B. Langmead, C. Trapnell, M. Pop, S. L. Salzberg, Ultrafast and memory-efficient alignment of short DNA sequences to the human genome. Genome Biology 10, R25 (2009). doi: 10.1186/gb-2009-10-3-r25. PubMed PMID: 19261174 PMCID: PMC2690996
  
* [cutadapt](http://journal.embnet.org/index.php/embnetjournal/article/view/200)
  > M. Martin, Cutadapt removes adapter sequences from high-throughput sequencing reads. 2011 17, 3 (2011). doi: 10.14806/ej.17.1.200
  
* [ClusterScan](https://www.ncbi.nlm.nih.gov/pubmed/29912285)
  > M. Volpe, M. Miralto, S. Gustincich, R. Sanges, ClusterScan: simple and generalistic identification of genomic clusters. Bioinformatics. 2018 Nov 15;34(22):3921-3923. doi: 10.1093/bioinformatics/bty486. PubMed PMID: 29912285
  
* [MACS](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2592715/)
  > Zhang Y, Liu T, Meyer CA, Eeckhoute J, Johnson DS, Bernstein BE, Nusbaum C, Myers RM, Brown M, Li W, Liu XS. Model-based analysis of ChIP-Seq (MACS). Genome Biol. 2008;9(9):R137. doi: 10.1186/gb-2008-9-9-r137. Epub 2008 Sep 17. PubMed PMID: 18798982. PubMed Central PMCID: PMC2592715
  
### Python Packages
* [pandas](https://pandas.pydata.org/docs/index.html)
  > Wes McKinney. Data Structures for Statistical Computing in Python, Proceedings of the 9th Python in Science Conference, 51-56 (2010)
