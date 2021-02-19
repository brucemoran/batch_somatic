# Batch Somatic-Germline DNAseq Analysis
## Call SNV and CNA from batches of somatic and matched germline sample data
### How to Setup
#### Dependencies:
[NextFlow](https://www.nextflow.io/index.html#GetStarted), [Singularity](https://sylabs.io/guides/3.0/user-guide/installation.html#)
#### Reference Generation:
See [somatic_n-of-1 repo's download-references.nf](https://github.com/brucemoran/somatic_n-of-1/blob/master/download-references.nf)
### Batch Somatic-Germline Pipeline
#### About the pipeline:
This pipeline was developed to analyse and report on clinical research. To this end we have tried to make a useful and readable output. This has been achieved largely on the back of [PCGR](https://github.com/sigven/pcgr)/[CPSR](https://github.com/sigven/cpsr) which provide really excellent HTML reports and annotation from multiple clinically relevant sources.

Variant calling uses an 'ensemble' approach with 3 callers (MuTect2, Manta/Strelka2 and Lancet) currently. More may be added in the future. The outputs are combined using our [somaticVariantConsensus](https://github.com/brucemoran/somaticVariantConsensus) method. This takes parsed VCF output and combines calls, requiring support from at least 2 callers. It also parses 'raw' (unfiltered) calls, testing those variants that do not have support from 2 callers to see if they are contained in any single callers filtered VCF. In this way we attempt to retain as much true-positive calls as possible.
#### To run the pipeline:
```
nextflow run brucemoran/batch_somatic -profile standard,singularity
```
N.B. that currently Manta/Strelka2 and Lancet are not available through Conda. Because of this Singularity is required to run the pipeline. This may change and we will update, but Singularity is good so we recommend it anyway!
#### Arguments:
```
nextflow run brucemoran/batch_somatic --help
```
#### Formats of --sampleCsv:
```
caseID,soma_sampleID,soma_read1,soma_read2,germ_sampleID,germ_read1,germ_read2
case_001,soma_001,/full/path/to/soma_001.R1.fastq.gz,/full/path/to/soma_001.R2.fastq.gz,germ_001,/full/path/to/germ_001.R1.fastq.gz,/full/path/to/germ_001.R2.fastq.gz
case_002,soma_002,/full/path/to/soma_002.R1.fastq.gz,/full/path/to/soma_002.R2.fastq.gz,germ_002,/full/path/to/germ_002.R1.fastq.gz,/full/path/to/germ_002.R2.fastq.gz
```
Headers of `sample.csv` file must match above exactly, and you should have only one somatic/germline pair of samples specified per line.
