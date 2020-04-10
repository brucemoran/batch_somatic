# Somatic-Germline n-of-1 DNAseq Analysis
## Call SNV and CNA from somatic and matched germline sample data
### How to Setup
#### Dependencies:
[NextFlow](https://www.nextflow.io/index.html#GetStarted), [Singularity](https://sylabs.io/guides/3.0/user-guide/installation.html#)
#### Reference Generation:
See [DNAseq_references repo](https://github.com/brucemoran/DNAseq_references/tree/master/GRCh_37-38)
### Somatic-Germline n-of-1 Pipeline
#### About the pipeline:
This pipeline was developed to analyse and report on clinical research. To this end we have tried to make a useful and readable output. This has been achieved largely on the back of [PCGR](https://github.com/sigven/pcgr)/[CPSR](https://github.com/sigven/cpsr) which provide really excellent HTML reports and annotation from multiple clinically relevant sources.

Variant calling uses an 'ensemble' approach with 3 callers (MuTect2, Manta/Strelka2 and Lancet) currently. More may be added in the future. The outputs are combined using our [somaticVariantConsensus](https://github.com/brucemoran/somaticVariantConsensus) method. This takes parsed VCF output and combines calls, requiring support from at least 2 callers. It also parses 'raw' (unfiltered) calls, testing those variants that do not have support from 2 callers to see if they are contained in any single callers filtered VCF. In this way we attempt to retain as much true-positive calls as possible.
#### To run the pipeline:
```
nextflow run brucemoran/somatic_n-of-1 -profile standard,singularity
```
N.B. that currently Manta/Strelka2 and Lancet are not available through Conda. Because of this Singularity is required to run the pipeline. This may change and we will update, but Singularity is good so we recommend it anyway!
#### Arguments:
```
nextflow run brucemoran/somatic_n-of-1 --help
```
#### Formats of --sampleCsv:
```
type,sampleID,meta,read1,read2
germline,germ1,whole_blood,/full/path/to/germ1.R1.fastq.gz,/full/path/to/germ1.R2.fastq.gz
somatic,soma1,primary_tumour,/full/path/to/soma1.R1.fastq.gz,/full/path/to/soma1.R2.fastq.gz
somatic,soma2,metastasis,/full/path/to/soma2.R1.fastq.gz,/full/path/to/soma2.R2.fastq.gz
```
The `meta` column is used for reporting in PCGR, CPSR where `sampleID` may include clinical/personal data, for example.

Headers of `sample.csv` file must match above exactly, and you should have only one germline/normal sample per run.

Column `type` must be `germline` for one sample only. In case of 2+ germline samples, run with each as `germline` in turn, specifying others as `somatic`.
