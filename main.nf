#!/usr/bin/env nextflow

def helpMessage() {
  log.info"""
  -----------------------------------------------------------------------
                          SOMATIC_BATCH PIPELINE
  -----------------------------------------------------------------------
  Usage:

  nextflow run brucemoran/batch_somatic

  Mandatory arguments:
    -profile        [str]       Configuration profile (required: standard,singularity)
    --sampleCsv      STRING      CSV format, headers: caseID, soma_sampleID, soma_read1, soma_read2, germ_sampleID, germ_read1, germ_read2'
    --runID         [str]       Name for run, used to tag outputs
    --refDir        [file]      Path of dir in which reference data are held; this should be created by download-references.nf and contain dir <assembly>
    --assembly      [str]       Either GRCh37 or GRCh38 (default), as per download-references.nf

  General Optional Arguments:
    --germline      [bool]      Run HaplotypeCaller on germline sample and annotate with CPSR (default: true)
    --scatGath      [int]       Number of pieces to divide intervalList into for scattering to variant calling processes (default: 20 for exome, 100 for WGS)
    --incOrder      [str]       In final plots, use this ordering of samples (if multiple somatic samples); comma-separated, no spaces (default: alphanumeric sort)4
    --multiqcConfig [str]       Config file for multiqc (default: bin/somatic_n-of-1.multiQC_config.yaml)
    --seqLevel      [str]       WGS, exome or panel (default: WGS)
    --exomeTag      [str]       Tag used for exome (or panel) kit when running download-references.nf
    --cosmic        [bool]      set this to specify output of COSMIC CGC genes only (somatic only; based on download and supply of CGC file in download_references.nf)
    --impacts       [str]       'IMPACTS from VEP to use, comma-separated one or set of HIGH, MODERATE, MODIFIER, LOW; multiple sets can be input separating by semi-colon (default: "HIGH,MODERATE,MODIFIER,LOW")'
    --agilentUMI   [bool]      set this to process reads using a UMI supplied in sampleCsv by adding 'soma_umi' and 'germ_umi' after 'soma', 'germ' read paths and supplying path to UMI fastq ; NB this bypasses bbduk, picard markDuplicates
    --agilentLib   [str]      which of 'halo, (haloplex)hs, (sureselect)xt (default), (sureselect xt hs2) v2, (sureselect)qxt' were used to design and prep libraries
    """.stripIndet()
}

if (params.help) exit 0, helpMessage()

//Test Mandatory Arguments
if(!Channel.from(params.sampleCsv, checkIfExists: true)){
  exit 1, "Please include --sampleCsv, see --help for format"
}

if(!Channel.from(params.runID, checkIfExists: true)){
    exit 1, "Please include --runID <your_runID>"
}

if(!Channel.from(params.refDir, checkIfExists: true)){
  exit 1, "Please include --refDir <path> see github.com/brucemoran/somatic_n-of-1/ for how to run download-references.nf"
}

if(!Channel.from(params.assembly, checkIfExists: true)){
    exit 1, "Please include --assembly <GRCh3x>"
}

// if(!params.email){
//     exit 1, "Please include --email your@email.com"
// }

//Globals
params.outDir = "batch_${params.seqLevel}/${params.runID}"
params.seqlevel = "${params.seqLevel}".toLowerCase()

//Java task memory allocation via task.memory
javaTaskmem = { it.replace(" GB", "g") }

//Reference data as value channels and reusable therefore
reference = [
    grchvers: false,
    fa: false,
    fai: false,
    dict: false,
    bwa: false,
    hc_dbs: false,
    dbsnp: false,
    gridss: false,
    pcgrbase: false,
    intlist: false,
    seqlevel: false,
    bbres: false
]

reference.grchvers  = Channel.fromPath("${params.refDir}/${params.assembly}/pcgr/data/*", type: 'dir').getVal()
reference.fa = Channel.value(file(params.genomes[params.assembly].fa))
reference.fai = Channel.value(file(params.genomes[params.assembly].fai))
reference.dict = Channel.value(file(params.genomes[params.assembly].dict))
reference.bwa = Channel.value(file(params.genomes[params.assembly].bwa))
reference.hc_dbs = Channel.value(file(params.genomes[params.assembly].hc_dbs))
reference.dbsnp = Channel.value(file(params.genomes[params.assembly].dbsnp))
reference.gridss = Channel.value(file(params.genomes[params.assembly].gridss))
reference.pcgrbase = Channel.value(file(params.genomes[params.assembly].pcgr))

//change panel to exome as we haven't drawn distinction before
//also add exomeTag to specify dir for files input
reference.seqlevel = Channel.value(file(params.genomes[params.assembly].exome))

//set cosmic
reference.cosmic = params.cosmic == true ? Channel.value(file(params.genomes[params.assembly].cosmic)) : null

//setting of intlist based on seqlevel and exomeTag
reference.intlist = params.seqlevel == "wgs" ? Channel.fromPath("${params.refDir}/${params.assembly}/${params.seqlevel}/wgs.bed.interval_list").getVal() : Channel.fromPath("${params.refDir}/${params.assembly}/exome/${params.exomeTag}.bed.interval_list").getVal()

//Java task memory allocation via task.memory
javaTaskmem = { it.replace(" GB", "g") }

/* -0.02: Input using sample.csv
*/
if(! params.agilentUMI){
  Channel.fromPath("$params.sampleCsv", type: 'file')
         .splitCsv( header: true )
         .map { row -> [row.caseID, row.soma_sampleID, file(row.soma_read1), file(row.soma_read2), row.germ_sampleID, file(row.germ_read1), file(row.germ_read2)] }
         .set { split_soma_germ }

  // -0.011: Input using sample.csv
  process splt_sg {

    label 'low_mem'

    publishDir "$params.outDir/cases/$caseID", mode: "copy", pattern: "*.csv"

    input:
    tuple val(caseID), val(soma_sampleID), file(soma_read1), file(soma_read2), val(germ_sampleID), file(germ_read1), file(germ_read2) from split_soma_germ

    output:
    file('*.csv') into splitcsv2

    when:
    !params.agilentUMI

    """
    SR1=\$(readlink -e ${soma_read1})
    SR2=\$(readlink -e ${soma_read2})
    GR1=\$(readlink -e ${germ_read1})
    GR2=\$(readlink -e ${germ_read2})
    echo "caseID,type,sampleID,read1,read2" > $caseID".csv"
    echo "${caseID},somatic,${soma_sampleID},\$SR1,\$SR2" >> ${caseID}".csv"
    echo "${caseID},germline,${germ_sampleID},\$GR1,\$GR2" >> ${caseID}".csv"
    """
  }

  splitcsv2.splitCsv( header: true )
           .map { row -> [row.caseID, row.type, row.sampleID, file(row.read1), file(row.read2)] }
           .set { bbduking }
}

if(params.agilentUMI){
  Channel.fromPath("$params.sampleCsv", type: 'file')
         .splitCsv( header: true )
         .map { row -> [row.caseID, row.soma_sampleID, file(row.soma_read1), file(row.soma_read2), file(row.soma_umi), row.germ_sampleID, file(row.germ_read1), file(row.germ_read2), file(row.germ_umi)] }
         .set { split_soma_germ_umi }

  // -0.012: Input using sample.csv with UMIs
  process splt_sg_umi {

    label 'low_mem'

    publishDir "$params.outDir/cases/$caseID", mode: "copy", pattern: "*.csv"

    input:
    tuple val(caseID), val(soma_sampleID), file(soma_read1), file(soma_read2), file(soma_umi), val(germ_sampleID), file(germ_read1), file(germ_read2), file(germ_umi) from split_soma_germ_umi

    output:
    file('*.csv') into splitcsv2_umi

    when:
    params.agilentUMI

    """
    SR1=\$(readlink -e ${soma_read1})
    SR2=\$(readlink -e ${soma_read2})
    SMI=\$(readlink -e ${soma_umi})
    GR1=\$(readlink -e ${germ_read1})
    GR2=\$(readlink -e ${germ_read2})
    GMI=\$(readlink -e ${germ_umi})
    echo "caseID,type,sampleID,read1,read2,umi" > $caseID".csv"
    echo "${caseID},somatic,${soma_sampleID},\$SR1,\$SR2,\$SMI" >> ${caseID}".csv"
    echo "${caseID},germline,${germ_sampleID},\$GR1,\$GR2,\$GMI" >> ${caseID}".csv"
    """
  }

  splitcsv2_umi.splitCsv( header: true )
               .map { row -> [row.caseID, row.type, row.sampleID, file(row.read1), file(row.read2), file(row.umi)] }
               .set { trimmer }

/*
================================================================================
                          0. PREPROCESS INPUT FASTQ
================================================================================
*/
  // 0.01: Input trimming
  process trimmer {

    label 'med_mem'

    publishDir path: "$params.outDir/cases/$caseID/bbduk", mode: "copy", pattern: "${sampleID}.trimmer.log.txt"

    input:
    tuple val(caseID), val(type), val(sampleID), file(read1), file(read2), file(umi) from trimmer

    output:
    file("${sampleID}.trimmer.log.txt") into log_trimmer
    tuple val(caseID), val(type), val(sampleID), file("${sampleID}.trimmer.R1.fastq.gz"), file("${sampleID}.trimmer.R2.fastq.gz") into bbduking
    tuple val(caseID), val(type), val(sampleID), file(umi) into locatiting

    script:
    def taskmem = task.memory == null ? "" : "-Xmx" + javaTaskmem("${task.memory}")
    """
    {
      java ${taskmem} \
        -jar /usr/share/agent/lib/trimmer-2.0.2.jar \
        -${params.agilentLib} \
        -out_loc ./ \
        -fq1 ${read1} \
        -fq2 ${read2}
      PR1=\$(ls | grep ${sampleID} | grep R1 | grep Cut)
      PR2=\$(ls | grep ${sampleID} | grep R3 | grep Cut)
      mv \$PR1 ${sampleID}.trimmer.R1.fastq.gz
      mv \$PR2 ${sampleID}.trimmer.R2.fastq.gz
    } 2>&1 | tee > ${sampleID}.trimmer.log.txt
    """
  }
}

// 0.1: Input trimming
process bbduk {

  label 'med_mem'

  publishDir path: "$params.outDir/cases/$caseID/bbduk", mode: "copy", pattern: "*.txt"

  input:
  tuple val(caseID), val(type), val(sampleID), file(read1), file(read2) from bbduking

  output:
  file('*.txt') into log_bbduk
  tuple val(caseID), val(type), val(sampleID), file("${sampleID}.bbduk.R1.fastq.gz"), file("${sampleID}.bbduk.R2.fastq.gz") into bwa_memming
  tuple val(caseID), val(type), val(sampleID), file(read1), file(read2) into fastping

  script:
  def taskmem = task.memory == null ? "" : "-Xmx" + javaTaskmem("${task.memory}")
  """
  {
  ##remove because of issues on readname including "1:xyz, 2:xyz for R1, R2"
  if [[ ${read1} =~ ".gz"\$ ]]; then
    gunzip -c ${read1} | perl -ane 'chomp; print "\$F[0]\\n";' > r1.fq
    gunzip -c ${read2} | perl -ane 'chomp; print "\$F[0]\\n";' > r2.fq
  else
    perl -ane 'chomp; print "\$F[0]\\n";' ${read1} > r1.fq
    perl -ane 'chomp; print "\$F[0]\\n";' ${read2} > r2.fq
  fi

  ##repair in case of disorder
  repair.sh in1=r1.fq in2=r2.fq out1=r1.f.fq out2=r2.f.fq repair

  sh bbduk.sh ${taskmem} \
    in1=r1.f.fq \
    in2=r2.f.fq \
    out1=${sampleID}".bbduk.R1.fastq.gz" \
    out2=${sampleID}".bbduk.R2.fastq.gz" \
    k=31 \
    mink=5 \
    hdist=1 \
    ktrim=r \
    trimq=20 \
    qtrim=rl \
    maq=20 \
    ref=/opt/miniconda/envs/somatic_n-of-1/opt/bbmap-adapters.fa \
    tpe \
    tbo \
    stats=${sampleID}".bbduk.adapterstats.txt" \
    overwrite=T
  } 2>&1 | tee > ${sampleID}.bbduk.runstats.txt
  rm r*.fq
  """
}

// 0.2: fastp QC of pre-, post-bbduk
process fastp {

  label 'low_mem'

  publishDir "$params.outDir/cases/$caseID/fastp", mode: "copy", pattern: "*.html"

  input:
  tuple val(caseID), val(type), val(sampleID), file(preread1), file(preread2) from fastping

  output:
  file('*.html') into fastp_html
  file('*.json') into fastp_multiqc

  script:
  """
  fastp -w ${task.cpus} -h ${sampleID}".fastp.html" -j ${sampleID}".fastp.json" --in1 ${preread1} --in2 ${preread2}
  """
}

// 1.0: Input alignment
process bwamem {

  label 'high_mem'

  input:
  tuple val(caseID), val(type), val(sampleID), file(read1), file(read2) from bwa_memming
  file(bwa) from reference.bwa

  output:
  tuple val(caseID), val(type), val(sampleID), file("${sampleID}.sort.bam"), file('*.bai') into (cramming, dup_marking)

  script:
  def fa = "${bwa}/*fasta"
  """
  DATE=\$(date +"%Y-%m-%dT%T")
  RGLINE="@RG\\tID:${sampleID}\\tPL:ILLUMINA\\tSM:${sampleID}\\tDS:${type}\\tCN:UCD\\tLB:LANE_X\\tDT:\$DATE"

  bwa mem \
    -t${task.cpus} \
    -M \
    -R \$RGLINE \
    ${fa} \
    ${read1} ${read2} | \
    samtools sort -T "tmp."${sampleID} -o ${sampleID}".sort.bam"
  samtools index ${sampleID}".sort.bam"
  """
}

// 1.1: CRAM alignment and output
process cram {

  label 'low_mem'
  publishDir path: "$params.outDir/cases/$caseID/bwa", mode: "copy", pattern: "*.cra*"

  input:
  tuple val(caseID), val(type), val(sampleID), file(bam), file(bai) from cramming
  file(bwa) from reference.bwa

  output:
  tuple file("${sampleID}.sort.cram"), file("${sampleID}.sort.cram.crai") into completedcram

  script:
  """
  samtools view -hC -T ${bwa}/*fasta ${sampleID}".sort.bam" > ${sampleID}".sort.cram"
  samtools index ${sampleID}".sort.cram"
  """
}

// 1.2: MarkDuplicates
//either mrkdup or locatit
if(! params.agilentUMI){
  process mrkdup {

    label 'high_mem'

    publishDir path: "$params.outDir/cases/$caseID/picard", mode: "copy", pattern: "*.txt"

    input:
    tuple val(caseID), val(type), val(sampleID), file(bam), file(bai) from dup_marking

    output:
    file('*.txt') into mrkdup_multiqc
    tuple val(caseID), val(type), val(sampleID), file('*.md.bam'), file('*.md.bam.bai') into gatk4recaling

    script:
    def taskmem = task.memory == null ? "" : "-Xmx" + javaTaskmem("${task.memory}")
    """
    OUTBAM=\$(echo ${bam} | sed 's/bam/md.bam/')
    OUTMET=\$(echo ${bam} | sed 's/bam/md.metrics.txt/')
    {
    picard ${taskmem} \
      MarkDuplicates \
      TMP_DIR=./ \
      INPUT=${bam} \
      OUTPUT=/dev/stdout \
      COMPRESSION_LEVEL=0 \
      QUIET=TRUE \
      METRICS_FILE=\$OUTMET \
      REMOVE_DUPLICATES=FALSE \
      ASSUME_SORTED=TRUE \
      VALIDATION_STRINGENCY=LENIENT \
      VERBOSITY=ERROR | samtools view -Shb - > \$OUTBAM

    samtools index \$OUTBAM
    } 2>&1 | tee > ${sampleID}.picard_markDuplicates.log.txt
    """
  }
}
if(params.agilentUMI){
  locatiting
    .map { it -> [it[0],it[1],it[2],it[3]] }
    .set { locatitjoin }
  dup_marking
    .map { it -> [it[0],it[1],it[2],it[3],it[4]] }
    .join(locatitjoin, by: [0,1,2])
    .set { locatitumi }
  process locatit {

    label 'high_mem'

    publishDir path: "$params.outDir/cases/$caseID/picard", mode: "copy", pattern: "*.txt"

    input:
    tuple val(caseID), val(type), val(sampleID), file(bam), file(bai), file(umi) from locatitumi
    file(exome) from reference.seqlevel

    output:
    tuple val(caseID), val(type), val(sampleID), file('*.md.bam'), file('*.md.bam.bai') into gatk4recaling
    file("${sampleID}.locatit.md.metrics.txt") into locatit_flagsqc

    script:
    def taskmem = task.memory == null ? "" : "-Xmx" + javaTaskmem("${task.memory}")
    def bedfile = "${exome}/${params.exomeTag}/${params.exomeTag}.bed"
    """
    {
      java ${taskmem} \
        -jar /usr/share/agent/lib/locatit-2.0.2.jar \
        -l $bedfile \
        -U \
        -i \
        -IB \
        -OB \
        -o ${sampleID}.md.bam \
        ${bam} \
        ${umi}
    } 2>&1 | tee > ${sampleID}.umi.log.txt

    ##allow above to fail (couldn't diagnose why it does...)
    if [[ ! -e ${sampleID}.md.bam ]]; then
      cp ${bam} ${sampleID}.md.bam
    fi

    samtools sort -T tmp.${sampleID} -o ${sampleID}.md.sort.bam ${sampleID}.md.bam
    mv ${sampleID}.md.sort.bam ${sampleID}.md.bam
    samtools index ${sampleID}.md.bam

    samtools flagstat ${bam} > ${sampleID}_flagstat_pre.txt
    samtools flagstat ${sampleID}.md.bam > ${sampleID}.flagstat_post.md.txt
    ##convert those to picard dup metric format

    RP=\$(grep "in total" ${sampleID}.flagstat_post.md.txt | perl -ane '\$t=\$F[0]+\$F[2]; print \$t/2;' | cut -d "." -f1)
    SP=\$(grep "secondary" ${sampleID}.flagstat_post.md.txt | perl -ane '\$t=\$F[0]+\$F[2]; print \$t/2;' | cut -d "." -f1)
    UM=\$(grep "singletons" ${sampleID}.flagstat_post.md.txt | perl -ane '\$t=\$F[0]+\$F[2]; print \$t/2;' | cut -d "." -f1)
    UPD=\$(grep "duplicates" ${sampleID}.flagstat_post.md.txt | perl -ane '\$t=\$F[2]; print \$t/2;' | cut -d "." -f1)
    PD=\$(grep "duplicates" ${sampleID}.flagstat_post.md.txt | perl -ane '\$t=\$F[0]; print \$t/2;' | cut -d "." -f1)
    RPO=0
    export UP PD RP
    if [[ \$PD == 0 ]]; then
      PC=0
    else
      PC=\$(perl -e '\$p=\$ENV{UP}+\$ENV{PD};print \$p/\$ENV{RP};')
    fi
    echo -e "## htsjdk.samtools.metrics.StringHeader\\n# MarkDuplicates INPUT=[${bam}] OUTPUT=/dev/stdout METRICS_FILE=${sampleID}.locatit.md.metrics.txt REMOVE_DUPLICATES=false ASSUME_SORTED=true TMP_DIR=[.] VERBOSITY=ERROR QUIET=true VALIDATION_STRINGENCY=LENIENT COMPRESSION_LEVEL=0    MAX_SEQUENCES_FOR_DISK_READ_ENDS_MAP=50000 MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=8000 SORTING_COLLECTION_SIZE_RATIO=0.25 TAG_DUPLICATE_SET_MEMBERS=false REMOVE_SEQUENCING_DUPLICATES=false TAGGING_POLICY=DontTag CLEAR_DT=true DUPLEX_UMI=false ADD_PG_TAG_TO_READS=true DUPLICATE_SCORING_STRATEGY=SUM_OF_BASE_QUALITIES PROGRAM_RECORD_ID=MarkDuplicates PROGRAM_GROUP_NAME=MarkDuplicates READ_NAME_REGEX=<optimized capture of last three ':' separated fields as numeric values> OPTICAL_DUPLICATE_PIXEL_DISTANCE=100 MAX_OPTICAL_DUPLICATE_SET_SIZE=300000 MAX_RECORDS_IN_RAM=500000 CREATE_INDEX=false CREATE_MD5_FILE=false GA4GH_CLIENT_SECRETS=client_secrets.json USE_JDK_DEFLATER=false USE_JDK_INFLATER=false\\n## htsjdk.samtools.metrics.StringHeader\\n# Started on: Sat Oct 03 10:49:00 UTC 2020\\n## METRICS CLASS\\tpicard.sam.DuplicationMetrics\\nLIBRARY\\tUNPAIRED_READS_EXAMINED\\tREAD_PAIRS_EXAMINED\\tSECONDARY_OR_SUPPLEMENTARY_RDS\\tUNMAPPED_READS\\tUNPAIRED_READ_DUPLICATES\\tREAD_PAIR_DUPLICATES\\tREAD_PAIR_OPTICAL_DUPLICATES\\tPERCENT_DUPLICATION\\tESTIMATED_LIBRARY_SIZE\\nLANE_X\\t0\\t\$RP\\t\$SP\\t\$UM\\t\$UPD\\t\$PD\\t\$RPO\\t\$PC\\t0" > ${sampleID}.locatit.md.metrics.txt
    """
  }
}
// 1.3: GATK4 BestPractices
process gtkrcl {

  label 'high_mem'

  publishDir path: "$params.outDir/cases/$caseID/gatk4/bqsr", mode: "copy", pattern: "*.GATK4_BQSR.log.txt"

  input:
  tuple val(caseID), val(type), val(sampleID), file(bam), file(bai) from gatk4recaling
  file(fasta) from reference.fa
  file(fai) from reference.fai
  file(dict) from reference.dict
  file(dbsnp_files) from reference.dbsnp
  file(intlist) from reference.intlist

  output:
  file('*.table') into gtkrcl_multiqc
  tuple val(caseID), val(type), val(sampleID), file('*.bqsr.bam'), file('*.bqsr.bam.bai') into ( germfiltering, gmultimetricing)
  tuple val(caseID), val(type), val(sampleID), file('*.bqsr.bam'), file('*.bqsr.bam.bai') into hc_germ
  tuple val(caseID), val(sampleID) into metas_pcgr
  file("${sampleID}.GATK4_BQSR.log.txt") into bqsr_log

  script:
  def dbsnp = "${dbsnp_files}/*gz"
  """
  {
  gatk BaseRecalibrator \
    -R ${fasta} \
    -I ${bam} \
    --known-sites \$(echo ${dbsnp}) \
    --use-original-qualities \
    -O ${sampleID}.recal_data.table \
    --disable-sequence-dictionary-validation true \
    -L ${intlist}

  #ApplyBQSR
  OUTBAM=\$(echo ${bam} | sed 's/bam/bqsr.bam/')
  gatk ApplyBQSR \
    -R ${fasta} \
    -I ${bam} \
    --bqsr-recal-file ${sampleID}.recal_data.table \
    --add-output-sam-program-record \
    --use-original-qualities \
    -O \$OUTBAM \
    -L ${intlist}

  samtools index \$OUTBAM \$OUTBAM".bai"
  } 2>&1 | tee >  ${sampleID}.GATK4_BQSR.log.txt
  """
}

// 1.31: scatter-gather implementation for mutect2, lancet
process scat_gath {

  label 'low_mem'

  input:
  file(intlist) from reference.intlist

  output:
  file('lancet.scatgath.*.bed') into lancet_bedding
  file('mutect2.scatgath.*.bed.interval_list') into mutect2_bedding
  file('hc.scatgath.*.bed.interval_list') into hc_bedding

  script:
  def sgcount = params.scatGath
  if (params.scatGath == null){
    if (params.seqlevel == "panel"){
      sgcount = 2
    }
    if (params.seqlevel == "exome"){
      sgcount = 20
    }
    if (params.seqlevel == "wgs"){
      sgcount = 100
    }
  }
  """
  ##strip out all but chromosomes in the interval_list (no decoys etc)
  CHRS=\$(grep -v "@" ${intlist} | cut -f 1 | uniq)
  for CHR in \$CHRS; do
    grep "SN:\$CHR\\s" ${intlist} >> used.interval_list
  done
  grep -v "@" ${intlist} >> used.interval_list

  ##generate scatters
  picard IntervalListTools \
    I=used.interval_list \
    SCATTER_COUNT=${sgcount} \
    O=\$(pwd)

  ##rename scatters and parse into appropriate format for tools
  ls temp*/* | while read FILE; do
    COUNTN=\$(dirname \$FILE | perl -ane '@s=split(/\\_/); print \$s[1];');
    mv \$FILE mutect2.scatgath.\${COUNTN}.bed.interval_list;
    cp mutect2.scatgath.\${COUNTN}.bed.interval_list hc.scatgath.\${COUNTN}.bed.interval_list
    grep -v @ mutect2.scatgath.\${COUNTN}.bed.interval_list | \
      cut -f 1,2,3,5 > lancet.scatgath.\${COUNTN}.bed
  done
  """
}

/*
================================================================================
                            2.  MUTATION CALLING
================================================================================
*/
// 2.0: GATK4 Germline Haplotypecaller
hcbedding = hc_bedding.flatten()
hc_germ
  .map { it -> [it[0],it[1],it[2],it[3],it[4]] }
  .combine(hcbedding)
  .set { hcgermbedding }

process haplotypecaller {

  label 'med_mem'

  input:
  tuple val(caseID), val(type), val(sampleID), file(bam), file(bai), file(intlist) from hcgermbedding
  file(fasta) from reference.fa
  file(fai) from reference.fai
  file(dict) from reference.dict
  file(dbsnp_files) from reference.dbsnp
  file(hc_dbs_files) from reference.hc_dbs

  output:
  tuple val(caseID), val(sampleID), file('*sort.hc.vcf') into hc_gt

  when:
  type == "germline" & params.germline != false

  script:
  def taskmem = task.memory == null ? "" : "--java-options \"-Xmx" + javaTaskmem("${task.memory}") + "\""
  def dbsnp = "${dbsnp_files}/*gz"
  def omni = "${hc_dbs_files}/KG_omni*.gz"
  def kgp1 = "${hc_dbs_files}/KG_phase1*.gz"
  def hpmp = "${hc_dbs_files}/hapmap*.gz"
  """
  SCATGATHN=\$(echo ${intlist} | perl -ane '@s=split(/\\./);print \$s[2];')
  gatk ${taskmem} HaplotypeCaller \
    -R ${fasta} \
    -I ${bam} \
    --dont-use-soft-clipped-bases \
    --standard-min-confidence-threshold-for-calling 20 \
    --dbsnp \$(echo ${dbsnp}) \
    --native-pair-hmm-threads ${task.cpus} \
    -O ${sampleID}".\${SCATGATHN}.hc.vcf" \
    --disable-sequence-dictionary-validation true \
    -L ${intlist}

  picard SortVcf \
    I=${sampleID}".\${SCATGATHN}.hc.vcf" \
    O=${sampleID}".\${SCATGATHN}.sort.hc.vcf" \
    SD=${dict}
  """
}

// 2.1: HC_merge
hc_gt
  .groupTuple()
  .map { it -> tuple(it[0], it[1][0], it[2][0..-1].flatten()) }
  .set { hc_fm }

process hc_merge {

  label 'high_mem'

  publishDir path: "$params.outDir/cases/$caseID/gatk4/haplotypecaller", mode: "copy", pattern: '*.hc.merge.vcf.gz'
  publishDir path: "$params.outDir/output/haplotypecaller", mode: "copy", pattern: '*.hc.merge.vcf.gz'

  input:
  tuple val(caseID), val(sampleID), file(rawvcfs) from hc_fm

  output:
  tuple val(caseID), val(sampleID), file("${sampleID}.hc.merge.vcf.gz"), file("${sampleID}.hc.merge.vcf.gz.tbi") into ( cpsr_vcf, vep_hc_vcf )

  script:
  """
  ls ${sampleID}.*.sort.hc.vcf > vcf.list
  picard MergeVcfs I=vcf.list O=${sampleID}".hc.merge.vcf"
  bgzip ${sampleID}".hc.merge.vcf"
  tabix ${sampleID}".hc.merge.vcf.gz"
  """
}

// 2.11: HC_merge VEP
process vepHC {

  label 'low_mem'

  publishDir path: "${params.outDir}/cases/${caseID}/vcf", mode: "copy"
  publishDir path: "${params.outDir}/output/haplotypecaller", mode: "copy"

  input:
  tuple val(caseID), val(sampleID), file(vcf), file(tbi) from vep_hc_vcf
  file(fasta) from reference.fa
  file(fai) from reference.fai
  file(dict) from reference.dict
  file(grchver) from reference.grchvers
  file(pcgrbase) from reference.pcgrbase

  output:
  file("${sampleID}.hc.merge.vep.vcf") into hc_vepd

  script:
  grch_vers = "${grchver}".split("\\/")[-1]
  """
  vep --dir_cache ${pcgrbase}/data/${grch_vers}/.vep \
    --offline \
    --assembly ${params.assembly} \
    --vcf_info_field ANN \
    --symbol \
    --species homo_sapiens \
    --check_existing \
    --cache \
    --fork ${task.cpus} \
    --af_1kg \
    --af_gnomad \
    --vcf \
    --input_file ${vcf} \
    --output_file ${sampleID}.hc.merge.vep.vcf \
    --format "vcf" \
    --fasta ${fasta} \
    --hgvs \
    --canonical \
    --ccds \
    --force_overwrite \
    --verbose
  """
}

// 2.12: VCF to CSV
process vepHCtsv {

  label 'low_mem'

  publishDir path: "${params.outDir}/output/haplotypecaller", mode: "copy"

  input:
  file(vcf_vep) from hc_vepd.collect()

  output:
  file("${params.runID}.haplotypecaller_all.tabvcf.tsv") into madetsv
  script:
  """
  perl ${workflow.projectDir}/bin/vepHCvcf_combine_tsv.pl "${params.runID}.haplotypecaller_all"
  """
}

// 2.2: CPSR annotation of GATK4 Germline
process cpsrreport {

  label 'low_mem'

  publishDir "$params.outDir/reports/cpsr", mode: "copy", pattern: "*.html"
  publishDir "$params.outDir/cases/$caseID/cpsr", mode: "copy"

  input:
  tuple val(caseID), val(sampleID), file(vcf), file(tbi) from cpsr_vcf
  file(grchver) from reference.grchvers
  file(pcgrbase) from reference.pcgrbase

  output:
  file('*') into cpsr_vcfs

  script:
  def grchv = "${grchver}".split("\\/")[-1]
  """
  {
  ##CPSR v0.6.0rc
  cpsr.py \
    --no-docker \
    --no_vcf_validate \
    --panel_id 0 \
    --query_vcf ${vcf} \
    --pcgr_dir ${pcgrbase} \
    --output_dir ./ \
    --genome_assembly ${grchv} \
    --conf ${pcgrbase}/data/${grchv}/cpsr_configuration_default.toml \
    --sample_id ${sampleID}

  } 2>&1 | tee > ${sampleID}.cpsr.log.txt
  """
}

// 2.3: filter germline channel, tap into somatic channels for all processes subsequent
def germfilter = branchCriteria {
                  germfiltered: it[1] == "germline"
                  return it

                  somafiltered: true
                  return it
                 }
germfiltering
    .branch(germfilter)
    .set { somagerm }

somagerm.somafiltered
    .map { [it[0], it[2..-1]] }
    .tap { somatap }

somagerm.germfiltered
    .map { [it[0], it[2..-1]] }
    .tap { germtap }

somatap.join(germtap).tap{ somagermtap }

//map caseID, soma_sampleID, soma_bam, soma_bai, germ_sampleID, germ_bam, germ_bai
somagermtap
    .map { it -> tuple(it[0],
                       it[1][0],
                       it[1][1..2],
                       it[2][0],
                       it[2][1..2]).flatten() }
    .into { mutect2somaticing; mutect2_contam; facetsomaing; mantastrelka2ing; lanceting }

// 2.4: PicardTools metrics suite for MultiQC HTML report
process mltmet {

  label 'med_mem'

  publishDir "$params.outDir/cases/$caseID/metrics", overwrite: 'true'

  input:
  tuple val(caseID), val(type), val(sampleID), file(bam), file(bai) from gmultimetricing
  file(fasta) from reference.fa
  file(fai) from reference.fai
  file(dict) from reference.dict
  file(intlist) from reference.intlist

  output:
  file('*.txt') into multimetrics_multiqc

  script:
  def taskmem = task.memory == null ? "" : "-Xmx" + javaTaskmem("${task.memory}")
  """
  {
  if [[ ${params.seqlevel} == "exome" || "panel" ]]; then
  picard ${taskmem} CollectHsMetrics \
    I=${bam} \
    O=${sampleID}".hs_metrics.txt" \
    TMP_DIR=./ \
    R=${fasta} \
    BAIT_INTERVALS=${intlist}  \
    TARGET_INTERVALS=${intlist}
  fi
  picard ${taskmem} CollectAlignmentSummaryMetrics \
    I=${bam} \
    O=${sampleID}".AlignmentSummaryMetrics.txt" \
    TMP_DIR=./ \
    R=${fasta}

  picard ${taskmem} CollectMultipleMetrics \
    I=${bam} \
    O=${sampleID}".CollectMultipleMetrics.txt" \
    TMP_DIR=./ \
    R=${fasta}

  picard ${taskmem} CollectSequencingArtifactMetrics \
    I=${bam} \
    O=${sampleID}".artifact_metrics.txt" \
    TMP_DIR=./ \
    R=${fasta}

  picard ${taskmem} CollectInsertSizeMetrics \
    I=${bam} \
    O=${sampleID}".insert_size_metrics.txt" \
    H=${bam}".histogram.pdf" \
    TMP_DIR=./

  } 2>&1 | tee > ${sampleID}.picard.metrics.log
  """
}

// 2.5: SCNA with facets CSV snp-pileup
process fctcsv {

  label 'med_mem'

  publishDir "$params.outDir/cases/$caseID/facets"

  input:
  tuple val(caseID), val(sampleID), file(tumourbam), file(tumourbai), val(germlineID), file(germlinebam), file(germlinebai) from facetsomaing
  file(dbsnp_files) from reference.dbsnp

  output:
  tuple file("${sampleID}.fit_cncf_jointsegs.tsv"), file("${sampleID}.fit_ploidy_purity.tsv") into facets_consensusing
  tuple val(sampleID), file("${sampleID}.cncf_jointsegs.pcgr.tsv"), file("${sampleID}.fit_ploidy_purity.pcgr.tsv") into facets_pcgr
  file("${sampleID}.facets.log.txt") into facets_log

  when:
  params.facets

  script:
  def dbsnp = "${dbsnp_files}/*gz"
  """
  { snp-pileup \
      \$(echo ${dbsnp}) \
      -r 10 \
      -p \
      ${sampleID}.facets.r10.csv \
      ${germlinebam} \
      ${tumourbam}

    Rscript -e "somenone::facets_cna_call(\\"${sampleID}.facets.r10.csv\\")"

    tail -n+2 ${sampleID}.fit_ploidy_purity.tsv > ${sampleID}.fit_ploidy_purity.pcgr.tsv
  } 2>&1 | tee > ${sampleID}.facets.log.txt
  """
}

// 2.6: SCNA consensus from facets
process fctcon {

  label 'med_mem'

  publishDir "$params.outDir/reports/scna/facets"

  input:
  file(filesn) from facets_consensusing.collect()
  file(cosmicbed) from reference.cosmic
  file(dict) from reference.dict

  output:
  file('*') into complete_facets

  when:
  params.facets

  script:
  if( !params.cosmic )
    """
    { Rscript -e "somenone::facets_cna_consensus(\\"fit_cncf_jointsegs.tsv\\", \\"${dict}\\", \\"${params.runID}\\")"
    } 2>&1 | tee > facets_cons.log.txt
    """
  else
    """
    { Rscript -e "somenone::facets_cna_consensus(\\"fit_cncf_jointsegs.tsv\\", \\"${dict}\\", \\"${params.runID}\\", \\"${cosmicbed}\\")"
    } 2>&1 | tee > facets_cons.log.txt
    """
}

mutect2bedding = mutect2_bedding.flatten()
mutect2somaticing
  .map { it -> [it[0],it[1],it[2],it[3],it[4],it[5],it[6]]}
  .combine(mutect2bedding)
  .set { mutect2somaticbedding }

// 2.71: MuTect2
// NB --germline-resource dollar-sign{dbsnp} removed as no AF causing error
process mutct2_sg {

  label 'med_mem'

  input:
  tuple val(caseID), val(sampleID), file(tumourbam), file(tumourbai), val(germlineID), file(germlinebam), file(germlinebai), file(intlist) from mutect2somaticbedding
  file(fasta) from reference.fa
  file(fai) from reference.fai
  file(dict) from reference.dict

  output:
  tuple val(caseID), val(sampleID), file('*.sort.mutect2.vcf') into mutect2_gt
  tuple val(caseID), val(sampleID), file('*.vcf.stats') into mutect2_st
  tuple val(caseID), val(sampleID), file('*.mutect2.f1r2.tar.gz') into mutect2_f1r2

  script:
  def taskmem = task.memory == null ? "" : "--java-options \"-Xmx" + javaTaskmem("${task.memory}") + "\""
  """
  SCATGATHN=\$(echo ${intlist} | perl -ane '@s=split(/\\./);print\$s[2];')
  gatk ${taskmem} \
    Mutect2 \
    --native-pair-hmm-threads ${task.cpus} \
    --reference ${fasta} \
    --input ${germlinebam} \
    --input ${tumourbam} \
    --normal-sample ${germlineID} \
    --tumor-sample ${sampleID} \
    --output ${sampleID}"."\${SCATGATHN}".mutect2.vcf" \
    --disable-sequence-dictionary-validation true \
    --f1r2-tar-gz ${sampleID}.\${SCATGATHN}".mutect2.f1r2.tar.gz" \
    -L ${intlist}

  picard SortVcf \
    I=${sampleID}"."\${SCATGATHN}".mutect2.vcf" \
    O=${sampleID}"."\${SCATGATHN}".sort.mutect2.vcf" \
    SD=${dict}
  """
}

// 2.72: MuTect2_merge
mutect2_gt
  .groupTuple()
  .map { it -> tuple(it[0], it[1][0], it[2][0..-1].flatten()) }
  .set { mutect2_fm }

process mutct2_concat {

  label 'med_mem'

  input:
  tuple val(caseID), val(sampleID), file(rawvcfs) from mutect2_fm

  output:
  tuple val(caseID), val(sampleID), file('*mutect2.merge.vcf') into mutect2_merge

  script:
  """
  ls *.sort.mutect2.vcf > vcf.list
  picard MergeVcfs I=vcf.list O=${sampleID}".mutect2.merge.vcf"
  """
}

mutect2_st
  .groupTuple()
  .map { it -> tuple(it[0], it[1][0], it[2][0..-1].flatten()) }
  .set { mutect2_sm }

// 2.73: MuTect2 Concatenate VCFs
process mutct2_concstat {

  label 'med_mem'

  input:
  tuple val(caseID), val(sampleID), file(stats) from mutect2_sm

  output:
  tuple val(caseID), val(sampleID), file('*mutect2.merge.vcf.stats') into mutect2_stats

  script:
  """
  STATS=\$(ls *stats | perl -ane 'foreach \$k (@F){print "--stats \$k ";}')
  gatk MergeMutectStats --output ${sampleID}".mutect2.merge.vcf.stats" \$STATS
  """
}

// 2.74: MuTect2 Concatenate VCFs
mutect2_f1r2.groupTuple()
            .map { it -> [it[0], it[1][0], it[2..-1].flatten()] }
            .set { mutect2_f1r2_set }
process mutct2_f1r2_comb {

  label 'med_mem'

  input:
  tuple val(caseID), val(sampleID), file(mutect2_ro) from mutect2_f1r2_set

  output:
  tuple val(caseID), val(sampleID), file("${sampleID}.mutect2.f1r2.tar.gz") into mutect2_f1r2_comb

  script:
  """
  ALL_F1R2_INPUT=\$(for x in *.mutect2.f1r2.tar.gz; do echo -n "-I \$x "; done)
  gatk LearnReadOrientationModel \$ALL_F1R2_INPUT -O ${sampleID}.mutect2.f1r2.tar.gz
  """
}

/* 2.75: MuTect2 Contamination
*/
mutect2_contam
  .join(mutect2_merge) //3
  .join(mutect2_stats) //3
  .join(mutect2_f1r2_comb) //3
  .groupTuple()
  .map { it -> [it[0], it[1][0], it[2][0], it[3][0], it[4][0], it[5][0], it[6][0], it[8][0], it[10][0], it[12][0]].flatten() }
  .set { mutect2_contam_merge }

process mutct2_contam_filter {

  label 'med_mem'

  publishDir path: "$params.outDir/cases/$caseID/mutect2", mode: "copy", overwrite: true
  publishDir path: "$params.outDir/output/mutect2", mode: "copy", pattern = "*.[mutect2.raw.vcf, mutect2.snv_indel.pass.vcf]"

  input:
  tuple val(caseID), val(sampleID), file(tumourbam), file(tumourbai), val(germlineID), file(germlinebam), file(germlinebai), file(mergevcf), file(statsvcf), file(readorient) from mutect2_contam_merge
  file(fasta) from reference.fa
  file(fai) from reference.fai
  file(dict) from reference.dict
  file(gps_files) from reference.seqlevel
  file(intlist) from reference.intlist

  output:
  tuple val(caseID), file("${sampleID}.mutect2.snv_indel.pass.vcf") into mutect2_veping
  tuple val(caseID), file("${sampleID}.mutect2.raw.vcf") into mutect2_rawVcf
  file('*') into completedmutect2call

  script:
  def taskmem = task.memory == null ? "" : "--java-options \"-Xmx" + javaTaskmem("${task.memory}") + "\""
  def gpsgz = params.seqlevel == "exome " || "panel" ? "${gps_files}/${params.exomeTag}/af-only-gnomad.${params.exomeTag}.hg*.noChr.vcf.gz" : "${gps_files}/af-only-gnomad.wgs.hg*.noChr.vcf.gz"
  """
  gatk ${taskmem} \
    GetPileupSummaries \
    -I ${tumourbam} \
    -V \$(echo ${gpsgz}) \
    -O ${sampleID}".getpileupsummaries.table" \
    -L ${intlist}

  gatk CalculateContamination \
    -I ${sampleID}".getpileupsummaries.table" \
    -O ${sampleID}".calculatecontamination.table"

  CONTAM=\$(tail -n+2 ${sampleID}.calculatecontamination.table | cut -f 2 | cut -d "." -f 1)
  if [[ \$CONTAM != 0 ]]; then
    touch ${sampleID}".CONTAMINATION.WARNING.txt"
  fi

  gatk IndexFeatureFile \
    --input ${mergevcf}

  gatk ${taskmem} \
    FilterMutectCalls \
    --reference ${fasta} \
    --contamination-table ${sampleID}".calculatecontamination.table" \
    --interval-padding 5 \
    --output ${sampleID}".mutect2.FilterMutectCalls.vcf" \
    --unique-alt-read-count 3 \
    --variant ${mergevcf} \
    --stats ${statsvcf} \
    --disable-sequence-dictionary-validation true \
    --ob-priors ${readorient} \
    -L ${intlist}

  perl ${workflow.projectDir}/bin/filter_Lancet_Mutect2_Manta-Strelka2_Format.pl \
    ID=${sampleID} \
    DP=14 \
    MD=2 \
    VCF=${sampleID}".mutect2.FilterMutectCalls.vcf"
  """
}

// 2.8: Manta output is a pre-req for Strelka2, so call both here
process mntstr {

  label 'high_mem'

  publishDir path: "$params.outDir/cases/$caseID/manta-strelka2", overwrite: 'true', mode: "copy"
  publishDir path: "${params.outDir}/output/manta-strelka2", mode: "copy", pattern: '*[.strelka2.snv_indel.raw.vcf, .strelka2.snv_indel.pass.vcf]'

  input:
  tuple val(caseID), val(sampleID), file(tumourbam), file(tumourbai), val(germlineID), file(germlinebam), file(germlinebai) from mantastrelka2ing
  file(fasta) from reference.fa
  file(fai) from reference.fai
  file(dict) from reference.dict
  file(bed_files) from reference.seqlevel

  output:
  tuple val(caseID), file("${sampleID}.strelka2.snv_indel.pass.vcf") into strelka2_veping
  tuple val(caseID), file("${sampleID}.strelka2.raw.vcf") into strelka2_rawVcf
  file('*.txt') into log_mantastrelka
  tuple val(caseID), val(germlineID) into vcfGraGermline

  script:
  def bedgz = params.seqlevel == "wgs" ? "${bed_files}/wgs.bed.gz" : "${bed_files}/${params.exomeTag}/${params.exomeTag}.bed.gz"
  def callRegions = params.seqlevel == "exome" || "panel" ? "--exome --callRegions ${bedgz}" : "--callRegions ${bedgz}"
  """
  {
    configManta.py ${callRegions} --referenceFasta=${fasta} --normalBam=${germlinebam} --tumourBam=${tumourbam} --runDir=manta

    manta/runWorkflow.py -m local -j ${task.cpus}

    configureStrelkaSomaticWorkflow.py ${callRegions} --referenceFasta=${fasta} --indelCandidates=manta/results/variants/candidateSmallIndels.vcf.gz --normalBam=${germlinebam} --tumorBam=${tumourbam} --runDir=strelka2

    strelka2/runWorkflow.py -m local -j ${task.cpus}

    ##merge into raw snv_indel
    gatk MergeVcfs -I strelka2/results/variants/somatic.snvs.vcf.gz -I strelka2/results/variants/somatic.indels.vcf.gz -O tmp.strelka2.snv_indel.vcf

    ${workflow.projectDir}/bin/manta_strelka2_rename_filter.sh  tmp.strelka2.snv_indel.vcf tmp2.strelka2.snv_indel.vcf ${sampleID} ${germlineID}

    perl ${workflow.projectDir}/bin/filter_Lancet_Mutect2_Manta-Strelka2_Format.pl \
        ID=${sampleID} \
        DP=14 \
        MD=2 \
        VCF=tmp2.strelka2.snv_indel.vcf

  } 2>&1 | tee > ${sampleID}.manta-strelka2.log.txt
  """
}

/* 2.91: Lancet
*/
lancetbedding = lancet_bedding.flatten()
lanceting
  .map { it -> [it[0],it[1],it[2],it[3],it[4],it[5],it[6]]}
  .combine(lancetbedding)
  .set { lancetsbedding }

process lancet_sg {

  label 'med_mem'

  input:
  tuple val(caseID), val(sampleID), file(tumourbam), file(tumourbai), val(germlineID), file(germlinebam), file(germlinebai), file(bed) from lancetsbedding
  file(fasta) from reference.fa
  file(fai) from reference.fai
  file(dict) from reference.dict

  output:
  tuple val(caseID), val(sampleID), file('*.sort.lancet.vcf') into lancet_gt

  when:
  params.seqlevel == "exome" || "panel"

  script:
  scatgathn = "${bed}".split("\\.")[2]
  """
  lancet \
    --num-threads ${task.cpus} \
    --ref ${fasta} \
    --bed ${bed} \
    --tumor ${tumourbam} \
    --normal ${germlinebam} | \
    perl -ane 'if(\$F[0]=~m/^\\#CHROM/){
      \$_=~s/TUMOR/${sampleID}/;
      \$_=~s/NORMAL/${germlineID}/;
      print \$_;}
    else{print \$_;}' > ${sampleID}"."${scatgathn}".lancet.vcf"

  picard SortVcf \
    I=${sampleID}"."${scatgathn}".lancet.vcf" \
    O=${sampleID}"."${scatgathn}".sort.lancet.vcf" \
    SD=${dict}
  """
}

/* 2.92: Lancet Merge
*/
lancet_gt
  .groupTuple()
  .map { it -> tuple(it[0], it[1][0], it[2][0..-1].flatten()) }
  .set { lancet_fm }

process lancet_concat {

  label 'med_mem'

  input:
  tuple val(caseID), val(sampleID), file(rawvcf) from lancet_fm

  output:
  tuple val(caseID), val(sampleID), file('*lancet.merge.vcf') into lancet_merge

  script:
  """
  ls *.sort.lancet.vcf > vcf.list
  picard MergeVcfs I=vcf.list O=${sampleID}".lancet.merge.vcf"
  """
}

/* 2.93: Lancet Filter
*/
process lancet_filter {

  label 'med_mem'

  publishDir path: "$params.outDir/cases/$caseID/lancet"
  publishDir path: "${params.outDir}/output/lancet", mode: "copy", pattern: '*.[raw.vcf, .lancet.snv_indel.pass.vcf]'

  input:
  tuple val(caseID), val(sampleID), file(mergevcf) from lancet_merge

  output:
  tuple val(caseID), file("${sampleID}.lancet.snv_indel.pass.vcf") into lancet_veping
  tuple val(caseID), file("${sampleID}.lancet.raw.vcf") into lancet_rawVcf
  file('*') into completedlancetcall

  script:
  """
  perl ${workflow.projectDir}/bin/filter_Lancet_Mutect2_Manta-Strelka2_Format.pl \
    ID=${sampleID} \
    DP=14 \
    MD=2 \
    VCF=${mergevcf}
  """
}

/*
================================================================================
                          3.  ANNOTATION AND REPORTING
================================================================================
*/
// 3.0: Annotate Vcfs

lancet_veping
  .join( mutect2_veping )
  .join( strelka2_veping )
  .groupTuple()
  .map { it -> [it[0], it[1][0], it[2][0], it[3][0]].flatten() }
  .set { case_veping }

process vepann {

  label 'med_mem'

  publishDir path: "${params.outDir}/cases/${caseID}/vcf", mode: "copy", pattern: '*.vcf'
  publishDir path: "${params.outDir}/output/VEP", mode: "copy", pattern: '*.vcf'

  input:
  tuple val(caseID), file(vcf1), file(vcf2), file(vcf3) from case_veping
  file(fasta) from reference.fa
  file(fai) from reference.fai
  file(dict) from reference.dict
  file(grchver) from reference.grchvers
  file(pcgrbase) from reference.pcgrbase

  output:
  tuple val(caseID), file('*vep.vcf') into runGRanges

  script:
  def sampleID = "${vcf1}".split("\\.")[0]
  def grch_vers = "${grchver}".split("\\/")[-1]
  """
  for VCF in $vcf1 $vcf2 $vcf3; do
    VCFANNO=\$(echo \$VCF | sed "s/.vcf/.vep.vcf/")
    vep --dir_cache ${pcgrbase}/data/${grch_vers}/.vep \
      --offline \
      --assembly ${params.assembly} \
      --vcf_info_field ANN \
      --symbol \
      --species homo_sapiens \
      --check_existing \
      --cache \
      --fork ${task.cpus} \
      --af_1kg \
      --af_gnomad \
      --vcf \
      --input_file \$VCF \
      --output_file \$VCFANNO \
      --format "vcf" \
      --fasta ${fasta} \
      --hgvs \
      --canonical \
      --ccds \
      --force_overwrite \
      --verbose
  done
  """
}

/* 3.1 RData GRanges from processed VCFs
* take publishDir and check for number of files therein
* each sample has 6 associated (raw, pass per caller)
* NB increment if adding callers!
*/
runGRanges
  .join(mutect2_rawVcf)
  .join(strelka2_rawVcf)
  .join(lancet_rawVcf)
  .groupTuple()
  .map { it -> [it[0], it[1][0], it[1][1], it[1][2], it[2], it[3], it[4]] }
  .join(vcfGraGermline)
  .into { cons_vcfs; cons_comb  }

Channel
  .from("${params.impacts}")
  .splitCsv(sep : ';')
  .flatten()
  .set {impactss}

process vcfGRa {

  label 'med_mem'

  publishDir "$params.outDir/cases/$caseID/consensus_vcfs"
  publishDir "${params.outDir}/output/consensus_variants/pdf", mode: "copy", pattern: '*.pdf'
  publishDir "${params.outDir}/output/consensus_variants/vcf", mode: "copy", pattern: '*.impacts.pcgr.tsv.vcf'
  publishDir "${params.outDir}/output/consensus_variants/data", mode: "copy", pattern: '*[.RData, .tsv]'

  input:
  tuple val(caseID), file(vvcf1), file(vvcf2), file(vvcf3), file(rvcf1), file(rvcf2), file(rvcf3), val(germlineID) from cons_vcfs
  each impact from impactss

  output:
  tuple val(caseID), val(sampleID), file("${sampleID}.*impacts.pcgr.tsv.vcf") into vcfs_pcgr
  file('*') into completedvcfGRangesConsensus

  script:
  def sampleID = "${vvcf1}".split("\\.")[0]
  def inc_ord = params.incOrder ? params.incOrder : "noord"
  def which_genome = params.assembly == "GRCh37" ? "hg19" : "hg38"
  """
  echo ${sampleID} > sId.txt && rm sID.txt
  Rscript -e "somenone::variant_consensus(germline_id = \\"${germlineID}\\", vep_vcf_pattern = \\"snv_indel.pass.vep.vcf\\", raw_vcf_pattern = \\"raw.vcf\\", tag = \\"${params.runID}\\", which_genome = \\"${which_genome}\\", included_order = \\"${inc_ord}\\", impacts = \\"${impact}\\")"
  """
}

// 3.2 Create VCF for PCGR from consensus
process pcgrVcf {

  label 'low_mem'

  input:
  tuple val(caseID), val(sampleID), file(cons_vcf) from vcfs_pcgr

  output:
  tuple val(sampleID), file("${pcgr_vcf}") into snvpass_pcgr

  when:
  cons_vcf =~ "${sampleID}.${params.runID}.HMML_impacts.pcgr.tsv.vcf"

  script:
  pcgr_vcf = "${cons_vcf}".replace("pcgr.tsv.vcf", "snv_indel.pass.pcgr.vcf")
  """
  cat ${workflow.projectDir}/assets/vcf42.head.txt > $pcgr_vcf
  head -n1 $cons_vcf >> $pcgr_vcf
  tail -n+2 $cons_vcf | sort -V >> $pcgr_vcf
  """
}

if(params.facets){
  snvpass_pcgr
    .join(facets_pcgr)
    .map { it -> [it[0], it[1], it[2], it[3]] }
    .set { pcgr_inputs }
}
if(! params.facets){
  snvpass_pcgr
    .map { it -> [it[0], it[1], null, null] }
    .set { pcgr_inputs }
}

/* 3.3 PCGR report
* take all mutations in consensus.tab from pass.vcfs into single VCF for PCGR
*/
process pcgrreport {

  label 'low_mem'

  publishDir "${params.outDir}/reports/pcgr", mode: "copy", pattern: "*html"
  publishDir "${params.outDir}/cases/${caseID}/pcgr", mode: "copy"

  input:
  tuple val(sampleID), file(vcf), file(jointsegs), file(ploidpur) from pcgr_inputs
  file(grchver) from reference.grchvers
  file(pcgrbase) from reference.pcgrbase

  output:
  file('*') into completedPCGR
  tuple val(sampleID), file(vcf), file("${sampleID}.pcgr_acmg.grch38.pass.vcf.gz") into cons_comb_pcgr

  script:
  caseID="${sampleID}".split("${params.tumourIDsplit}")[0]
  grch_vers = "${grchver}".split("\\/")[-1]
  ploidy = ploidpur =~ "input." ? "" : "--tumor_ploidy \$(cut -f 1 ${ploidpur})"
  purity = ploidpur =~ "input." ? "" : "--tumor_purity \$(cut -f 2 ${ploidpur})"
  jointseg = jointsegs =~ "input." ? "" : "--input_cna ${jointsegs}"
  config = params.seqlevel == "exome" || "panel" ? "${pcgrbase}/data/${grch_vers}/pcgr_configuration_${params.exomeTag}.toml" : "${pcgrbase}/data/${grch_vers}/pcgr_configuration_wgs.toml"
  """
  {
    ##PCGR 0.9.0rc
    pcgr.py \
      --pcgr_dir ${pcgrbase} \
      --output_dir ./ \
      --genome_assembly ${grch_vers} \
      --conf ${config} \
      --sample_id ${sampleID} \
      --input_vcf ${vcf} ${jointseg} ${ploidy} ${purity} \
      --no-docker \
      --force_overwrite \
      --no_vcf_validate

  } 2>&1 | tee > ${sampleID}.pcgr.log.txt
  """
}

//  3.4 PCGR + consensus VCF to make a good VCF for combining
process prepvepSomtsv {

  label 'low_mem'

  input:
  tuple val(sampleID), file(vcf), file(vep) from cons_comb_pcgr

  output:
  file("${sampleID}.${params.runID}.combine.vcf") into vep_som_tsv

  script:
  """
  perl ${workflow.projectDir}/bin/combine_select_elements.pl \
    $vep \
    $vcf \
    ${sampleID}.${params.runID}.combine.vcf
  """
}

process vepSomtsv {

  label 'low_mem'

  publishDir path: "${params.outDir}/output/consensus_variants/combine_vcf", mode: "copy"

  input:
  file(vcf) from vep_som_tsv.collect()

  output:
  file("${params.runID}.consensus_variants.tabvcf.tsv") into madetsv2

  script:
  """
  perl ${workflow.projectDir}/bin/vepHCvcf_combine_tsv.pl "${params.runID}.consensus_variants"
  """
}

/*
================================================================================
                          4.  MULTIQC AND CLOSEOUT
================================================================================
*/
// 4.0 Run multiQC to finalise report
if(! params.agilentUMI){
  dedup_multiqc_all = mrkdup_multiqc.collect()
}
if(params.agilentUMI) {
  dedup_multiqc_all = locatit_flagsqc.collect()
}

process mltiQC {

  label 'low_mem'
  publishDir path: "$params.outDir/reports", mode: "copy", pattern: "*html"

  input:
  file(fastps) from fastp_multiqc.collect()
  file(gtkrcls) from gtkrcl_multiqc.collect()
  file(multimetrics) from multimetrics_multiqc.collect()
  file(mrkdups) from dedup_multiqc_all

  output:
  file('*') into completedmultiqc

  script:
  """
  multiqc . -i ${params.runID}".somatic_n-of-1" --tag DNA -f -c ${params.multiqcConfig}
  """
}

// 4.1.1: somatic_n-of-1 container software versions
process somenone_software_vers {

  label 'low_mem'
  publishDir "${params.outDir}/../pipeline_info", mode: 'copy'

  output:
  file 'somenone_software_versions.yaml' into ch_somenone_software_vers

  script:
  """
  conda env export > somenone_software_versions.yaml
  """
}
