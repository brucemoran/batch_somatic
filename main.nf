#!/usr/bin/env nextflow

params.help = ""

if (params.help) {
  log.info ''
  log.info '-----------------------------------------------------------------------'
  log.info 'NEXTFLOW 19.10 FASTQ QC, TRIM, ALIGN, SOMATIC+GERMLINE SNV, CNA, REPORT BATCH IMPLEMENTATION'
  log.info '-----------------------------------------------------------------------'
  log.info ''
  log.info 'Usage: '
  log.info 'nextflow run brucemoran/batch_somatic'
  log.info ''
  log.info 'Mandatory arguments:'
  log.info '    -profile    Configuration profile (required: standard,singularity)'
  log.info '    --sampleCsv      STRING      CSV format, headers: caseID, soma_sampleID, soma_read1, soma_read2, germ_sampleID, germ_read1, germ_read2'
  log.info '    --refDir        STRING      dir in which reference data and required indices are held; if not specified, this is created by brucemoran/DNAseq_references/GRCh_37-38 (default: work/GRCh38)'
  log.info ''
  log.info 'Optional arguments:'
  log.info '    --exometag        STRING      if multiple exomes have reference data, specify a tag for matching (should have been used at reference generation )'
  log.info '    --germline      STRING      run HaplotypeCaller on germline sample and annotate with CPSR (true/false, default: true)'
  log.info '    --scatGath      NUM         number of pieces to divide intervalList into for scattering to variant calling processes (default: 20 for exome, 100 for WGS)'
  log.info '    --multiqcConfig      STRING      config file for multiqc (default: bin/somatic_n-of-1.multiQC_config.yaml)'
  log.info '    --seqLevel      STRING      WGS or exome (default: WGS)'
  log.info ''
  exit 1
}

/* -2 Test if refDir is defined, if not run DNAseq_references pipeline under defaults
*/
if(!params.refDir){
  exit 1, "Please run: nextflow run brucemoran/DNAseq_references --outDir work -profile standard,singularity, then specify: nextflow run brucemoran/batch_somatic --refDir work/GRCh38"
}

/* -1: Global Variables
*/
params.outDir = "analysis/batch/${params.seqLevel}"

// Reference data as params, and reusable therefore
params.fasta = Channel.fromPath("$params.refDir/*fasta").getVal()
params.fai = Channel.fromPath("$params.refDir/*fasta.fai").getVal()
params.dict = Channel.fromPath("$params.refDir/*dict").getVal()

params.amb = Channel.fromPath("$params.refDir/*fasta.amb").getVal()
params.ann = Channel.fromPath("$params.refDir/*fasta.ann").getVal()
params.bwt = Channel.fromPath("$params.refDir/*fasta.bwt").getVal()
params.pac = Channel.fromPath("$params.refDir/*fasta.pac").getVal()
params.sa = Channel.fromPath("$params.refDir/*fasta.sa").getVal()

params.twobit = Channel.fromPath("$params.refDir/*fasta.2bit").getVal()

params.seqlevel = "$params.seqLevel".toLowerCase()
if(params.seqlevel == "wgs"){
  params.exometag = "wgs"
}
params.intlist = Channel.fromPath("$params.refDir/${params.seqlevel}/${params.exometag}*bed.interval_list").getVal()
params.bed = Channel.fromPath("$params.refDir/${params.seqlevel}/${params.exometag}*bed").getVal()
params.bedgz = Channel.fromPath("$params.refDir/${params.seqlevel}/${params.exometag}*bed.gz").getVal()
params.bedgztbi = Channel.fromPath("$params.refDir/${params.seqlevel}/${params.exometag}*bed.gz.tbi").getVal()

params.dbsnp = Channel.fromPath("$params.refDir/dbsnp*.gz").getVal()
params.dbsnptbi = Channel.fromPath("$params.refDir/dbsnp*.tbi").getVal()

params.omni = Channel.fromPath("$params.refDir/KG_omni*.gz").getVal()
params.otbi = Channel.fromPath("$params.refDir/KG_omni*.gz.tbi").getVal()
params.kgp1 = Channel.fromPath("$params.refDir/KG_phase1*.gz").getVal()
params.ktbi = Channel.fromPath("$params.refDir/KG_phase1*.gz.tbi").getVal()
params.hpmp = Channel.fromPath("$params.refDir/hapmap*.gz").getVal()
params.htbi = Channel.fromPath("$params.refDir/hapmap*.gz.tbi").getVal()

params.gps = Channel.fromPath("$params.refDir/${params.seqlevel}/af-only-gnomad*${params.exometag}*vcf.gz").getVal()
params.gpstbi = Channel.fromPath("$params.refDir/${params.seqlevel}/af-only-gnomad*${params.exometag}*vcf.gz.tbi").getVal()

//PCGR, CPSR version and base data dir
params.grchvers = Channel.fromPath("$params.refDir/pcgr/data").getVal()
params.pcgrdir = Channel.fromPath("$params.refDir/pcgr").getVal()
params.pcgrtoml = Channel.fromPath("$params.refDir/pcgr/data/${params.grchvers}/pcgr_configuration_${params.exometag}.toml").getVal()

//somaticVariantConsensus scripts
process cons_scripts {

  publishDir "${workflow.projectDir}/bin"

  output:
  tuple file('variants_GRanges_consensus_plot_batch.call.R'), file('variants_GRanges_consensus_plot_batch.func.R') into vcfGRa_Scripts

  script:
  """
  wget "https://raw.githubusercontent.com/brucemoran/somaticVariantConsensus/master/scripts/variants_GRanges_consensus_plot_batch.func.R"
  wget "https://raw.githubusercontent.com/brucemoran/somaticVariantConsensus/master/scripts/variants_GRanges_consensus_plot_batch.call.R"
  """
}

//Java task memory allocation via task.memory
javaTaskmem = { it.replace(" GB", "g") }

/* 0.00: Input using sample.csv
*/
Channel.fromPath("$params.sampleCsv", type: 'file')
       .splitCsv( header: true )
       .map { row -> [row.caseID, row.soma_sampleID, file(row.soma_read1), file(row.soma_read2), row.germ_sampleID, file(row.germ_read1), file(row.germ_read2)] }
       .set { split_soma_germ }

/* 0.00: Input using sample.csv
*/
process splt_sg {
  executor 'local'
  publishDir "$params.outDir/cases/$caseID", mode: "copy", pattern: "*.csv"

  input:
  set val(caseID), val(soma_sampleID), file(soma_read1), file(soma_read2), val(germ_sampleID), file(germ_read1), file(germ_read2) from split_soma_germ

  output:
  file('*.csv') into splitcsv2

  """
  SR1=\$(readlink -e $soma_read1)
  SR2=\$(readlink -e $soma_read2)
  GR1=\$(readlink -e $germ_read1)
  GR2=\$(readlink -e $germ_read2)
  echo "caseID,type,sampleID,read1,read2" > $caseID".csv"
  echo "$caseID,somatic,$soma_sampleID,\$SR1,\$SR2" >> $caseID".csv"
  echo "$caseID,germline,$germ_sampleID,\$GR1,\$GR2" >> $caseID".csv"
  """
}

/*two set channels need to be processed the same; use
*/
splitcsv2.splitCsv( header: true )
         .map { row -> [row.caseID, row.type, row.sampleID, file(row.read1), file(row.read2)] }
         .set { bbduking }

/* 0.01: Input trimming
*/
process bbduk {

  label 'med_mem'

  publishDir path: "$params.outDir/cases/$caseID/bbduk", mode: "copy", pattern: "*.txt"

  input:
  set val(caseID), val(type), val(sampleID), file(read1), file(read2) from bbduking

  output:
  file('*.txt') into log_bbduk
  tuple val(caseID), val(type), val(sampleID), file('*.bbduk.R1.fastq.gz'), file('*.bbduk.R2.fastq.gz') into bwa_memming
  tuple val(caseID), val(type), val(sampleID), file(read1), file(read2) into fastping

  script:
  taskmem = javaTaskmem("${task.memory}")
  """
  {
  ##remove because of issues on readname including "1:xyz, 2:xyz for R1, R2"
  if [[ $read1 =~ ".gz"\$ ]]; then
    gunzip -c $read1 | perl -ane 'chomp; print "\$F[0]\\n";' > r1.fq
    gunzip -c $read2 | perl -ane 'chomp; print "\$F[0]\\n";' > r2.fq
  else
    perl -ane 'chomp; print "\$F[0]\\n";' $read1 > r1.fq
    perl -ane 'chomp; print "\$F[0]\\n";' $read2 > r2.fq
  fi

  ##repair in case of disorder
  repair.sh in1=r1.fq in2=r2.fq out1=r1.f.fq out2=r2.f.fq repair

  sh bbduk.sh -Xmx$taskmem \
    in1=r1.f.fq \
    in2=r2.f.fq \
    out1=$sampleID".bbduk.R1.fastq.gz" \
    out2=$sampleID".bbduk.R2.fastq.gz" \
    k=31 \
    mink=5 \
    hdist=1 \
    ktrim=r \
    trimq=20 \
    qtrim=rl \
    maq=20 \
    ref=/opt/miniconda/envs/somatic_n-of-1/opt/bbmap-38.57-0/resources/adapters.fa \
    tpe \
    tbo \
    stats=$sampleID".bbduk.adapterstats.txt" \
    overwrite=T
  } 2>&1 | tee > ${sampleID}.bbduk.runstats.txt
  rm r*.fq
  """
}

/* 0.2: fastp QC of pre-, post-bbduk
*/
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
  fastp -w ${task.cpus} -h $sampleID".fastp.html" -j $sampleID".fastp.json" --in1 $preread1 --in2 $preread2
  """
}

/* 1.0: Input alignment
*/
process bwamem {

  label 'high_mem'

  input:
  tuple val(caseID), val(type), val(sampleID), file(read1), file(read2) from bwa_memming
  tuple file(fasta), file(fai), file(dict) from Channel.value([params.fasta, params.fai, params.dict])
  tuple file(amb), file(ann), file(bwt), file(pac), file(sa) from Channel.value([params.amb, params.ann, params.bwt, params.pac, params.sa])

  output:
  tuple val(caseID), val(type), val(sampleID), file('*.bam'), file('*.bai') into (cramming, dup_marking)

  script:
  """
  DATE=\$(date +"%Y-%m-%dT%T")
  RGLINE="@RG\\tID:$sampleID\\tPL:ILLUMINA\\tSM:$sampleID\\tDS:$type\\tCN:UCD\\tLB:LANE_X\\tDT:\$DATE"

  bwa mem \
    -t${task.cpus} \
    -M \
    -R \$RGLINE \
    $fasta \
    $read1 $read2 | \
    samtools sort -T "tmp."$sampleID -o $sampleID".sort.bam"
  samtools index $sampleID".sort.bam"
  """
}

/* 1.1: CRAM alignment and output
* TODO: output upload schema for ENA/EGA
*/
process cram {

  label 'low_mem'

  publishDir path: "$params.outDir/cases/$caseID/bwa", mode: "copy", pattern: "*.cra*"

  input:
  tuple val(caseID), val(type), val(sampleID), file(bam), file(bai) from cramming
  tuple file(fasta), file(fai), file(dict) from Channel.value([params.fasta, params.fai, params.dict])

  output:
  tuple file("${sampleID}.sort.cram"), file("${sampleID}.sort.cram.crai") into completedcram

  script:
  """
  samtools view -hC -T $fasta $sampleID".sort.bam" > $sampleID".sort.cram"
  samtools index $sampleID".sort.cram"
  """
}

/* 1.2: MarkDuplicates
*/
process mrkdup {

  label 'high_mem'

  publishDir path: "$params.outDir/cases/$caseID/picard", mode: "copy", pattern: "*.txt"

  input:
  tuple val(caseID), val(type), val(sampleID), file(bam), file(bai) from dup_marking

  output:
  file('*.txt') into mrkdup_multiqc
  tuple val(caseID), val(type), val(sampleID), file('*.md.bam'), file('*.md.bam.bai') into gatk4recaling

  script:
  taskmem = javaTaskmem("${task.memory}")
  """
  OUTBAM=\$(echo $bam | sed 's/bam/md.bam/')
  OUTMET=\$(echo $bam | sed 's/bam/md.metrics.txt/')
  {
  picard -Xmx$taskmem \
    MarkDuplicates \
    TMP_DIR=./ \
    INPUT=$bam \
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

/* 1.3: GATK4 BestPractices
*/
process gtkrcl {

  label 'high_mem'

  publishDir path: "$params.outDir/cases/$caseID/gatk4/bqsr", mode: "copy", pattern: "*.GATK4_BQSR.log.txt"

  input:
  tuple val(caseID), val(type), val(sampleID), file(bam), file(bai) from gatk4recaling
  tuple file(fasta), file(fai), file(dict) from Channel.value([params.fasta, params.fai, params.dict])
  tuple file(dbsnp), file(dbsnptbi) from Channel.value([params.dbsnp, params.dbsnptbi])
  file(intlist) from Channel.value(params.intlist)

  output:
  file('*.table') into gtkrcl_multiqc
  tuple val(caseID), val(type), val(sampleID), file('*.bqsr.bam'), file('*.bqsr.bam.bai') into ( germfiltering, gmultimetricing)
  tuple val(caseID), val(type), val(sampleID), file('*.bqsr.bam'), file('*.bqsr.bam.bai') into hc_germ
  tuple val(caseID), val(sampleID) into metas_pcgr
  file("${sampleID}.GATK4_BQSR.log.txt") into bqsr_log

  script:
  """
  {
  gatk BaseRecalibrator \
    -R $fasta \
    -I $bam \
    --known-sites $dbsnp \
    --use-original-qualities \
    -O ${sampleID}.recal_data.table \
    --disable-sequence-dictionary-validation true \
    -L $intlist

  #ApplyBQSR
  OUTBAM=\$(echo $bam | sed 's/bam/bqsr.bam/')
  gatk ApplyBQSR \
    -R $fasta \
    -I $bam \
    --bqsr-recal-file ${sampleID}.recal_data.table \
    --add-output-sam-program-record \
    --use-original-qualities \
    -O \$OUTBAM \
    -L $intlist

  samtools index \$OUTBAM \$OUTBAM".bai"
  } 2>&1 | tee >  ${sampleID}.GATK4_BQSR.log.txt
  """
}

/* 1.31: scatter-gather implementation for mutect2, lancet
*/
process scat_gath {

  label 'low_mem'

  input:
  file(intlist) from Channel.value(params.intlist)

  output:
  file('lancet.scatgath.*.bed') into lancet_bedding
  file('mutect2.scatgath.*.bed.interval_list') into mutect2_bedding
  file('hc.scatgath.*.bed.interval_list') into hc_bedding

  script:
  def sgcount = params.scatGath
  if (params.scatGath == null){
    if (params.seqlevel == "exome"){
      sgcount = 20
    }
    else {
      sgcount = 100
    }
  }
  """
  picard IntervalListTools \
    I=$intlist \
    SCATTER_COUNT=$sgcount \
    O=./
  ls temp*/* | while read FILE; do
    COUNTN=\$(dirname \$FILE | perl -ane '@s=split(/\\_/); print \$s[1];');
    mv \$FILE mutect2.scatgath.\${COUNTN}.bed.interval_list;
    cp mutect2.scatgath.\${COUNTN}.bed.interval_list hc.scatgath.\${COUNTN}.bed.interval_list
    grep -v @ mutect2.scatgath.\${COUNTN}.bed.interval_list | \
      cut -f 1,2,3,5 > lancet.scatgath.\${COUNTN}.bed
  done
  """
}

/* 1.4: GATK4 Germline Haplotypecaller
*/
hcbedding = hc_bedding.flatten()
hc_germ
  .map { it -> [it[0],it[1],it[2],it[3],it[4]] }
  .combine(hcbedding)
  .set { hcgermbedding }

process haplotypecaller {

  label 'med_mem'

  input:
  tuple val(caseID), val(type), val(sampleID), file(bam), file(bai), file(intlist) from hcgermbedding
  tuple file(fasta), file(fai), file(dict) from Channel.value([params.fasta, params.fai, params.dict])
  tuple file(dbsnp), file(dbsnptbi) from Channel.value([params.dbsnp, params.dbsnptbi])
  tuple file(omni), file(otbi), file(kgp1), file(ktbi), file(hpmp), file(htbi) from Channel.value([params.omni, params.otbi, params.kgp1, params.ktbi, params.hpmp, params.htbi])

  output:
  tuple val(caseID), val(sampleID), file('*sort.hc.vcf') into hc_gt

  when:
  type == "germline" & params.germline != false

  script:
  taskmem = javaTaskmem("${task.memory}")
  """
  SCATGATHN=\$(echo $intlist | perl -ane '@s=split(/\\./);print \$s[2];')
  gatk --java-options -Xmx$taskmem HaplotypeCaller \
    -R $fasta \
    -I $bam \
    --dont-use-soft-clipped-bases \
    --standard-min-confidence-threshold-for-calling 20 \
    --dbsnp $dbsnp \
    --native-pair-hmm-threads ${task.cpus} \
    -O $sampleID".\${SCATGATHN}.hc.vcf" \
    --disable-sequence-dictionary-validation true \
    -L $intlist

  picard SortVcf \
    I=$sampleID".\${SCATGATHN}.hc.vcf" \
    O=$sampleID".\${SCATGATHN}.sort.hc.vcf" \
    SD=$dict
  """
}

// 2.42: HC_merge
hc_gt
  .groupTuple()
  .map { it -> tuple(it[0], it[1][0], it[2][0..-1].flatten()) }
  .set { hc_fm }

process hc_merge {

  label 'low_mem'

  publishDir path: "$params.outDir/cases/$caseID/gatk4/haplotypecaller", mode: "copy", pattern: '*.hc.merge.vcf'

  input:
  tuple val(caseID), val(sampleID), file(rawvcfs) from hc_fm

  output:
  tuple val(caseID), val(sampleID), file("${sampleID}.hc.merge.vcf") into cpsr_vcf

  script:
  """
  ls *.sort.hc.vcf > vcf.list
  picard MergeVcfs I=vcf.list O=$sampleID".hc.merge.vcf"
  """
}

/* 1.25: CPSR annotation of GATK4 Germline
*/
process cpsrreport {

  label 'low_mem'

  publishDir "$params.outDir/reports/cpsr", mode: "copy", pattern: "*.html"
  publishDir "$params.outDir/cases/$caseID/cpsr", mode: "copy"

  input:
  tuple val(caseID), val(sampleID), file(vcf) from cpsr_vcf

  output:
  file('*') into cpsr_vcfs

  script:
  """
  {
  VERS=\$(ls ${params.pcgrdir}/data)
  CONFIG=\$(readlink -e ${params.pcgrdir}/data/*/cpsr_configuration_default.toml)

  ##CPSR v0.5.2.2
  cpsr.py \
    --no-docker \
    --no_vcf_validate \
    --panel_id 0 \
    $vcf \
    ${params.pcgrdir} \
    ./ \
    \$VERS \
    \$CONFIG \
    $sampleID

  } 2>&1 | tee > ${sampleID}.cpsr.log.txt
  """
}

/* 2.0: filter germline channel, tap into somatic channels for all processes subsequent
*/
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

/* 2.1: PicardTools metrics suite for MultiQC HTML report
*/
process mltmet {

  label 'med_mem'

  publishDir "$params.outDir/cases/$caseID/metrics", overwrite: 'true'

  input:
  tuple val(caseID), val(type), val(sampleID), file(bam), file(bai) from gmultimetricing
  tuple file(fasta), file(fai), file(dict) from Channel.value([params.fasta, params.fai, params.dict])
  file(intlist) from Channel.value(params.intlist)

  output:
  file('*.txt') into multimetrics_multiqc

  script:
  taskmem = javaTaskmem("${task.memory}")
  """
  {
  if [[ ${params.seqlevel} == "exome" ]]; then
  picard -Xmx$taskmem CollectHsMetrics \
    I=$bam \
    O=$sampleID".hs_metrics.txt" \
    TMP_DIR=./ \
    R=$fasta \
    BAIT_INTERVALS=$intlist  \
    TARGET_INTERVALS=$intlist
  fi
  picard -Xmx$taskmem CollectAlignmentSummaryMetrics \
    I=$bam \
    O=$sampleID".AlignmentSummaryMetrics.txt" \
    TMP_DIR=./ \
    R=$fasta

  picard -Xmx$taskmem CollectMultipleMetrics \
    I=$bam \
    O=$sampleID".CollectMultipleMetrics.txt" \
    TMP_DIR=./ \
    R=$fasta

  picard -Xmx$taskmem CollectSequencingArtifactMetrics \
    I=$bam \
    O=$sampleID".artifact_metrics.txt" \
    TMP_DIR=./ \
    R=$fasta

  picard -Xmx$taskmem CollectInsertSizeMetrics \
    I=$bam \
    O=$sampleID".insert_size_metrics.txt" \
    H=$bam".histogram.pdf" \
    TMP_DIR=./

  } 2>&1 | tee > ${sampleID}.picard.metrics.log
  """
}

/*2.21: SCNA with facets CSV snp-pileup
*/

process fctcsv {

  label 'med_mem'

  publishDir "$params.outDir/cases/$caseID/facets"

  input:
  tuple val(caseID), val(sampleID), file(tumourbam), file(tumourbai), val(germlineID), file(germlinebam), file(germlinebai) from facetsomaing
  tuple file(dbsnp), file(dbsnptbi) from Channel.value([params.dbsnp, params.dbsnptbi])

  output:
  file("${sampleID}.cncf-jointsegs.pcgr.tsv") into facets_consensusing
  tuple val(sampleID), file("${sampleID}.cncf-jointsegs.pcgr.tsv"), file("${sampleID}.fit_ploidy-purity.pcgr.tsv") into facets_pcgr
  file("${sampleID}.facets.log.txt") into log_facets

  when:
  params.facets != false

  script:
  """
  {
  snp-pileup \
    $dbsnp \
    -r 10 \
    -p \
    ${sampleID}.facets.r10.csv \
    $germlinebam \
    $tumourbam

  WCTEST=\$(wc -l ${sampleID}.facets.r10.csv | perl -ane 'print \$F[0];')

  if [[ \$WCTEST > 1 ]]; then
    Rscript --vanilla ${workflow.projectDir}/bin/facets_cna.call.R ${sampleID}.facets.r10.csv

    echo -e "Chromosome\\tStart\\tEnd\\tSegment_Mean" > $sampleID".cncf-jointsegs.pcgr.tsv"
    tail -n+2 $sampleID".fit_cncf-jointsegs.tsv" | awk '{print \$1"\\t"\$10"\\t"\$11"\\t"\$5}' >> $sampleID".cncf-jointsegs.pcgr.tsv"

    tail -n+2 $sampleID".fit_ploidy-purity.tab" > $sampleID".fit_ploidy-purity.pcgr.tsv"
  else
    echo -e "Chromosome\\tStart\\tEnd\\tSegment_Mean" > $sampleID".cncf-jointsegs.pcgr.tsv"
    echo -e "NA\tNA" > $sampleID".fit_ploidy-purity.pcgr.tsv"
  fi
  } 2>&1 | tee > ${sampleID}.facets.log.txt
  """
}

/* 2.22: SCNA consensus from facets
*/
process fctcon {

  label 'med_mem'

  publishDir "$params.outDir/reports/scna/facets"

  input:
  file(filesn) from facets_consensusing.collect()
  file(dict) from Channel.value(params.dict)

  output:
  file('*') into complete_facets

  script:
  """
  {
  Rscript --vanilla ${workflow.projectDir}/bin/facets_cna_consensus.call.R \
    $dict \
    ./ \
    ${workflow.projectDir}/bin/facets_cna_consensus.func.R
  } 2>&1 | tee > facets_cons.log.txt
  """
}

mutect2bedding = mutect2_bedding.flatten()
mutect2somaticing
  .map { it -> [it[0],it[1],it[2],it[3],it[4],it[5],it[6]]}
  .combine(mutect2bedding)
  .set { mutect2somaticbedding }

/* 2.41: MuTect2
* NB --germline-resource dollar-sign{dbsnp} removed as no AF causing error
*/
process mutct2_sg {

  label 'med_mem'

  input:
  tuple val(caseID), val(sampleID), file(tumourbam), file(tumourbai), val(germlineID), file(germlinebam), file(germlinebai), file(intlist) from mutect2somaticbedding
  tuple file(fasta), file(fai), file(dict) from Channel.value([params.fasta, params.fai, params.dict])

  output:
  tuple val(caseID), val(sampleID), file('*sort.mutect2.vcf') into mutect2_gt
  tuple val(caseID), val(sampleID), file('*.vcf.stats') into mutect2_st
  tuple val(caseID), val(sampleID), file('*mutect2.f1r2.tar.gz') into mutect2_f1r2

  script:
  taskmem = javaTaskmem("${task.memory}")
  """
  SCATGATHN=\$(echo $intlist | perl -ane '@s=split(/\\./);print\$s[2];')
  gatk --java-options -Xmx$taskmem \
    Mutect2 \
    --native-pair-hmm-threads ${task.cpus} \
    --reference $fasta \
    --input $germlinebam \
    --input $tumourbam \
    --normal-sample $germlineID \
    --tumor-sample $sampleID \
    --output $sampleID"."\${SCATGATHN}".mutect2.vcf" \
    --disable-sequence-dictionary-validation true \
    --f1r2-tar-gz \${SCATGATHN}".mutect2.f1r2.tar.gz" \
    -L $intlist

  picard SortVcf \
    I=$sampleID"."\${SCATGATHN}".mutect2.vcf" \
    O=$sampleID"."\${SCATGATHN}".sort.mutect2.vcf" \
    SD=$dict
  """
}

// 2.42: MuTect2_merge
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
  picard MergeVcfs I=vcf.list O=$sampleID".mutect2.merge.vcf"
  """
}

mutect2_st
  .groupTuple()
  .map { it -> tuple(it[0], it[1][0], it[2][0..-1].flatten()) }
  .set { mutect2_sm }

/* 2.43: MuTect2 Concatenate VCFs
*/
process mutct2_concstat {

  label 'med_mem'

  input:
  tuple val(caseID), val(sampleID), file(stats) from mutect2_sm

  output:
  tuple val(caseID), val(sampleID), file('*mutect2.merge.vcf.stats') into mutect2_stats

  script:
  """
  STATS=\$(ls *stats | perl -ane 'foreach \$k (@F){print "--stats \$k ";}')
  gatk MergeMutectStats --output $sampleID".mutect2.merge.vcf.stats" \$STATS
  """
}

/* 2.44: MuTect2 Concatenate VCFs
*/
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

/* 2.44: MuTect2 Contamination
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

  input:
  tuple val(caseID), val(sampleID), file(tumourbam), file(tumourbai), val(germlineID), file(germlinebam), file(germlinebai), file(mergevcf), file(statsvcf), file(readorient) from mutect2_contam_merge
  tuple file(fasta), file(fai), file(dict), file(gps), file(gpstbi), file(intlist) from Channel.value([params.fasta, params.fai, params.dict, params.gps, params.gpstbi, params.intlist])

  output:
  tuple val(caseID), file("${sampleID}.mutect2.snv_indel.pass.vcf") into mutect2_veping
  tuple val(caseID), file("${sampleID}.mutect2.raw.vcf") into mutect2_rawVcf
  file('*') into completedmutect2call

  script:
  taskmem = javaTaskmem("${task.memory}")
  """
  gatk --java-options -Xmx$taskmem \
    GetPileupSummaries \
    -I $tumourbam \
    -V $gps \
    -O $sampleID".getpileupsummaries.table" \
    -L $intlist

  gatk CalculateContamination \
    -I $sampleID".getpileupsummaries.table" \
    -O $sampleID".calculatecontamination.table"

  Rscript --vanilla ${workflow.projectDir}/bin/MuTect2_contamination.call.R $sampleID".calculatecontamination.table" $sampleID

  gatk IndexFeatureFile \
    --feature-file $mergevcf

  if [[ ! -e ${sampleID}.contamination.NA-issue.table ]]; then
    CALCONTAM="--contamination-table ${sampleID}.calculatecontamination.table"
  else
    CALCONTAM=""
  fi

  gatk --java-options -Xmx$taskmem \
    FilterMutectCalls \
    --reference $fasta \
    \$CALCONTAM \
    --interval-padding 5 \
    --output $sampleID".mutect2.FilterMutectCalls.vcf" \
    --unique-alt-read-count 3 \
    --variant $mergevcf \
    --stats $statsvcf \
    --disable-sequence-dictionary-validation true \
    --ob-priors $readorient \
    -L $intlist

  perl ${workflow.projectDir}/bin/filter_Lancet_Mutect2_Manta-Strelka2_Format.pl \
    ID=$sampleID \
    DP=14 \
    MD=2 \
    VCF=$sampleID".mutect2.FilterMutectCalls.vcf"
  """
}

/* 2.5: Manta output is a pre-req for Strelka2, so call both here
*/
process mntstr {

  label 'high_mem'

  publishDir path: "$params.outDir/cases/$caseID/manta-strelka2", overwrite: 'true', mode: "copy"

  input:
  tuple val(caseID), val(sampleID), file(tumourbam), file(tumourbai), val(germlineID), file(germlinebam), file(germlinebai) from mantastrelka2ing
  tuple file(fasta), file(fai), file(dict) from Channel.value([params.fasta, params.fai, params.dict])
  tuple file(bedgz), file(bedgztbi) from Channel.value([params.bedgz, params.bedgztbi])

  output:
  tuple val(caseID), file("${sampleID}.strelka2.snv_indel.pass.vcf") into strelka2_veping
  tuple val(caseID), file("${sampleID}.strelka2.raw.vcf") into strelka2_rawVcf
  file('*.txt') into log_mantastrelka

  script:
  """
  {
  if [[ ${params.seqlevel} == "exome" ]];then
    CR="--exome --callRegions $bedgz"
  else
    CR="--callRegions $bedgz"
  fi

  configManta.py \$CR --referenceFasta=$fasta --normalBam=$germlinebam --tumourBam=$tumourbam --runDir=manta

  manta/runWorkflow.py -m local -j ${task.cpus}

  configureStrelkaSomaticWorkflow.py \$CR --referenceFasta=$fasta --indelCandidates=manta/results/variants/candidateSmallIndels.vcf.gz --normalBam=$germlinebam --tumorBam=$tumourbam --runDir=strelka2

  strelka2/runWorkflow.py -m local -j ${task.cpus}

  ##merge into raw snv_indel
  gatk MergeVcfs -I strelka2/results/variants/somatic.snvs.vcf.gz -I strelka2/results/variants/somatic.indels.vcf.gz -O tmp.strelka2.snv_indel.vcf

  ${workflow.projectDir}/bin/manta_strelka2_rename_filter.sh  tmp.strelka2.snv_indel.vcf tmp2.strelka2.snv_indel.vcf ${sampleID} ${germlineID}

  perl ${workflow.projectDir}/bin/filter_Lancet_Mutect2_Manta-Strelka2_Format.pl \
      ID=$sampleID \
      DP=14 \
      MD=2 \
      VCF=tmp2.strelka2.snv_indel.vcf

  } 2>&1 | tee > ${sampleID}.manta-strelka2.log.txt
  """
}

/* 2.61: Lancet
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
  tuple file(fasta), file(fai), file(dict) from Channel.value([params.fasta, params.fai, params.dict])

  output:
  tuple val(caseID), val(sampleID), file('*.sort.lancet.vcf') into lancet_gt

  when:
  params.seqlevel == "exome"

  script:
  """
  SCATGATHN=\$(echo $bed | perl -ane '@s=split(/\\./);print\$s[2];')
  lancet \
    --num-threads ${task.cpus} \
    --ref $fasta \
    --bed $bed \
    --tumor $tumourbam \
    --normal $germlinebam | \
    perl -ane 'if(\$F[0]=~m/^\\#CHROM/){
      \$_=~s/TUMOR/$sampleID/;
      \$_=~s/NORMAL/$germlineID/;
      print \$_;}
    else{print \$_;}' > $sampleID"."\${SCATGATHN}".lancet.vcf"

  picard SortVcf \
    I=$sampleID"."\${SCATGATHN}".lancet.vcf" \
    O=$sampleID"."\${SCATGATHN}".sort.lancet.vcf" \
    SD=$dict
  """
}

/* 2.62: Lancet Merge
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
  picard MergeVcfs I=vcf.list O=$sampleID".lancet.merge.vcf"
  """
}

/* 2.63: Lancet Filter
*/
process lancet_filter {

  label 'med_mem'

  publishDir path: "$params.outDir/cases/$caseID/lancet"

  input:
  tuple val(caseID), val(sampleID), file(mergevcf) from lancet_merge

  output:
  tuple val(caseID), file("${sampleID}.lancet.snv_indel.pass.vcf") into lancet_veping
  tuple val(caseID), file("${sampleID}.lancet.raw.vcf") into lancet_rawVcf
  file('*') into completedlancetcall

  script:
  """
  perl ${workflow.projectDir}/bin/filter_Lancet_Mutect2_Manta-Strelka2_Format.pl \
    ID=$sampleID \
    DP=14 \
    MD=2 \
    VCF=$mergevcf
  """
}

/* 3.0: Annotate Vcfs
*/
lancet_veping
  .join( mutect2_veping )
  .join( strelka2_veping )
  .groupTuple()
  .map { it -> [it[0], it[1][0], it[2][0], it[3][0]].flatten() }
  .set { case_veping }

process vepann {

  label 'med_mem'

  publishDir path: "$params.outDir/cases/$caseID/vcf", mode: "copy", pattern: '*.vcf'

  input:
  tuple val(caseID), file(vcf1), file(vcf2), file(vcf3) from case_veping
  tuple file(fasta), file(fai), file(dict) from Channel.value([params.fasta, params.fai, params.dict])

  output:
  tuple val(caseID), file('*vep.vcf') into vepGRanges

  script:
  sampleID = "${vcf1}".split("\\.")[0]
  """
  GRCHVER=\$(ls ${params.grchvers})
  VEPVERS=\$(ls ${params.pcgrdir}/data/*/.vep/homo_sapiens/ | cut -d "_" -f2)

  for VCF in $vcf1 $vcf2 $vcf3; do
    VCFANNO=\$(echo \$VCF | sed "s/.vcf/.vep.vcf/")
    vep --dir_cache ${params.pcgrdir}/data/\$GRCHVER/.vep \
      --offline \
      --assembly \$VEPVERS \
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
      --fasta $fasta \
      --hgvs \
      --canonical \
      --ccds \
      --sift b \
      --polyphen b \
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
vepGRanges
  .join(mutect2_rawVcf)
  .join(strelka2_rawVcf)
  .join(lancet_rawVcf)
  .groupTuple()
  .map { it -> [it[0], it[1][0], it[1][1], it[1][2], it[2], it[3], it[4]] }
  .set { cons_vcfs }

process vcfGRa {

  label 'med_mem'

  publishDir "$params.outDir/cases/$caseID/consensus_vcfs"

  input:
  tuple val(caseID), file(vvcf1), file(vvcf2), file(vvcf3), file(rvcf1), file(rvcf2), file(rvcf3) from cons_vcfs
  tuple file(callR), file(funcR) from vcfGRa_Scripts

  output:
  tuple val(caseID), val(sampleID), file("${sampleID}.ALL.consensus.tab.pcgr.vcf") into pcgr_vcfs

  script:
  sampleID = "${vvcf1}".split("\\.")[0]
  """
  Rscript --vanilla $callR \
    $funcR \
    $sampleID \
    "snv_indel.pass.vep.vcf" \
    ""
  """
}

/* 3.2 Create VCF for PCGR from consensus
*/
process pcgrVcf {

  label 'low_mem'

  input:
  tuple val(caseID), val(sampleID), file(cons_vcf) from pcgr_vcfs

  output:
  tuple val(sampleID), file("${sampleID}.snv_indel.pass.pcgr.vcf") into snvpass_pcgr

  script:
  """
  NVCF="${sampleID}.snv_indel.pass.pcgr.vcf"
  cat ${workflow.projectDir}/bin/vcf42.head.txt > \$NVCF
  head -n1 $cons_vcf >> \$NVCF
  tail -n+2 $cons_vcf | sort -V >> \$NVCF
  """
}

snvpass_pcgr
  .join(facets_pcgr)
  .set { pcgr_inputs }

/* 3.3 PCGR report
* take all mutations in consensus.tab from pass.vcfs into single VCF for PCGR
*/
process pcgrreport {

  label 'low_mem'

  publishDir "$params.outDir/reports/pcgr", mode: "copy", pattern: "*html"
  publishDir "$params.outDir/cases/$caseID/pcgr", mode: "copy"

  input:
  tuple val(sampleID), file(vcf), file(jointsegs), file(ploidpur) from pcgr_inputs
  file(toml) from Channel.value([params.pcgrtoml])

  output:
  file('*') into completedPCGR

  script:
  caseID="${sampleID}".split("T_")[0]
  """
  {
  PLOIDY=""; PURITY="";
  if [[ \$(cut -f 1 $ploidpur) != "NA" ]]; then
    PLOIDY="--tumor_ploidy \$(cut -f 1 $ploidpur)"
  fi
  if [[ \$(cut -f 2 $ploidpur) != "NA" ]]; then
    PURITY="--tumor_purity \$(cut -f 2 $ploidpur)"
  fi

  VERS=\$(ls ${params.pcgrdir}/data)

  ##8 lines in an empty VCF
  LINETEST=\$(wc -l $vcf | perl -ane 'print \$F[0];')

  if [[ \$LINETEST != 8 ]]; then
    pcgr.py ${params.pcgrdir} \
      --input_vcf $vcf \
      --input_cna $jointsegs \$PLOIDY \$PURITY \
      --no-docker \
      --force_overwrite \
      ./ \
      \$VERS \
      $toml \
      $sampleID
  else
    echo "No variants"
  fi
  } 2>&1 | tee > ${sampleID}.pcgr.log.txt
  """
}

/* 4.0 Run multiQC to finalise report
*/
process mltiQC {

  label 'low_mem'

  publishDir path: "$params.outDir/reports", mode: "copy", pattern: "*html"

  input:
  file(fastps) from fastp_multiqc.collect()
  file(gtkrcls) from gtkrcl_multiqc.collect()
  file(multimetrics) from multimetrics_multiqc.collect()
  file(mrkdups) from mrkdup_multiqc.collect()

  output:
  file('*') into completedmultiqc

  script:
  """
  OUTID=\$(basename ${workflow.launchDir})".batch_somatic"
  multiqc . -i \$OUTID --tag DNA -f -c ${params.multiqcConfig}
  """
}
