#!/bin/bash
 
#####################################  Soft directory ###################################
read1="/xxx/Reseq/Test/sen.guo/lee/exome-data/A-1_1.fq.gz"
read2="/xxx/Reseq/Test/sen.guo/lee/exome-data/A-1_2.fq.gz"
outdir="/xxx/Reseq/Test/sen.guo/lee"
sample="A-1"
adaptor1="GATCGGAAGAGCACACGTCT"
adaptor2="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
cutadapt="/xxx/app/python-package/bin"
bwa="/xxx/app/bwa-0.7.12"
refgenome="/xxx/Public/Database/hg19/ucsc.hg19.fasta"
samtools="/xxx/app/samtools-1.5/bin"
picard="/xxx/app/picard-tools-2.10.3"
GATK="/xxx/app/gatk-3.7.0"
GATK_db="/xxx/Public/Database/gatk_hg19"
java="/xxx/app/jdk18/bin"
bed="/xxx/script/Exome_Reseq/Exome_BED/SeqCapEZ_Exome_v3.0_Design_Annotation_files/exome.bed"
outbed="/xxx/script/Exome_Reseq/Exome_BED/SeqCapEZ_Exome_v3.0_Design_Annotation_files/out.exome.bed"
Ref_prefix="hg19"
annovar="/xxx/app/annovar-latest"
Annovar_db="/xxx/Public/Database/annovardb/humandb"
 
 
#################################### set up log file #####################################################
 
log="$outdir/WES.log";
echo "**************************" >> $log;
echo "command $0 \n" >> $log;
echo "Run begin: `date`" >> $log;
echo "**************************" >> $log; 
 
 
############################################## Remove adapter by Cutadapt ###################################
echo "cutadapt begin: `date`" >> $log;
mkdir -p $outdir/01_cutadapt;
$cutadapt/cutadapt --format=fastq -O 1 -a $adaptor1 -A $adaptor2 -e 0.1 -q 20,20 -m 75 --max-n=0.1 -o $outdir/01_cutadapt/A-1_1.clean.fq.gz -p $outdir/01_cutadapt/A-1_2.clean.fq.gz $read1 $read2;
echo "cutadapt end: `date`" >> $log;
 
############################################# Mapping by BWA ####################################
echo "bwa begin: `date`" >> $log;
mkdir -p $outdir/02_bwa;
$bwa/bwa mem  -R '@RG\tID:IDa\tSM:A-1\tPL:ILLUMINA' $refgenome $outdir/01_cutadapt/A-1_1.clean.fq.gz $outdir/01_cutadapt/A-1_2.clean.fq.gz > $outdir/02_bwa/A-1.sam;
 
$samtools/samtools view -bt $refgenome.fai $outdir/02_bwa/A-1.sam  > $outdir/02_bwa/A-1.bam;
 
$samtools/samtools sort $outdir/A-1.bam -o $outdir/02_bwa/A-1.sorted.bam;
echo "bwa end: `date`" >> $log;
############################################# Remove Duplicates by Picard  ##################################
echo "picard begin: `date`" >> $log;
 
mkdir -p $outdir/03_picard;
$java/java -Xmx20g -jar $picard/picard.jar MarkDuplicates ASSUME_SORTED=true REMOVE_SEQUENCING_DUPLICATES=true INPUT=$outdir/02_bwa/A-1.sorted.bam OUTPUT=$outdir/03_picard/A-1.sorted.rmdup.bam METRICS_FILE=$outdir/03_picard/A-1.sorted.rmdup.metrics;
 
$samtools/samtools index $outdir/03_picard/A-1.sorted.rmdup.bam ;
echo "picard end: `date`" >> $log;
 
############################################ GATK call VCF ##################################################
echo "GATK begin: `date`" >> $log;
 
mkdir -p $outdir/04_vcf;
$java/java -Xmx10g -jar $GATK/GenomeAnalysisTK.jar -nt 5  -T RealignerTargetCreator -R $refgenome  -I $outdir/03_picard/A-1.sorted.rmdup.bam -o $outdir/04_vcf/A-1.sorted.rmdup.bam.intervals -known $GATK_db/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf -known $GATK_db/1000G_phase1.indels.hg19.sites.vcf;
 
$java/java -Xmx10g -jar $GATK/GenomeAnalysisTK.jar  -T IndelRealigner -R $refgenome -targetIntervals $outdir/04_vcf/A-1.sorted.rmdup.bam.intervals -I $outdir/03_picard/A-1.sorted.rmdup.bam -o $outdir/04_vcf/A-1.realigned.bam -known $GATK_db/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf -known $GATK_db/1000G_phase1.indels.hg19.sites.vcf;
 
$java/java -Xmx10g -jar $GATK/GenomeAnalysisTK.jar -nct 5  -T BaseRecalibrator -R $refgenome  -I $outdir/04_vcf/A-1.realigned.bam -o $outdir/04_vcf/A-1.recal_data.grp -knownSites $GATK_db/dbsnp_138.hg19.vcf -knownSites $GATK_db/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf     -knownSites $GATK_db/1000G_phase1.indels.hg19.sites.vcf;
 
$java/java -Xmx10g -jar $GATK/GenomeAnalysisTK.jar -nct 5 -T PrintReads -R $refgenome  -I $outdir/04_vcf/A-1.realigned.bam  -BQSR $outdir/04_vcf/A-1.recal_data.grp -o $outdir/04_vcf/A-1.recal.bam;
 
 
$java/java -Xmx50g -jar $GATK/GenomeAnalysisTK.jar -nct 8 -T HaplotypeCaller -R $refgenome -I $outdir/04_vcf/A-1.recal.bam -stand_call_conf 30 -D $GATK_db/dbsnp_138.hg19.vcf -o $outdir/04_vcf/A-1.raw.vcf;
echo "GATK end : `date`" >> $log;
 
############################################# Recalibration of the vcf #######################################
echo "GATK VCF recalibration begin: `date`" >> $log;
 
mkdir -p $outdir/05_vcf_recal;
 
java -Xmx50g -jar $GATK/GenomeAnalysisTK.jar -T VariantRecalibrator -R $refgenome -mode SNP --maxGaussians 4 -input $outdir/04_vcf/A-1.raw.vcf -resource:hapmap,known=false,training=true,truth=true,prior=15.0 $GATK_db/hapmap_3.3.hg19.sites.vcf -resource:omni,known=false,training=true,truth=true,prior=12.0 $GATK_db/1000G_omni2.5.hg19.sites.vcf -resource:1000G,known=false,training=true,truth=false,prior=10.0 $GATK_db/1000G_phase1.snps.high_confidence.hg19.sites.vcf -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $GATK_db/dbsnp_138.hg19.vcf  -an QD -an FS -an MQ  -an MQRankSum -an ReadPosRankSum -an SOR -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 -recalFile $outdir/05_vcf_recal/A-1.snp.recal -tranchesFile $outdir/05_vcf_recal/A-1.snp.tranches -rscriptFile $outdir/05_vcf_recal/A-1.snp.plots.R;
 
java -Xmx50g -jar $GATK/GenomeAnalysisTK.jar  -T VariantRecalibrator -R $refgenome -mode INDEL --maxGaussians 4 -input $outdir/04_vcf/A-1.raw.vcf -resource:mills,known=true,training=true,truth=true,prior=12.0 $GATK_db/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $GATK_db/dbsnp_138.hg19.vcf  -an QD -an FS -an SOR -an MQ -an MQRankSum -an ReadPosRankSum -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 -recalFile $outdir/05_vcf_recal/A-1.indel.recal -tranchesFile $outdir/05_vcf_recal/A-1.indel.tranches    -rscriptFile $outdir/05_vcf_recal/A-1.indel.plots.R;
 
java -Xmx50g -jar $GATK/GenomeAnalysisTK.jar  -T ApplyRecalibration -R $refgenome -mode SNP -input $outdir/04_vcf/A-1.raw.vcf -tranchesFile $outdir/05_vcf_recal/A-1.snp.tranches -recalFile $outdir/05_vcf_recal/A-1.snp.recal -o $outdir/05_vcf_recal/A-1.t90.SNPrecal.vcf --ts_filter_level 99.0;    
 
 
java -Xmx50g -jar $GATK/GenomeAnalysisTK.jar  -T ApplyRecalibration -R $refgenome -mode INDEL -input $outdir/05_vcf_recal/A-1.t90.SNPrecal.vcf -tranchesFile $outdir/05_vcf_recal/A-1.indel.tranches -recalFile $outdir/05_vcf_recal/A-1.indel.recal -o $outdir/05_vcf_recal/A-1.t90.bothrecal.vcf --ts_filter_level 99.0;
 
 
java -Xmx50g -jar $GATK/GenomeAnalysisTK.jar -T VariantFiltration -R $refgenome --variant $outdir/05_vcf_recal/A-1.t90.bothrecal.vcf -o $outdir/05_vcf_recal/A-1.filter.vcf --filterExpression " DP<=6 || QD<2.00 || MQ<40.00" --filterName "lowQual"; 
echo "GATK VCF recalibration end: `date`" >> $log;
 
##################################################### Annovar annotation of the VCF file ######################
echo "annovar  begin: `date`" >> $log;
 
mkdir -p $outdir/06_annovar;
$annovar/convert2annovar.pl -format vcf4 $outdir/05_vcf_recal/A-1.filter.vcf -includeinfo > $outdir/06_annovar/A-1.exome.vcf4;
 
$annovar/table_annovar.pl $outdir/06_annovar/A-1.exome.vcf4  $Annovar_db -buildver $Ref_prefix -remove -out $outdir/06_annovar/A-1 -protocol refGene,cytoBand,esp6500siv2_all,1000g2015aug_all,1000g2015aug_afr,1000g2015aug_eas,1000g2015aug_eur,avsnp147,dbnsfp30a,exac03,gnomad_genome,clinvar_20170130,cosmic70 -operation g,r,f,f,f,f,f,f,f,f,f,f,f -nastring . ;
echo "annovar  end: `date`" >> $log;
echo "Run end: `date`" >> $log;
