#####################################################
#######           PIPELINE SCRIPT             #######
#####################################################


# Downloading data
echo -e '\033[1;0;33m DOWNLOADING DATA \033[0m'

# Downloading patient data
if [ ! -e patient7.tar.gz ] || [ ! -e patient7.exome/TCRBOA7-N-WEX-chr16_r1F.fastq ] || [ ! -e patient7.exome/TCRBOA7-N-WEX-chr16_r2F.fastq ] || [ ! -e patient7.exome/TCRBOA7-T-WEX-chr16_r1F.fastq ] || [ ! -e patient7.exome/TCRBOA7-T-WEX-chr16_r2F.fastq ];then     
    wget --load-cookies /tmp/cookies.txt "https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://docs.google.com/uc?export=download&id=1DM9g8OulE1ScBk-HoaREfUZs6gurtkBe' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=1DM9g8OulE1ScBk-HoaREfUZs6gurtkBe" -O patient7.tar.gz && rm -rf /tmp/cookies.txt  # downloading data
    tar -zxvf patient7.tar.gz    # data decompression
    gunzip patient7.exome/*      # .fastq decompression
else
    echo "RNA-Seq data is already downloaded."
fi

# Downloading the reference genome
if [ ! -e index/chr16.fa ];then 
    wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/chr16.fa.gz  # download the reference chromosome
    gunzip chr16.fa.gz             # reference chromosome file decompression
    mv chr16.fa index/chr16.fa     # moving the reference chromosome file to index
else
    echo "Reference genome is already downloaded."
fi

# Downloading the genome annotation
if [ ! -e annotation/gencode.v24lift37.basic.annotation.gtf ];then 
    wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_24/GRCh37_mapping/gencode.v24lift37.basic.annotation.gtf.gz   # download the annotation
    gunzip *.gtf.gz            # annotation file decompression
    mv *.gtf annotation/   # moving the annotation file to index
else
    echo "Annotation is already downloaded."
fi

echo -e '\033[1;0;33m INDEXING THE REFERENCE GENOME \033[0m'
# Indexing of the reference genome
if [ ! -e index/chr16.fa.amb ] || [ ! -e index/chr16.fa.ann ] || [ ! -e index/chr16.fa.bwt ] || [ ! -e index/chr16.fa.fai ] || [ ! -e index/chr16.fa.pac ] || [ ! -e index/chr16.fa.sa ];then 
    bwa index -a bwtsw index/chr16.fa   # indexing of the reference chromosome with bwa index
else
    echo " Indexing of the reference genome has already been performed."
fi

# Trimming data with trimmomatic
echo -e '\033[1;0;33m DATA TRIMMING WITH TRIMMOMATIC  \033[0m'
chmod +x patient7.exome/*  # gives permissions to files in the directory
# Retrieving files for analysis
files=()  # list containing the names of the files
j=0     # increment index
for i in `find patient7.exome/*.fastq`    # Loop saving each .fastq file name in the list
do
name=` echo $i|cut -d"/" -f2 `
files[j]=$name
j=$(($j + 1))
done
# Analysis with Trimmomatic
h=0
for i in `seq 0 $(($j-1))`  # Trimmomatic analysis loop on each file
do
name=` echo ${files[h]}|cut -d"_" -f1 `
trimmomatic PE patient7.exome/${files[h]} patient7.exome/${files[h+1]} -baseout trimmomatic_results/$name.fastq LEADING:20 TRAILING:20 MINLEN:50
h=$(($h + 2))
done


# Mapping with BWA
echo -e '\033[1;0;33m MAPPING WITH BWA  \033[0m'
files=()  # list containing the names of the files
j=0     # increment index
for i in `find trimmomatic_results/*_1P.fastq`    # Loop saving each 1P.fastq file name in the list
do
name=` echo $i|cut -d"/" -f2 ` 
name=` echo $name|cut -d"_" -f1 `
files[j]=$name
j=$(($j + 1))
done

for i in `seq 0 $(($j-1))`  # Loop on each file
do
bwa mem -M -t 2 -A 2 -E 1 index/chr16.fa trimmomatic_results/${files[$i]}_1P.fastq trimmomatic_results/${files[$i]}_2P.fastq > output_results/${files[$i]}.sam
done

# Processing SAM files
echo -e '\033[1;0;33m PROCESSING SAM FILES WITH SAMTOOLS \033[0m'

files=()  # list containing the names of the files
j=0     # increment index
for i in `find output_results/*.sam`    # Loop saving each 1P.fastq file name in the list
do
name=` echo $i|cut -d"." -f1 ` 
files[j]=$name
j=$(($j + 1))
done


for i in `seq 0 $(($j-1))`  # Loop on each file
do
file=${files[$i]}

#Sam 2 Bam
samtools view -S -b $file.sam  > $file.bam
# flagstats
samtools flagstat $file.bam
#Sort Bam
samtools sort $file.bam > $file.sorted_bam
#Index bam file
samtools index $file.bam

#Convert to Mpileup
samtools mpileup -B -A -f index/chr16.fa  $file.sorted_bam > $file.mpileup

done


# Calling somatic variants with Varscan
echo -e '\033[1;0;33m CALLING SOMATIC VARIANTS WITH VARSCAN  \033[0m'

path_to_normal_mpileup=`find output_results/*-N-*.mpileup`
path_to_tumor_mpileup=`find output_results/*-T-*.mpileup`
output_name='TCRBOA7'  # A CHANGER


varscan somatic path_to_normal_mpileup path_to_tumor_mpileup varscan_results/$output_name


# Basic VCF Annotation
echo -e '\033[1;0;33m BASIC VCF ANNOTATION  \033[0m'

for i in `find varscan_results/*`
do

name=`echo $i|cut -d"/" -f2 `

grep -i 'somatic' $i > filtered_vcf/$name
awk '{OFS="\t"; if (!/^#/){print $1,$2-1,$2,$4"/"$5,"+"}}' \
   filtered_vcf/$name > bed_results/$name

bedtools intersect -a annotation/*.gtf -b bed_results/$name > intersect_results/$name 
grep '\sgene\s' intersect_results/$name  | awk '{print " " $1 " " $4 " " $5 " " $16}'

done

