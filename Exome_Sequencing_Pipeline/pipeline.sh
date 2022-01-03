#####################################################
#######           PIPELINE SCRIPT             #######
#####################################################


# Downloading data
echo -e '\033[1;0;33m DOWNLOADING DATA \033[0m'

# Downloading patient data
if [ ! -e initial_datas/patient7.tar.gz ];then     
    wget --load-cookies /tmp/cookies.txt "https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://docs.google.com/uc?export=download&id=1DM9g8OulE1ScBk-HoaREfUZs6gurtkBe' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=1DM9g8OulE1ScBk-HoaREfUZs6gurtkBe" -O patient7.tar.gz && rm -rf /tmp/cookies.txt  # downloading RNA-Seq data
    tar -zxvf patient7.tar.gz    # data decompression
    gunzip patient7.exome/*  # .fastq decompression
    mv *.fastq initial_datas/
    mv *.tar.gz initial_datas/
else
    echo "RNA-Seq data is already downloaded."
fi

# Downloading the reference genome
if [ ! -e index/chr16.fa ];then 
    wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/chr16.fa.gz
    gunzip chr16.fa.gz
    mv chr16.fa index/chr16.fa
else
    echo "Reference genome is already downloaded."
fi

# Downloading the genome annotation
if [ ! -e gencode.v24lift37.basic.annotation.gtf ];then 
    wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_24/GRCh37_mapping/gencode.v24lift37.basic.annotation.gtf.gz
    gunzip gencode.v24lift37.basic.annotation.gtf.gz
else
    echo "Annotation is already downloaded."
fi

# Indexing of the reference genome
if [ ! -e index/chr16.fa.amb ] || [ ! -e index/chr16.fa.ann ] || [ ! -e index/chr16.fa.bwt ] || [ ! -e index/chr16.fa.fai ] || [ ! -e index/chr16.fa.pac ] || [ ! -e index/chr16.fa.sa ];then 
    bwa index -a bwtsw index/chr16.fa
else
    echo " Indexing of the reference genome has already been performed."
fi


# Trimming data
echo -e '\033[1;0;33m FILES ANALYSIS WITH TRIMMOMATIC \033[0m'
chmod +x patient7.exome/*
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
trimmomatic PE patient7.exome/${files[h]} patient7.exome/${files[h+1]} -baseout $name.fastq LEADING:20 TRAILING:20 MINLEN:50
h=$(($h + 2))
done

# Creating BWA Index



# Mapping with BWA
files=()  # list containing the names of the files
j=0     # increment index
for i in `find *_1P.fastq`    # Loop saving each 1P.fastq file name in the list
do
name=` echo $i|cut -d"_" -f1 ` 
files[j]=$name
j=$(($j + 1))
done

for i in `seq 0 $(($j-1))`  # Loop on each file
do
bwa mem -M -t 2 -A 2 -E 1 index/chr16.fa ${files[$i]}_1P.fastq ${files[$i]}_2P.fastq > outputs/${files[$i]}
done

# Processing SAM files


files=()  # list containing the names of the files
j=0     # increment index
for i in `find outputs/*.sam`    # Loop saving each 1P.fastq file name in the list
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
samtools sort $file.bam > $file.bai
#Index bam file
samtools index $file.bam

#Convert to Mpileup
samtools mpileup -B -A -f index/chr16.fa  $file.bai > $file.mpileup

done


# Calling somatic variants with Varscan


path_to_normal_mpileup=`find outputs/*-N-*.mpileup`
path_to_tumor_mpileup=`find outputs/*-T-*.mpileup`
output_name='TCRBOA7'  # A CHANGER


varscan somatic path_to_normal_mpileup path_to_tumor_mpileup varscan_results/$output_name


# Basic VCF Annotation

for i in `find varscan_results/*`
do

name=`echo $i|cut -d"/" -f2 `

grep -i 'somatic' $i > filtered/$name
awk '{OFS="\t"; if (!/^#/){print $1,$2-1,$2,$4"/"$5,"+"}}' \
   filtered/$name > bed/$name

bedtools intersect -a *.gtf -b bed/$name > intersect/$name 
grep '\sgene\s' intersect/$name  | awk '{print " " $1 " " $4 " " $5 " " $16}'

done

