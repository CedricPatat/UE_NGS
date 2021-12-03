#####################################################
#######       RNA-Seq PIPELINE SCRIPT         #######
#####################################################

# Download data
echo -e '\033[1;0;33m DOWNLOADING DATA \033[0m'
if [ ! -e initial_datas/TPrnaseq.tar.gz ];then     
    wget  http://rssf.i2bc.paris-saclay.fr/X-fer/AtelierNGS/TPrnaseq.tar.gz  # downloading RNA-Seq data
    tar -zxvf TPrnaseq.tar.gz    # data decompression
    mv *.fastq initial_datas/
    mv *.tar.gz initial_datas/
else
    echo "RNA-Seq data is already downloaded."
fi

if [ ! -e genome_data/chr18.fa ];then 
    wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/chr18.fa.gz    # downloading human genome
    gunzip chr18.fa.gz  # Decompressing the .fasta file
    mv *.fa genome_data/
else
    echo "The reference genome is already downloaded."
fi

if [ ! -e genome_data/gencode.v24lift37.basic.annotation.gtf ];then 
    wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_24/GRCh37_mapping/gencode.v24lift37.basic.annotation.gtf.gz # downloading reference genome annotation 
    gunzip gencode.v24lift37.basic.annotation.gtf.gz  # Decompressing the .gtf file
    mv *.gtf genome_data/
else
    echo "The annotation of the reference genome is already downloaded."
fi

# Indexing of the reference genome
echo -e '\033[1;0;33m INDEXING THE REFERENCE GENOME WITH STAR \033[0m'
if [ ! -e star_index/Genome ] || [ ! -e star_index/Log.out ] || [ ! -e star_index/SA ] || [ ! -e star_index/SAindex ] || [ ! -e star_index/chrLength.txt ] || [ ! -e star_index/chrName.txt ] || [ ! -e star_index/chrNameLength.txt ] || [ ! -e star_index/chrStart.txt ] || [ ! -e star_index/exonGeTrInfo.tab ] || [ ! -e star_index/exonInfo.tab ] || [ ! -e star_index/geneInfo.tab ] || [ ! -e star_index/genomeParameters.txt ] || [ ! -e star_index/sjdbInfo.txt ] || [ ! -e star_index/sjdbList.fromGTF.out.tab ] || [ ! -e star_index/sjdbList.out.tab ] || [ ! -e star_index/transcriptInfo.tab ]; then
    STAR \
        --runMode genomeGenerate \
        --runThreadN 4 \
        --genomeDir star_index \
        --genomeFastaFiles genome_data/chr18.fa \
        --sjdbGTFfile genome_data/gencode.v24lift37.basic.annotation.gtf \
        --genomeSAindexNbases 12
else
    echo "Indexing of the reference genome has already been performed."
fi

# Trimming
echo -e '\033[1;0;33m FILES ANALYSIS WITH TRIMMOMATIC \033[0m'
# Retrieving files for analysis
files=()  # list containing the names of the files
j=0     # increment index
for i in `find initial_datas/*.fastq`    # Loop saving each .fastq file name in the list
do
name=` echo $i|cut -d"/" -f2 `
files[j]=$name
j=$(($j + 1))
done
# Analysis with Trimmomatic
h=0
for i in `seq 0 $(($j-1))`  # Trimmomatic analysis loop on each file
do
name=` echo ${files[h]}|cut -d"." -f1 `
trimmomatic PE initial_datas/${files[h]} initial_datas/${files[h+1]} -baseout $name.fastq LEADING:20 TRAILING:20 MINLEN:50
h=$(($h + 2))
done
mv *P.fastq trimming_results/
mv *U.fastq trimming_results/

# Fastqc
echo -e '\033[1;0;33m QUALITY CONTROL \033[0m'
fastqc trimming_results/*P.fastq fastqc_results/  # Fastqc quality control
mv trimming_results/*fastqc* fastqc_results/

# Mapping with STAR
echo -e '\033[1;0;33m GENOME MAPPING WITH STAR \033[0m'
files_1=()  # list containing the names of the files
j=0     # increment index
for i in `find trimming_results/*1P.fastq`    # Loop saving each .fastq file name in the list
do
files_1[j]=$i
j=$(($j + 1))
done

files_2=()  # list containing the names of the files
j=0     # increment index
for i in `find trimming_results/*2P.fastq`    # Loop saving each .fastq file name in the list
do
files_2[j]=$i
j=$(($j + 1))
done
for i in `seq 0 $(($j-1))`
do
name=` echo ${files_1[i]}|cut -d"c" -f1`
name=` echo $name|cut -d"/" -f2 `
STAR \
    --runThreadN 4 \
    --outFilterMultimapNmax 1 \
    --genomeDir star_index \
    --outSAMattributes All \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix $name \
    --readFilesIn ${files_1[$i]} ${files_2[$i]}
done
mv *.bam mapping/ && mv *.out mapping/ && mv *.tab mapping/
# Sorting and indexing BAM files
for i in `find mapping/*.bam`  # Loop through all .bam files
do
samtools index $i
done

echo -e '\033[1;0;33m COMPTAGE \033[0m'
featureCounts -p -t exon -g gene_id -a genome_data/*.gtf \
    -o counts.txt mapping/*.bam

echo -e '\033[1;0;33m ENCODAGE TO HUGO \033[0m'
perl -ne 'print "$1 $2\n" if /gene_id \"(.*?)\".*gene_name \"(.*?)\"/' \
   genome_data/*.gtf | sort | uniq > encode-to-hugo.tab

# Replacement of gene names
sort counts.txt > temp1
sort encode-to-hugo.tab > temp2
join temp1 temp2 |grep "chr18" > temp3

# Final matrix generation
echo -e '\033[1;0;33m MATRIX GENERATION \033[0m'
awk '{print $13 " " $6 " " $7 " " $8 " " $9 " " $10 " " $11 " " $12}' temp3 > counts.mtx  # matrix whit HUGO nomenclature, length and the counts for the 6 individuals
mv *.txt* counts/ && mv *.tab counts/ && mv temp* counts/

echo -e '\033[1;0;32m EXECUTION COMPLETED \033[0m'
echo "The count matrix corresponds to the file 'counts.mtx'."
echo "Columns: HUGO nomenclature, length, 6 * counts."
