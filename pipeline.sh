# Pipeline ateliers NGS

echo -e '\033[1;0;33m CONNEXION A LA VM ET CREATION DU RÉPERTOIRE \033[0m'
#ssh [user@yourVM-IP-address]    # connexion ssh à la machine virtuelle
#mkdir TPRNAseq                  # Création du répertoire
#cd TPRNAseq                     # ouverture du répertoire


# Téléchargement des données 
echo -e '\033[1;0;33m RÉCUPÉRATION ET DÉCOMPRESSION DES DONNÉES \033[0m'
wget  http://rssf.i2bc.paris-saclay.fr/X-fer/AtelierNGS/TPrnaseq.tar.gz  # récupération des données
tar -zxvf TPrnaseq.tar.gz    # décompression des données


# Vérification des packages
echo -e '\033[1;0;33m VÉRIFICATION DES PACKAGES \033[0m'

dpkg -s fastqc &> /dev/null  # Fastqc
if [ $? -eq 0 ]; then
    echo "Le package fastqc est déjà installé"
else
    echo "Installation du package fastqc"
    sudo apt-get install -y fastqc  # installation de fastqc
fi

dpkg -s trimmomatic &> /dev/null  # Trimmomatic
if [ $? -eq 0 ]; then
    echo "Le package trimmomatic est déjà installé"
else
    echo "Installation du package trimmomatic"
    sudo apt-get install -y trimmomatic  # installation de trimmomatic
fi

dpkg -s rna-star &> /dev/null  # Trimmomatic
if [ $? -eq 0 ]; then
    echo "Le package RNA-STAR est déjà installé"
else
    echo "Installation du package RNA-STAR"
    sudo apt-get install -y rna-star  # installation de trimmomatic
fi


# Récupération du nom des fichiers d'intérêts
echo -e '\033[1;0;33m RÉCUPÉRATION DU NOM DES FICHIERS À ANALYSER \033[0m'
FILES=()  # liste contenant les noms des fichiers
j=0     # indice d'incrémentation
for i in `find *.fastq`    # Boucle enregistrant chaque nom fichier fastq dans la liste
do
FILES[j]=$i
j=$(($j + 1))
done


# Trimming sur les fichiers fastq
echo -e '\033[1;0;33m ANALYSE DES FICHIERS AVEC TRIMMOMATIC \033[0m'
h=0
for i in `seq 0 $(($j-1))`  # Boucle d'analyse trimmomatic sur chacun des fichiers
do
name=` echo ${FILES[h]}|cut -d"." -f1 `
java -jar Trimmomatic-0.39/trimmomatic-0.39.jar PE ${FILES[h]} ${FILES[h+1]} -baseout $name.fastq LEADING:20 TRAILING:20 MINLEN:50
h=$(($h + 2))
done

# Récupération du nom des fichiers d'intérêt après trimming
echo -e '\033[1;0;33m RÉCUPÉRATION DU NOM DES FICHIERS À ANALYSER \033[0m'
FILES=()  # liste contenant les noms des fichiers
j=0     # indice d'incrémentation
for i in `find *P.fastq`    # Boucle enregistrant chaque nom fichier issu du trimming dans la liste
do
FILES[j]=$i
j=$(($j + 1))
done

# Analyse des fichiers avec fastqc
echo -e '\033[1;0;33m ANALYSE DES FICHIERS AVEC FASTQC \033[0m'
for i in `seq 0 $(($j-1))`  # Boucle d'analyse fastqc sur chacun des fichiers
do
fastqc ${FILES[i]} 
done


# Téléchargement du génome humain et de son annotation
echo -e '\033[1;0;33m TÉLÉCHARGEMENT DU GÉNOME HUMAIN ET DE SON ANNOTATION \033[0m'
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/chr18.fa.gz    # téléchargement du génome humain 
gunzip chr18.fa.gz  # Decompression du fichier fasta

wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_24/GRCh37_mapping/gencode.v24lift37.basic.annotation.gtf.gz # téléchargement de l'annotation du génome humain
gunzip gencode.v24lift37.basic.annotation.gtf.gz  # Décompression de l'archive



# Indexation du génome avec STAR
echo -e '\033[1;0;33m INDEXATION DU GENOME AVEC STAR \033[0m'
mkdir star_index
STAR \
    --runMode genomeGenerate \
    --runThreadN 4 \
    --genomeDir star_index \
    --genomeFastaFiles chr18.fa \
    --sjdbGTFfile gencode.v24lift37.basic.annotation.gtf \
    --genomeSAindexNbases 12


# Mapping avec STAR
echo -e '\033[1;0;33m Récupération des préfixes \033[0m'

FILES_R1=()  # liste contenant les noms des fichiers
j=0     # indice d'incrémentation
for i in `find *_1P.fastq`    # Boucle enregistrant chaque nom fichier fastq dans la liste
do
FILES_R1[j]=$i
j=$(($j + 1))
done

FILES_R2=()  # liste contenant les noms des fichiers
j=0     # indice d'incrémentation
for i in `find *_2P.fastq`    # Boucle enregistrant chaque nom fichier fastq dans la liste
do
FILES_R2[j]=$i
j=$(($j + 1))
done

echo -e '\033[1;0;33m MAPPING DU GENOME AVEC STAR \033[0m'
for i in `seq 0 $(($j-1))`
do
name=` echo ${FILES_R1[i]}|cut -d"c" -f1`
STAR \
    --runThreadN 4 \
    --outFilterMultimapNmax 1 \
    --genomeDir star_index \
    --outSAMattributes All \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix $name \
    --readFilesIn ${FILES_R1[i]} ${FILES_R2[i]}
done



# Tri et indexation du fichier BAM
echo -e '\033[1;0;33m TRI ET INDEXATION DU FICHIER BAM \033[0m'

FILES_BAM=()  # liste contenant les noms des fichiers
j=0     # indice d'incrémentation
for i in `find *.bam`    # Boucle enregistrant chaque nom fichier bamls dans la liste
do
FILES_BAM[j]=$i
j=$(($j + 1))
done


for i in `seq 0 $(($j-1))`
do
samtools index ${FILES_BAM[i]}
done

# Comptage
echo -e '\033[1;0;33m COMPTAGE \033[0m'
sudo apt-get install subread
featureCounts -p -t exon -g gene_id -a *.gtf \
    -o counts.txt *.bam

# Encodage to HUGO
echo -e '\033[1;0;33m ENCODAGE TO HUGO \033[0m'
perl -ne 'print "$1 $2\n" if /gene_id \"(.*?)\".*gene_name \"(.*?)\"/' \
   *.gtf | sort | uniq > encode-to-hugo.tab

#Remplacement des noms de gènes
sort counts.txt > temp1
sort encode-to-hugo.tab > temp2
join temp1 temp2 |grep "chr18" > temp3

awk '{print $13 " " $6 " " $7 " " $8 " " $9 " " $10 " " $11 " " $12}' temp3 > counts.mtx

#Organisation des résultats et des données dans des répertoires pour une meilleur lisibilité
# echo -e '\033[1;0;33m RANGEMENT DES FICHIERS \033[0m'

# mkdir trimming_results  # Répertoire qui contiendra les fichiers issus du trimming
# mkdir fastqc_results    # Répertoire qui contiendra les fichiers de contrôle qualité
# mkdir initial_datas     # Répertoire qui contiendra les fichiers de données initiaux
# mkdir genome_datas      # Répertoire qui contiendra les fichiers de données du génome humain (chromosome 18)
# mkdir star_mapping      # Répertoire qui contiendra les fichiers du mapping STAR
# mkdir samtools_index    # Répertoire qui contiendra les fichiers .bai

# # Déplacement des fichiers issus du trimming dans le répertoire trimming_results
# for i in `find *P.fastq`    
# do
# mv $i trimming_results/$i
# done
# for i in `find *U.fastq`    
# do
# mv $i trimming_results/$i
# done

# # Déplacement des fichiers issus de l'analyse qualité dans le répertoire inital_datas
# for i in `find *fastqc.html`    
# do
# mv $i fastqc_results/$i
# done
# for i in `find *fastqc.zip`    
# do
# mv $i fastqc_results/$i
# done


# # Déplacement des fichiers de données initiales dans le répertoire inital_datas
# for i in `find *.fastq`    
# do
# mv $i initial_datas/$i
# done

# # Déplacement des fichiers du génome + annotation dans le répertoire genome_datas
# for i in `find *.fa`    
# do
# mv $i genome_datas/$i
# done
# for i in `find *.gtf`    
# do
# mv $i genome_datas/$i
# done


# # Déplacement des fichiers du mapping STAR dans le répertoire star_mapping
# for i in `find *.bam`    
# do
# mv $i star_mapping/$i
# done
# for i in `find *.out`    
# do
# mv $i star_mapping/$i
# done
# for i in `find *.tab`    
# do
# mv $i star_mapping/$i
# done


# # Déplacement des fichiers de l'indexation samtools dans le répertoire samtools_index
# for i in `find *.bai`    
# do
# mv $i samtools_index/$i
# done


echo -e '\033[1;0;32m EXECUTION TERMINÉE \033[0m'


