# Pipeline ateliers NGS

echo -e '\033[1;0;33m CONNEXION A LA VM ET CREATION DU RÉPERTOIRE \033[0m'
#ssh [user@yourVM-IP-address]    # connexion ssh à la machine virtuelle
#mkdir TPRNAseq                  # Création du répertoire
#cd TPRNAseq                     # ouverture du répertoire


# Téléchargement des données 
echo -e '\033[1;0;33m RÉCUPÉRATION ET DÉCOMPRESSION DES DONNÉES \033[0m'
#wget  http://rssf.i2bc.paris-saclay.fr/X-fer/AtelierNGS/TPrnaseq.tar.gz  # récupération des données
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


# Récupération du nom des fichiers d'intérêts
echo -e '\033[1;0;33m RÉCUPÉRATION DU NOM DES FICHIERS À ANALYSER \033[0m'
FILES=()  # liste contenant les noms des fichiers
j=0     # indice d'incrémentation
for i in `find *.fastq`    # Boucle enregistrant chaque nom fichier fastq dans la liste
do
FILES[j]=$i
j=$(($j + 1))
done


# trimming sur les fichiers fastq
echo -e '\033[1;0;33m ANALYSE DES FICHIERS AVEC TRIMMOMATIC \033[0m'
h=0
for i in `seq 0 $(($j-1))`  # Boucle d'analyse trimmomatic sur chacun des fichiers
do
name=` echo ${FILES[h]}|cut -d"." -f1 `
java -jar /home/cpatat/Telechargements/Trimmomatic-0.39/trimmomatic-0.39.jar PE ${FILES[h]} ${FILES[h+1]} -baseout $name.fastq LEADING:20 TRAILING:20 MINLEN:50
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
rm -r chr18.fa.gz   # Suppression de l'archive

wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_24/GRCh37_mapping/gencode.v24lift37.basic.annotation.gtf.gz # téléchargement de l'annotation du génome humain
gunzip gencode.v24lift37.basic.annotation.gtf.gz  # Décompression de l'archive
rm -r gencode.v24lift37.basic.annotation.gtf.gz   # Suppression de l'archive



# Organisation des résultats et des données dans des répertoires pour une meilleur lisibilité
echo -e '\033[1;0;33m RANGEMENT DES FICHIERS \033[0m'

mkdir trimming_results  # Répertoire qui contiendra les fichiers issus du trimming
mkdir fastqc_results    # Répertoire qui contiendra les fichiers de contrôle qualité
mkdir inital_datas      # Répertoire qui contiendra les fichiers de données initiaux
mkdir genome_datas      # Répertoire qui contiendra les fichiers de données du génome humain (chromosome 18)

# Déplacement des fichiers issus du trimming dans le répertoire trimming_results
for i in `find *P.fastq`    
do
mv $i trimming_results/$i
done
for i in `find *U.fastq`    
do
mv $i trimming_results/$i
done

# Déplacement des fichiers issus de l'analyse qualité dans le répertoire inital_datas
for i in `find *fastqc.html`    
do
mv $i fastqc_results/$i
done
for i in `find *fastqc.zip`    
do
mv $i fastqc_results/$i
done


# Déplacement des fichiers de données initiales dans le répertoire inital_datas
for i in `find *.fastq`    
do
mv $i inital_datas/$i
done

# Déplacement des fichiers du génome + annotation dans le répertoire genome_datas
for i in `find *.fa`    
do
mv $i genome_datas/$i
done
for i in `find *.gtf`    
do
mv $i genome_datas/$i
done

echo -e '\033[1;0;32m EXECUTION TERMINÉE \033[0m'

