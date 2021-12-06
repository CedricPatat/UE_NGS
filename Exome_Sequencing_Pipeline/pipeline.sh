#####################################################
#######           PIPELINE SCRIPT             #######
#####################################################

# Downloading data
echo -e '\033[1;0;33m DOWNLOADING DATA \033[0m'
if [ ! -e patient7.tar.gz ];then     
    wget --load-cookies /tmp/cookies.txt "https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://docs.google.com/uc?export=download&id=1DM9g8OulE1ScBk-HoaREfUZs6gurtkBe' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=1DM9g8OulE1ScBk-HoaREfUZs6gurtkBe" -O patient7.tar.gz && rm -rf /tmp/cookies.txt  # downloading RNA-Seq data
    tar -zxvf patient7.tar.gz    # data decompression
    gunzip patient7.exome/*  # .fastq decompression
    #mv *.fastq initial_datas/
    #mv *.tar.gz initial_datas/
else
    echo "RNA-Seq data is already downloaded."
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


