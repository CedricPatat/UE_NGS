#####################################################
#######           PIPELINE SCRIPT             #######
#####################################################

# Downloading data
echo -e '\033[1;0;33m DOWNLOADING DATA \033[0m'
if [ ! -e patient7.tar.gz ];then     
    wget --load-cookies /tmp/cookies.txt "https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://docs.google.com/uc?export=download&id=1DM9g8OulE1ScBk-HoaREfUZs6gurtkBe' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=1DM9g8OulE1ScBk-HoaREfUZs6gurtkBe" -O patient7.tar.gz && rm -rf /tmp/cookies.txt  # downloading RNA-Seq data
    tar -zxvf patient7.tar.gz    # data decompression
    #mv *.fastq initial_datas/
    #mv *.tar.gz initial_datas/
else
    echo "RNA-Seq data is already downloaded."
fi

# Trimming data
