#####################################################
#######       Virtual machine settings        #######
#####################################################

# Giving authorisations to scripts
echo -e '\033[1;0;33m SETTING UP PIPELINE SCRIPTS \033[0m'
chmod +x pipeline.sh  # allow pipeline.sh to be executed as programm

# Tools installation
echo -e '\033[1;0;33m TOOLS INSTALLATION \033[0m'
conda install -y trimmomatic bwa samtools varscan bedtools    # install tools into conda

# Create files 
echo -e '\033[1;0;33m ORGANIZATION OF ENVIRONMENT  \033[0m'
if [ ! -f initial_datas ];then     # initial_datas file containing RNA-Seq data
  mkdir initial_datas
else
    echo "initial_datas already exists"
fi
if [ ! -f index ];then     # index file containing RNA-Seq data
  mkdir index
else
    echo "trimmomatic_results already exists"
fi
if [ ! -f trimmomatic_results ];then     # trimmomatic_results file containing RNA-Seq data
  mkdir trimmomatic_results
else
    echo "trimmomatic_results already exists"
fi
if [ ! -f varscan_results ];then     # varscan_results file containing RNA-Seq data
  mkdir varscan_results
else
    echo "varscan_results already exists"
fi
