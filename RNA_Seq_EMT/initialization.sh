#####################################################
#######       Virtual machine settings        #######
#####################################################

# Giving authorisations to scripts
echo -e '\033[1;0;33m SETTING UP PIPELINE SCRIPTS \033[0m'
chmod +x pipeline.sh  # allow pipeline.sh to be executed as programm

# Tools installation
echo -e '\033[1;0;33m TOOLS INSTALLATION \033[0m'
conda install -y fastqc samtools trimmomatic subread star        # install tools into conda

# Create files 
echo -e '\033[1;0;33m ORGANIZATION OF ENVIRONMENT  \033[0m'
if [ ! -f initial_datas ];then     # initial_datas file containing RNA-Seq data
  mkdir initial_datas
else
    echo "initial_datas already exists"
fi
if [ ! -f genome_data ];then     # genome_data file containing reference genome and annotation
  mkdir genome_data
else
    echo "genome_data already exists"
fi
if [ ! -f trimming_results ];then     # trimming_results file containing trimming results
  mkdir trimming_results
else
    echo "trimming_results already exists"
fi
if [ ! -f fastqc_results ];then     # fastqc_results file containing fastqc results
  mkdir fastqc_results
else
    echo "fastqc_results already exists"
fi
if [ ! -f mapping ];then     # mapping file containing mapping and datas index 
  mkdir mapping
else
    echo "mapping already exists"
fi
if [ ! -f star_index ];then     # star_index file containing genome index
  mkdir star_index
else
    echo "star_index already exists"
fi

if [ ! -f counts ];then     # counts file containing counts matrix before HUGO encoding
  mkdir counts
else
    echo "counts already exists"
fi

 
