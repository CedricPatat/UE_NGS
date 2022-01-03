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
if [ ! -f index ];then     # index file containing RNA-Seq data
  mkdir index
else
    echo "index already exists"
fi
if [ ! -f annotation ];then     # annotation file containing RNA-Seq data
  mkdir annotation
else
    echo "annotation already exists"
fi
if [ ! -f trimmomatic_results ];then     # trimmomatic_results file containing RNA-Seq data
  mkdir trimmomatic_results
else
    echo "trimmomatic_results already exists"
fi
if [ ! -f output_results ];then     # output_results file containing RNA-Seq data
  mkdir output_results
else
    echo "output_results already exists"
fi
if [ ! -f varscan_results ];then     # varscan_results file containing RNA-Seq data
  mkdir varscan_results
else
    echo "varscan_results already exists"
fi
if [ ! -f filtered_vcf ];then     # filtered_vcf file containing RNA-Seq data
  mkdir filtered_vcf
else
    echo "filtered_vcf already exists"
fi
if [ ! -f bed_results ];then     # bed_results file containing RNA-Seq data
  mkdir bed_results
else
    echo "bed_results already exists"
fi
if [ ! -f intersect_results ];then     # intersect_results file containing RNA-Seq data
  mkdir intersect_results
else
    echo "intersect_results already exists"
fi
