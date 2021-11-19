# Directory deletion script

# This script removes all directories created by the "pipeline.sh" script. 

echo ========================================
echo -e '\033[1;0;32m      Deletion of all directories \033[0m'
echo ========================================

rm -fr fastqc_results
rm -fr genome_datas
rm -fr initial_datas
rm -fr trimming_results
rm -fr TPrnaseq.tar.gz


echo ========================================
echo -e '\033[1;0;32m          Deletion completed \033[0m'
echo ========================================
