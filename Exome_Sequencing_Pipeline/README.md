<div align="center"><h1>Exome Sequencing pipeline</h1></div>
<br>

<div align="justify">
  <p>
    The scripts contained in this repository are used to identify somatic mutation.
  </p>
</div>

<div align="left"><h2>How to use the pipeline?</h2></div>
<div align="justify">
  <p>
    
   The proposed pipeline uses Conda (Bioconda) for different processes. You must therefore have installed conda beforehand. <br>
   To retrieve the scripts, clone the repository: `$ git clone https://github.com/CedricPatat/UE_NGS.git`. <br>
   Open the folder containing the scripts: `$ cd Exome_Sequencing_Pipeline`. <br>
   Then activate conda: `$ conda activate`. <br>
   And allow initialization.sh to be executed: `$ chmod +x initialization.sh`. <br>
   Initialize the machine: `$ ./initialization.sh`. <br>
   Run the pipeline: `$ ./pipeline.sh`.
      
  </p>
</div>

<div align="left"><h2>Functioning</h2></div>
<div align="justify">
  <p>
    
   The `initialization.sh` script allows the pipeline script to be run as a program. It also creates the directories in which the files created during the 
    execution of the pipeline will be stored. <br>
    <br>
   The `pipeline.sh` script runs the pipeline. It downloads the genome of two patients, one with a tumor and a control. an analysis with TRIMMOMATIC is then
   carried out then the mapping of the data is done with BWA. The SAMTOOLS tool allows you to convert the files to the right format and then analyze them with 
   BEDTOOLS in order to identify somatic mutations. The output is a file of 4 columns corresponding to the genes exhibiting somatic variations.
 
    
  To modify the working data, just change the links in the `pipeline.sh` script.
    
  </p>
</div>
