<div align="center"><h1>RNA-Seq analysis pipeline</h1></div>
<br>
<div align="justify">
  <p>
    The scripts contained in this repository are used to produce a counting matrix from RNA-Seq .fastq files
  </p>
</div>

<div align="left"><h2>How to use the pipeline?</h2></div>
<div align="justify">
  <p>
    
   The proposed pipeline uses Conda (Bioconda) for different processes. You must therefore have installed conda beforehand. <br>
   To retrieve the scripts, clone the repository: `$ git clone https://github.com/CedricPatat/UE_NGS.git`. <br>
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
   The `pipeline.sh` script runs the pipeline. It downloads RNA-Seq data and reference data and indexes the reference genome only once. An analysis with 
    Trimmomatic is then carried out followed by a quality control. Finally, the STAR and SUBREAD tools are used to map the reads and index them. 
    FEATURECOUNTS allows to obtain a matrix of counts which will then be modified so as to obtain the file 'counts.mtx' containing 8 columns:
   * name of the genes (HUGO format)
   * gene length
   * sample 1 counts
   * sample 2 counts
   * sample 3 counts
   * sample 4 counts
   * sample 5 counts
   * sample 6 counts
 
    
  To modify the working data, just change the links in the `pipeline.sh` script.
    
  </p>
</div>
