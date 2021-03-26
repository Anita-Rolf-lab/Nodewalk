# Nodewalk
nodewalk pipeline to analyze the data using our in lab developed protocol
The Nodewalk.sh is containing bash script to align your multiple compressed fastq.gz file and then analyze the Nodewalk analysis. 
you may run the code using the following default settings:

  sh Nodewalk.sh

For complete input of the arguments please write the following command:

  sh Nodewalk.sh help
  
   Use of the Nodewalk pipeline

-i : Input fastq.gz directory
-s : Output Sam directory
-a : Reference Genome in .fa extension
-p : Prob defination in .py extension
-r : Restricted Enzyme
-o : Output Stats Directory


Run with the detault settings: nw.sh
Detail of folder structure in detault settings
-i: UPLOADS/NW_ST1
-s: UPLOADS/NW_ST2
-a: Genomes/HG19_BDGP5_PhiX_K12.fa
-p: PROBE_DEFS/NWProbeDef171116.py
-r: AAGCTT
-o: UPLOADS/STATS

  The pipeline will not create any directory, please create your exact directory structure to execute the pipeline.
  
  For help and assistance:
  rashid.mehmood@ki.se
  
