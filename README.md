
**Nodewalk Analysis Pipeline**

Nodewalk Analysis Pipeline (NAP) provides tools to map, analyse and identify interactors from the results generated from Nodewalk Experiment. 

Requirement: 
a.	Python (version)
**c.**	


Input Datasets: 
a.	The FASTQ files generated from the Nodewalk wet-lab protocol are the starting inputs for the Nodewalk Analysis Pipeline. 
Supportive file formats : Compressed FASTQ files (fastq.gz)
To detect the pair-end files, a naming convention is followed: The pairs should have suffix
_R1 and _R2 followed by the extension .fastq.gz.

For eg.: 
example1_R1.fastq.gz
example1_R2.fastq.gz

b.	Reference Genome 
c.	Sequence of recognition site for restriction enzyme used. For eg., AAGCTT is the recognition site for HindIII
d.	Probe definition file. The format is as follows:
"Probe Name":	["Probe Name","chr#", "strand", probestart, probe-end, "Probe-Sequence", "GenomeVersion (eg., HG19)"]

Output Files:
Five files are generated as output:
1.	Example.FragStats.tab
2.	Example.ProbeStats.tab
3.	Example.ReadStats.tab
4.	Example.UmiStats.tab
5.	Example.ValidCovStats.tab

Fragstats.tab has the following collumns:
SAMP: Sample name
Frag: Interactor name: Format:
GENOME: 
CHROMO	
FragStart	
FragEnd	
PROBE 
ProbeChromo	
ProbePosition	
ProbeGenome	
LT10	
LT30	
GT30	
ctpos	
ctPos5	
ctPos10	
ctumi	
ctUmi5	
ctUmi10	
maxPrbMatches	
maxPrbOverPos	
maxPrbTp	
Positions	
Umis


1)	**Mapping and Pre-processing **

Paired-end reads are independently mapped using bwa version (##) to the reference genome in .fa format. 
The original Nodewalk pipeline utilizes a merged reference genome composed of phiX174, Drosophila (BDGP5.65), Escherichia coli K12 and human (GRCh37.75) genomes. 

The resulting sam files are compressed and provided as an input into the pre-processing python script (PreProcess.py) along with the definition of the probe coordinates. Briefly, alignments with mapping quality greater than 10 of the second read were used to determine probe position. An extension region (extending from the probe’s end to the first restriction site) was used to discriminate valid from mis-annealing events. The total number of alignments in the first read with the proper probe extension in the second read was reported by restriction fragment.


2)	**Analysis**

The pipeline provides different measurements, ctTot, ctPos and ctUMI to quantify the original ligation event count. These parameters are present in the FragStat File. 
ctTot provides the measure of the total number of valid reads. This measurement may get affected by PCR amplification. ctPos provides the measurement of distinct digestion sites per fragment. 



nodewalk pipeline to analyze the data using our in lab developed protocol The Nodewalk.sh is containing bash script to align your multiple compressed fastq.gz file and then analyze the Nodewalk analysis. you may run the code using the following default settings:


**USAGE**

sh Nodewalk.sh

For complete input of the arguments please write the following command:

sh Nodewalk.sh help

Use of the Nodewalk pipeline

-i : Input fastq.gz directory 
-s : Output Sam directory 
-a : Reference Genome in .fa extension 
-p : Prob defination in .py extension 
-r : Restriction Enzyme 
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

For help and assistance: rashid.mehmood@ki.se
