
**Nodewalk Analysis Pipeline**

Nodewalk Analysis Pipeline (NAP) provides tools to map, analyse and identify interactors from the results generated from Nodewalk Experiment. 

Requirement: 
a.	Python3 and its packagaes and modules

Packages: Bio.SeqIO, Levenshtein

Module Bio.Seq 

Modules Time, gzip, itertools, sys, json, pysam, distance
 

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

"Probe Name":["Probe Name","chr#", "strand", probestart, probe-end, "Probe-Sequence", "GenomeVersion (eg., HG19)"]

Different probe definitions are separated by comma.

Output Files:
Five files are generated as output:
1.	Example.FragStats.tab
2.	Example.ProbeStats.tab
3.	Example.ReadStats.tab
4.	Example.UmiStats.tab
5.	Example.ValidCovStats.tab

Fragstats.tab has 24 columns. The details of the columns mainly used for analysing the interactors are:

1. SAMP: Sample name 
2. Frag: Interactor name: Format: **Genome:chr:FragStart+Length** (Eg: HG19:chr8:128207867+18243)
3. GENOME: For eg., HG19
4. CHROMO: Eg., chr1
5. FragStart: start position of interactor
6. FragEnd: end position of interactor
7. PROBE: Probe name 
8. ProbeChromo	
9. ProbePosition	
10. ProbeGenome	
11. ctpos: Approximate number of digested reads (Column No: 14)


1)	**Mapping and Pre-processing **

Paired-end reads are independently mapped using BWA to the reference genome in .fa format (the executable of BWA is provided along with the code).
The original Nodewalk pipeline utilizes a merged reference genome composed of phiX174, Drosophila (BDGP5.65), Escherichia coli K12 and human (GRCh37.75) genomes. 

The resulting SAM files are compressed and provided as an input into the pre-processing python script (PreProcess.py) along with the definition of the probe coordinates. Briefly, alignments with mapping quality greater than 10 of the second read were used to determine probe position. An extension region (extending from the probe’s end to the first restriction site) was used to discriminate valid from mis-annealing events. The total number of alignments in the first read with the proper probe extension in the second read was reported by restriction fragment.


2)	**Analysis**

The pipeline provides different measurements, ctTot, ctPos and ctUMI to quantify the original ligation event count. These parameters are present in the FragStat File. 
ctTot provides the measure of the total number of valid reads. This measurement may get affected by PCR amplification. ctPos provides the measurement of distinct digestion sites per fragment. 



The Nodewalk.sh calls BWA to align fastq.gz file and then call interactors. You may run the code using the following default settings:


**USAGE**

sh Nodewalk_V2.sh

For complete input of the arguments please write the following command:

sh Nodewalk_V2.sh help

Use of the Nodewalk pipeline

-i : Input fastq.gz directory 
-s : Output Sam directory 
-a : Reference Genome in .fa extension 
-p : Probe definition in .py extension 
-r : Restriction Enzyme 
-o : Output Stats Directory



For help and assistance: rashid.mehmood@ki.se, anita.gondor@ki.se
