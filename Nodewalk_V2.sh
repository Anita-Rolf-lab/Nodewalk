while getopts "i:s:a:o:h:p:r:" option; do
  case $option in
    i ) inputDir=$OPTARG
    ;;
    s ) outputDirSam=$OPTARG
    ;;
 a ) reffa=$OPTARG
    ;;
 o ) outputDirStats=$OPTARG
    ;;
 h ) hep=$OPTARG
    ;;
 p ) pdef=$OPTARG
    ;;
 r ) renzyme=$OPTARG
    ;;
  esac
done

if [ -z "$inputDir" ]; then inputDir="UPLOADS/NW_ST1"; fi

if [ -z "$outputDirSam" ]; then outputDirSam="$inputDir/NW_ST2"; fi
if [ -z "$reffa" ]; then reffa="Genomes/HG19_BDGP5_PhiX_K12.fa"; fi

if [ -z "$outputDirStats" ]; then  outputDirStats="$inputDir/STATS"; fi


if [ -z "$pdef" ]; then pdef="PROBE_DEFS/NWProbeDef171116.py"; fi
if [ -z "$renzyme" ]; then renzyme="AAGCTT"; fi

if [  -n "$1" ]
 then
	VAR1="Linusdxize"
	VAR2="Linuxize"
	if [ "$1" = "help" ]
	then
		echo "\n Use of the Nodewalk pipeline\n\n Mandatory "
		echo "\t-i : Input fastq.gz directory"
		echo "\t-a : Reference Genome in .fa extension"

		echo "\n Optional"
		echo "\t-s : Output Sam directory, Optional"
		echo "\t-r : Restricted Enzyme"
		echo "\t-p : Prob defination in .py extension"
		echo "\t-o : Output Stats Directory, Optional"	
			
		
		echo "\n\nRun with the detault settings: nw.sh"
		echo "Detail of folder structure in detault settings"
		echo "-i: $inputDir"
		echo "-s: $outputDirSam"
		echo "-a: $reffa"		
		echo "-p: $pdef"
		echo "-r: $renzyme"
		echo "-o: $outputDirStats"
		exit
	
	fi

fi

if [ -d "$outputDirSam" ]; then mkdir "$inputDir/NW_ST2"; outputDirSam="$inputDir/NW_ST2"; fi
if [ -d "$outputDirStats" ]; then mkdir "$inputDir/STATS"; outputDirStats="$inputDir/STATS"; fi

if [ ! -d $outputDirSam ] 
then
  mkdir $outputDirSam
fi
if [ ! -d $outputDirStats ] 
then
  mkdir $outputDirStats

fi


echo "\n*************Provided/detaul Settings************************\n"
echo "Input fastq.gz directory: $inputDir"
echo "Output Sam directory: $outputDirSam"
echo "Input reference genome .fa : $reffa"

echo "Prob Definition: $pdef"
echo "Restricted Enzyme : $renzyme"

echo "Output Stats Directory :$outputDirStats"

echo "\n\nValidating Provided Directory and files\n\n";
if [ ! -d "$inputDir" ]; then echo "Error, Directory of input fastq.gz files  '$inputDir' DOES NOT exists.\n";  exit 9999; fi
if [ ! -d "$outputDirSam" ]; then echo "Error, Directory of output sam files  '$outputDirSam' DOES NOT exists.\n";  exit 9999; fi
if [ ! -d "$outputDirStats" ]; then echo "Error, Directory of output STATs   '$outputDirStats' DOES NOT exists.\n";  exit 9999; fi
if [ ! -f "$reffa" ]; then echo "Reference genome $reffa does not exist on your filesystem.\n"; exit 9999; fi
if [ ! -f "$pdef" ]; then echo "Prob definition $pdef does not exist on your filesystem.\n"; exit 9999; fi



echo "Validation is Done!";
#exit

for f in ls "$inputDir"/*.fastq.gz
do 
if [ $echo $f != 'ls' ]
then
S=$(basename $f)
variable=$(echo $S | awk -F"_R" '{print $1;}')

#args+=("variable")
var1=$var1' '$variable
#echo $variable
#variable=$variable
fi
done
#echo $var1 
unqpre=$(echo "$var1" | xargs -n1 | sort -u | xargs)
#echo"${args[@]}"
echo $unqpre
echo "total number of Libraries are " ; echo "$unqpre" | wc -w;

flag=0

for f in $unqpre
do
	echo "using prefix $f checking the existance of both pair-end files"
	filep1=$inputDir/"$f"_R1.fastq.gz
	filep2=$inputDir/"$f"_R2.fastq.gz
	if [ -f "$filep1" ]; then
		if [ -f "$filep2" ]; then
		echo "\nFound status of pair-end : OK\n\nStart of Mapping $filep1 using bwa-bwasw"
		./bwa bwasw -t 32 $reffa $filep1 | gzip - > $outputDirSam/$f"_R1.sam.gz"
		echo "Start of Mapping $filep2 using bwa-bwasw\n"
		./bwa bwasw -t 32 $reffa $filep2 | gzip - > $outputDirSam/$f"_R2.sam.gz"
		echo "Computing Nodewalk\n"
		python nodewalk_github_v1.py $f $renzyme $pdef $reffa $inputDir $outputDirSam $outputDirStats
		#python ProcessNWSam_chimeric.py $f $renzyme $pdef $reffa $inputDir $outputDirSam $outputDirStats
		echo "\nProcessing of $f is finished\n"
		flag=$((flag+1))

	fi
	else 
		echo " \n\nfile do not exist"
	fi
	
	#head UPLOADS/NW_ST1/$f"_R1.fastq.gz"
	#printf 'start of second file\n\n\n\n****************************************************************************************************************'
done

echo "$flag Libraries have been processed, finished"


