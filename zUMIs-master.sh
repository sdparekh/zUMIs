#!/bin/bash
# LMU Munich. AG Enard
# Pipeline to run UMI-seq analysis from fastq to read count tables.
# Authors: Swati Parekh, Christoph Ziegenhain, Beate Vieth & Ines Hellmann
# Contact: sparekh@age.mpg.de or christoph.ziegenhain@ki.se
vers=0.2.0
function check_opts() {
    value=$1
    name=$2
    flag=$3

    if [[ -z "$value" ]]
    then failure "No $name!! One can not run this pipeline without $flag option."
    fi
}

function failure() {
	echo -e "\n There seems to be a problem. Please check the usage: \n $1 \n\n"
	usage
	exit 1
}

zumis=$0

function usage () {
    cat >&2 <<EOF

  USAGE: $zumis [options]
	-h  Print the usage info.

## Required parameters ##

	-y  <YAML config file> : Path to the YAML config file. Required.

## Program paths ##
	-s  <STAR-executable>	 : path to STAR executable in your system. Default: STAR
	-t  <samtools-executable>: path to samtools executable in your system. Default: samtools
  	-p <pigz-executable> 	 : path to pigz executable in your system. Default: pigz
  	-r <Rscript-executable>  : path to Rscript executable in your system. Default: Rscript
	-d  <zUMIs-dir>   	 : Directory containing zUMIs scripts.  Default: path to this script.

zUMIs version $vers

EOF
}

# Define the default variables #
starexc=STAR
samtoolsexc=samtools
pigzexc=pigz
Rexc=Rscript
zumisdir=$(dirname `readlink -f $0`)

while getopts ":y:s:t:p:r:d:h" options; do #Putting <:> between keys implies that they can not be called without an argument.
  case $options in
  y ) yaml=$OPTARG;;
  s ) starexc=$OPTARG;;
  t ) samtoolsexc=$OPTARG;;
  p ) pigzexc=$OPTARG;;
  d ) zumisdir=$OPTARG;;
  r ) Rexc=$OPTARG;;
  h ) usage
          exit 1;;
  \? ) echo -e "\n This key is not available! Please check the usage again: -$OPTARG"
  	usage
  	exit 1;;
  esac
done

if [[ $OPTIND -eq 1 ]] ; then
    usage
    exit 1
fi

check_opts "$yaml" "YAML" "-y"

#now get some variables from YAML
num_threads=`grep 'num_threads' $yaml | awk '{print $2}'`
project=`grep 'project:' $yaml | awk '{print $2}'`
whichStage=`grep 'which_Stage:' $yaml | awk '{print $2}'`
outdir=`grep 'out_dir' $yaml | awk '{print $2}'`
isslurm=`grep 'use_SLURM:' $yaml | awk '{print $2}'`
genomedir=`grep 'STAR_index:' $yaml | awk '{print $2}'`
mem_limit=`grep 'mem_limit:' $yaml | awk '{print $2}'`
isstats=`grep 'make_stats:' $yaml | awk '{print $2}'`

echo -e "\n\n You provided these parameters:
 YAML file:	$yaml
 zUMIs directory:		$zumisdir
 STAR executable		$starexc
 samtools executable		$samtoolsexc
 pigz executable		$pigzexc
 Rscript executable		$Rexc
 RAM limit:   $mem_limit
 zUMIs version $vers \n\n" | tee "$outdir/zUMIs_runlog.txt"

#create output folders
outdir=`grep 'out_dir' $yaml | awk '{print $2}'`
[ -d $outdir/zUMIs_output/ ] || mkdir $outdir/zUMIs_output/
[ -d $outdir/zUMIs_output/expression ] || mkdir $outdir/zUMIs_output/expression
[ -d $outdir/zUMIs_output/stats ] || mkdir $outdir/zUMIs_output/stats


if
[[ "$whichStage" == "Filtering" ]]
then
  echo "Starting Filtering..."
  echo "samtools_exec: $samtoolsexc" >> $yaml
  echo "pigz_exec: $pigzexc" >> $yaml
  echo "STAR_exec: $starexc" >> $yaml
  echo "zUMIs_directory: $zumisdir" >> $yaml
  if [[ "$isslurm" == "yes" ]]; then
    sbatch --cpus-per-task=$num_threads --wrap=="perl $zumisdir/fqfilter_v2.pl $yaml $samtoolsexc $Rexc $pigzexc $zumisdir" > $outdir/$project.slurmjobid.txt
  else
    perl $zumisdir/fqfilter_v2.pl $yaml $samtoolsexc $Rexc $pigzexc $zumisdir
  fi
fi

if
[[ "$whichStage" == "Filtering" ]] ||
[[ "$whichStage" == "Mapping" ]]
then
  echo "Starting Mapping..."
  if [[ "$isslurm" == "yes" ]]; then
    memory=`du -sh $genomedir | cut -f1` #STAR genome index size
    j=`cat $outdir/$project.slurmjobid.txt | cut -f4 -d' '`
    sbatch --dependency=afterok:'$j' --mem=$memory --cpus-per-task=$num_threads --wrap=="$Rexc $zumisdir/zUMIs-mapping.R $yaml" > $outdir/$project.slurmjobid.txt
  else
    $Rexc $zumisdir/zUMIs-mapping.R $yaml
  fi
fi

if
[[ "$whichStage" == "Filtering" ]] ||
[[ "$whichStage" == "Mapping" ]] ||
[[ "$whichStage" == "Counting" ]]
then
  echo "Starting Counting..."
  yamlnew=$outdir/$project.postmap.yaml  #note! mapping creates a temporary new yaml because of custom GTF!
  if [[ "$isslurm" == "yes" ]]; then
    if [[ $mem_limit == "null" ]]; then
      mem_limit=`du -sh $genomedir | cut -f1`
    else
      mem_limit=`expr $mem_limit \* 1000`
    fi
    j=`cat $outdir/$project.slurmjobid.txt | cut -f4 -d' '`
    sbatch --dependency=afterok:'$j' --mem=$mem_limit --cpus-per-task=$num_threads --wrap=="$Rexc $zumisdir/zUMIs-dge2.R $yamlnew" > $outdir/$project.slurmjobid.txt
   else
     $Rexc $zumisdir/zUMIs-dge2.R $yamlnew
   fi
fi

if
[[ "$whichStage" == "Filtering" ]] ||
[[ "$whichStage" == "Mapping" ]] ||
[[ "$whichStage" == "Counting" ]] ||
[[ "$whichStage" == "Summarising" ]]
then
  yamlnew=$outdir/$project.postmap.yaml  #note! mapping creates a temporary new yaml because of custom GTF!
  if [[ "$isstats" == "yes" ]]; then
    echo "Starting descriptive statistics..."
    if [[ "$isslurm" == "yes" ]]; then
      j=`cat $outdir/$project.slurmjobid.txt | cut -f4 -d' '`
      sbatch --dependency=afterok:'$j' --mem=5000 --cpus-per-task=$num_threads --wrap=="$Rexc $zumisdir/zUMIs-stats2.R $yamlnew" > $outdir/$project.slurmjobid.txt
    else
      $Rexc $zumisdir/zUMIs-stats2.R $yamlnew
    fi
  fi
fi
