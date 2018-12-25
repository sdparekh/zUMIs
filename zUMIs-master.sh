#!/bin/bash
# LMU Munich. AG Enard
# Pipeline to run UMI-seq analysis from fastq to read count tables.
# Authors: Swati Parekh, Christoph Ziegenhain, Beate Vieth & Ines Hellmann
# Contact: sparekh@age.mpg.de or christoph.ziegenhain@ki.se
vers=2.2.2
currentv=`curl -s https://raw.githubusercontent.com/sdparekh/zUMIs/master/zUMIs-master.sh | grep '^vers=' | cut -f2 -d "="`
if [ "$currentv" != "$vers" ]; then echo -e "------------- \n\n Good news! A newer version of zUMIs is available at https://github.com/sdparekh/zUMIs \n\n-------------"; fi

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
genomedir=`grep 'STAR_index:' $yaml | awk '{print $2}'`
mem_limit=`grep 'mem_limit:' $yaml | awk '{print $2}'`
isstats=`grep 'make_stats:' $yaml | awk '{print $2}'`
fqfiles=`grep 'name:' $yaml | awk '{print $2}'`

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
[ -d $outdir/zUMIs_output/.tmpMerge ] || mkdir $outdir/zUMIs_output/.tmpMerge

if grep -q 'samtools_exec:' $yaml
  then
    sed -i "s|samtools_exec:.*|samtools_exec: $samtoolsexc|" $yaml
  else
    echo "samtools_exec: $samtoolsexc" >> $yaml
  fi

  if grep -q 'pigz_exec:' $yaml
  then
    sed -i "s|pigz_exec:.*|pigz_exec: $pigzexc|" $yaml
  else
    echo "pigz_exec: $pigzexc" >> $yaml
  fi

  if grep -q 'STAR_exec:' $yaml
  then
    sed -i "s|STAR_exec:.*|STAR_exec: $starexc|" $yaml
  else
    echo "STAR_exec: $starexc" >> $yaml
  fi

  if grep -q 'zUMIs_directory:' $yaml
  then
    sed -i "s|zUMIs_directory:.*|zUMIs_directory: $zumisdir|" $yaml
  else
    echo "zUMIs_directory: $zumisdir" >> $yaml
  fi


if
[[ "$whichStage" == "Filtering" ]]
then
  echo "Filtering..."

  f=`cut -d' ' -f1 <(echo $fqfiles)` # the first fastq file to determine gzip status

  tmpMerge=$outdir/zUMIs_output/.tmpMerge/

    if [[ $f =~ \.gz$ ]]; then
      for i in $fqfiles;do bash $zumisdir/splitfq.sh $i $pigzexc $num_threads $tmpMerge splitfqgz $project $f & done
      wait
      pref=`basename $f .gz`
      l=`ls $tmpMerge$pref* | sed "s|$tmpMerge$pref||" | sed 's/.gz//'`
    else
      for i in $fqfiles;do bash $zumisdir/splitfq.sh $i $pigzexc $num_threads $tmpMerge splitfq $project $f & done
      wait
      pref=`basename $f`
      l=`ls $tmpMerge$pref* | sed "s|$tmpMerge$pref||"`
    fi

    for x in $l; do perl $zumisdir/fqfilter_v2.pl $yaml $samtoolsexc $Rexc $pigzexc $zumisdir $x & done
    wait
    bash $zumisdir/mergeBAM.sh $zumisdir $tmpMerge $num_threads $project $outdir $yaml $samtoolsexc
    for i in $fqfiles; do
      pref=`basename $i | sed 's/.fastq.gz//' | sed 's/.fq.gz//'`
      rm $tmpMerge$pref*gz
    done
fi

if
[[ "$whichStage" == "Filtering" ]] ||
[[ "$whichStage" == "Mapping" ]]
then
  echo "Mapping..."
    $Rexc $zumisdir/zUMIs-mapping.R $yaml
fi

if
[[ "$whichStage" == "Filtering" ]] ||
[[ "$whichStage" == "Mapping" ]] ||
[[ "$whichStage" == "Counting" ]]
then
  echo "Counting..."
  yamlnew=$outdir/$project.postmap.yaml  #note! mapping creates a temporary new yaml because of custom GTF!
  $Rexc $zumisdir/zUMIs-dge2.R $yamlnew
fi

if
[[ "$whichStage" == "Filtering" ]] ||
[[ "$whichStage" == "Mapping" ]] ||
[[ "$whichStage" == "Counting" ]] ||
[[ "$whichStage" == "Summarising" ]]
then
  yamlnew=$outdir/$project.postmap.yaml  #note! mapping creates a temporary new yaml because of custom GTF!
  if [[ "$isstats" == "yes" ]]; then
    echo "Descriptive statistics..."
      $Rexc $zumisdir/zUMIs-stats2.R $yamlnew
  fi
fi
