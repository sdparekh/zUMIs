#!/bin/bash
# LMU Munich. AG Enard
# Pipeline to run UMI-seq analysis from fastq to read count tables.
# Authors: Swati Parekh, Christoph Ziegenhain, Beate Vieth & Ines Hellmann
# Contact: sparekh@age.mpg.de or christoph.ziegenhain@ki.se
vers=2.5.4
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

## Program path ##
	-d  <zUMIs-dir>   	 : Directory containing zUMIs scripts.  Default: path to this script.

zUMIs version $vers

EOF
}

# Define the default variables #
zumisdir=$(dirname `readlink -f $0`)

while getopts ":y:d:h" options; do #Putting <:> between keys implies that they can not be called without an argument.
  case $options in
  y ) yaml=$OPTARG;;
  d ) zumisdir=$OPTARG;;
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
velo=`grep 'velocyto:' $yaml | awk '{print $2}'`


if grep -q 'samtools_exec:' $yaml
  then
    samtoolsexc=`grep 'samtools_exec' $yaml | awk '{print $2}'`
  else
    samtoolsexc=samtools
    echo "samtools_exec: $samtoolsexc" >> $yaml
fi

if grep -q 'pigz_exec:' $yaml
  then
    pigzexc=`grep 'pigz_exec' $yaml | awk '{print $2}'`
  else
    pigzexc=pigz
    echo "pigz_exec: $pigzexc" >> $yaml
fi

if grep -q 'STAR_exec:' $yaml
  then
    starexc=`grep 'STAR_exec' $yaml | awk '{print $2}'`
  else
    starexc=STAR
    echo "STAR_exec: $starexc" >> $yaml
fi

if grep -q 'Rscript_exec:' $yaml
  then
    Rexc=`grep 'Rscript_exec' $yaml | awk '{print $2}'`
  else
    Rexc=Rscript
    echo "Rscript_exec: $Rexc" >> $yaml
fi

if grep -q 'zUMIs_directory:' $yaml
  then
    sed -i "s|zUMIs_directory:.*|zUMIs_directory: $zumisdir|" $yaml
  else
    echo "zUMIs_directory: $zumisdir" >> $yaml
fi

$Rexc $zumisdir/checkyaml.R $yaml > $project.zUMIs_YAMLerror.log
iserror=`tail $project.zUMIs_YAMLerror.log -n1 | awk '{print $2}'`

if [[ $iserror -eq 1 ]] ; then
    echo "YAML file has an error. Look at the zUMIs_YAMLerror.log or contact developers."
    exit 1
fi


echo -e "\n\n You provided these parameters:
 YAML file:	$yaml
 zUMIs directory:		$zumisdir
 STAR executable		$starexc
 samtools executable		$samtoolsexc
 pigz executable		$pigzexc
 Rscript executable		$Rexc
 RAM limit:   $mem_limit
 zUMIs version $vers \n\n" | tee "$outdir/$project.zUMIs_runlog.txt"
date

#check for executables
sam_exc_check=`which $samtoolsexc`
pigz_exc_check=`which $pigzexc`
r_exc_check=`which $Rexc`
star_exc_check=`which $starexc`

if [[ -z "$sam_exc_check" ]] ||
   [[ -z "$pigz_exc_check" ]] ||
   [[ -z "$r_exc_check" ]] ||
   [[ -z "$star_exc_check" ]]
 then
  echo "One or more of your executables were not found. Please check back."
  exit 1
fi

#create output folders
outdir=`grep 'out_dir' $yaml | awk '{print $2}'`
#[ -d $outdir ] || mkdir $outdir
[ -d $outdir/zUMIs_output/ ] || mkdir $outdir/zUMIs_output/
[ -d $outdir/zUMIs_output/expression ] || mkdir $outdir/zUMIs_output/expression
[ -d $outdir/zUMIs_output/stats ] || mkdir $outdir/zUMIs_output/stats
[ -d $outdir/zUMIs_output/.tmpMerge ] || mkdir $outdir/zUMIs_output/.tmpMerge

if [[ ! -d $outdir ]]; then
  mkdir $outdir
  if [ $? -ne 0 ] ; then
      echo "Please provide a valide output directory path."
      exit 1
  fi
fi

if
[[ "$whichStage" == "Filtering" ]]
then
  echo "Filtering..."

  f=`cut -d' ' -f1 <(echo $fqfiles)` # the first fastq file to determine gzip status
  fullsize=`stat -L --printf="%s" $f`

  tmpMerge=$outdir/zUMIs_output/.tmpMerge/

    if [[ $f =~ \.gz$ ]]; then
      $pigzexc -dc $f | head -n 4000000 | $pigzexc > $tmpMerge/$project.1mio.check.fq.gz
      smallsize=`stat --printf="%s" $tmpMerge/$project.1mio.check.fq.gz`
      rm $tmpMerge/$project.1mio.check.fq.gz
      nreads=`expr $fullsize \* 1000000 / $smallsize`

      for i in $fqfiles;do bash $zumisdir/splitfq.sh $i $pigzexc $num_threads $tmpMerge splitfqgz $project $nreads & done
      wait
      pref=`basename $f .gz`
      l=`ls $tmpMerge$pref* | sed "s|$tmpMerge$pref||" | sed 's/.gz//'`
    else
      cat $f | head -n 4000000 > $tmpMerge/$project.1mio.check.fq
      smallsize=`stat --printf="%s" $tmpMerge/$project.1mio.check.fq`
      rm $tmpMerge/$project.1mio.check.fq
      nreads=`expr $fullsize \* 1000000 / $smallsize`

      for i in $fqfiles;do bash $zumisdir/splitfq.sh $i $pigzexc $num_threads $tmpMerge splitfq $project $nreads & done
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
    date

    #run barcode detection
    $Rexc $zumisdir/zUMIs-BCdetection.R $yaml

    #check if BC correction should be performed!
    BCbinTable=$outdir/zUMIs_output/"$project".BCbinning.txt
    if [[ -f "$BCbinTable" ]]; then
      for x in $l; do
        rawbam="$tmpMerge/$project.$x.raw.tagged.bam"
        fixedbam="$tmpMerge/$project.$x.filtered.tagged.bam"
        mv $fixedbam $rawbam
        perl $zumisdir/correct_BCtag.pl $rawbam $fixedbam $BCbinTable $samtoolsexc &
      done
      wait
    fi

fi

if
[[ "$whichStage" == "Filtering" ]] ||
[[ "$whichStage" == "Mapping" ]]
then
  echo "Mapping..."
    $Rexc $zumisdir/zUMIs-mapping.R $yaml
  date
fi

if
[[ "$whichStage" == "Filtering" ]] ||
[[ "$whichStage" == "Mapping" ]] ||
[[ "$whichStage" == "Counting" ]]
then
  echo "Counting..."
  $Rexc $zumisdir/zUMIs-dge2.R $yaml
  date
  if [[ "$velo" == "yes" ]]; then
    $Rexc $zumisdir/runVelocyto.R $yaml
  fi
fi

if
[[ "$whichStage" == "Filtering" ]] ||
[[ "$whichStage" == "Mapping" ]] ||
[[ "$whichStage" == "Counting" ]] ||
[[ "$whichStage" == "Summarising" ]]
then
  if [[ "$isstats" == "yes" ]]; then
    echo "Descriptive statistics..."
      $Rexc $zumisdir/zUMIs-stats2.R $yaml
  fi
  date
fi

#convenience function
#if grep -q 'find_pattern: ATTGCGCAATG' $yaml &&
#   [[ "$whichStage" == "Filtering" ]]
#  then
#    cp $yaml $outdir/$project.allReads.yaml #copy yaml
#    sed -i "s/project: "$project"/project: "$project"_allReads/" $outdir/$project.allReads.yaml
#    sed -i '/find_pattern/d' $outdir/$project.allReads.yaml
#    sed -i '/UMI(/d' $outdir/$project.allReads.yaml
#    bash $zumisdir/zUMIs-master.sh -y $outdir/$project.allReads.yaml
#fi
