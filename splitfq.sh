#!/bin/bash
# LMU Munich. AG Enard
# splitting fastq files to filter using multiple processors.
# Authors: Swati Parekh, Christoph Ziegenhain, Beate Vieth & Ines Hellmann
# Contact: sparekh@age.mpg.de or christoph.ziegenhain@ki.se

function splitfq() {
	fqfile=$1
  pexc=$2
  nthreads=$3
  nreads=$4
  t=$5
	project=$6

  n=`expr $nreads / $nthreads`
	n=`expr $n + 1`
	nl=`expr $n \* 4`
  pref=`basename $fqfile`
  d=`dirname $fqfile`
  split --lines=$nl --filter=''$pexc' -p '$nthreads' > $FILE.gz' $fqfile $t$pref$project

	ls $t$pref$project* | sed "s|$t$pref$project||" > $t/$project.listPrefix.txt

	exit 1
}
function splitfqgz() {
	fqfile=$1
  pexc=$2
  nthreads=$3
  nreads=$4
  t=$5
	project=$6

  n=`expr $nreads / $nthreads`
  n=`expr $n + 1`
  nl=`expr $n \* 4`
  pref=`basename $fqfile .gz`
  d=`dirname $fqfile`
  $pexc -dc -p $nthreads $d/$pref.gz | split --lines=$nl --filter=''$pexc' -p '$nthreads' > $FILE.gz' - $t$pref$project

	ls $t$pref$project* | sed "s|$t$pref$project||" | sed 's/.gz//' > $t/$project.listPrefix.txt

	exit 1
}

i=$1
pigzexc=$2
num_threads=$3
tmpMerge=$4
fun=$5
project=$6
nreads=$7



$fun "$i" "$pigzexc" "$num_threads" "$nreads" "$tmpMerge" "$project"
