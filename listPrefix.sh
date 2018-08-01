#!/bin/bash
# LMU Munich. AG Enard
# getting list of prefix from tmpMerge
# Authors: Swati Parekh, Christoph Ziegenhain, Beate Vieth & Ines Hellmann
# Contact: sparekh@age.mpg.de or christoph.ziegenhain@ki.se

zumisdir=$1
tmpMerge=$2
f=$3
project=$4
yaml=$5
samtoolsexc=$6
Rexc=$7
pigzexc=$8

if [[ $f =~ \.gz$ ]]; then
	pref=`basename $f .gz`
	ls $tmpMerge$pref* | sed "s|$tmpMerge$pref||" | sed 's/.gz//' > $tmpMerge/$project.listPrefix.txt
else
	pref=`basename $f`
	ls $tmpMerge$pref* | sed "s|$tmpMerge$pref||" > $tmpMerge/$project.listPrefix.txt
fi

for x in `cat $tmpMerge/$project.listPrefix.txt`; do perl $zumisdir/fqfilter_v2.pl $yaml $samtoolsexc $Rexc $pigzexc $zumisdir $x & done
wait
