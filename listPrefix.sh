#!/bin/bash
# LMU Munich. AG Enard
# getting list of prefix from tmpMerge
# Authors: Swati Parekh, Christoph Ziegenhain, Beate Vieth & Ines Hellmann
# Contact: sparekh@age.mpg.de or christoph.ziegenhain@ki.se

zumisdir=$1
tmpMerge=$2
f=$3
project=$4

if [[ $f =~ \.gz$ ]]; then
	pref=`basename $f .gz`
	ls $tmpMerge$pref* | sed "s|$tmpMerge$pref||" | sed 's/.gz//' > $tmpMerge/$project.listPrefix.txt
else
	pref=`basename $f`
	ls $tmpMerge$pref* | sed "s|$tmpMerge$pref||" > $tmpMerge/$project.listPrefix.txt
fi
