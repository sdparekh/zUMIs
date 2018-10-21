#!/bin/bash
# LMU Munich. AG Enard
# getting list of prefix from tmpMerge
# Authors: Swati Parekh, Christoph Ziegenhain, Beate Vieth & Ines Hellmann
# Contact: sparekh@age.mpg.de or christoph.ziegenhain@ki.se

zumisdir=$1
tmpMerge=$2
nthreads=$3
project=$4
outdir=$5
yaml=$6
samtoolsexc=$7

ls $tmpMerge/$project.*.filtered.tagged.bam > $tmpMerge/$project.bamlist.txt
#samtools merge -f -@ $nthreads -b $tmpMerge/$project.bamlist.txt $outdir/$project.filtered.tagged.bam > /dev/null 2>&1

cat $tmpMerge/$project.*.BCstats.txt > $outdir/$project.BCstats.txt

f=`head -n1 $tmpMerge/$project.bamlist.txt`
flag=`$samtoolsexc view $f | head -n1 | cut -f2`

if [[ $flag == 4 ]]; then
	if grep -q 'read_layout:' $yaml
	then
		sed -i "s|read_layout:.*|read_layout: SE|" $yaml
	else
		echo "read_layout: SE" >> $yaml
	fi
else
	if grep -q 'read_layout:' $yaml
	then
		sed -i "s|read_layout:.*|read_layout: PE|" $yaml
	else
		echo "read_layout: PE" >> $yaml
	fi
fi
