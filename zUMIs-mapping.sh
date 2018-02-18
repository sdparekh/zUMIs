#!/bin/bash

sn=$1
o=$2
g=$3
gtf=$4
t=$5
r=$6
j=`cat $o/$sn.preparejobid.txt | cut -f4 -d' '`
m=`du -B 1000000 -s $g | cut -f1`
x=$7
starexc=$8
samtoolsexc=$9
xm="${10}"
bt="${11}"
isstrt="${12}"

if [[ "$isstrt" == "no" ]] ; then
	rl=`expr $r - 1`
else
	c=`echo $xm | cut -f2 -d '-'`
	rl=`expr $c + $bt - 1`
fi

# MAKING THE HEADER
echo '#!/bin/bash' >$o/$sn.map.sh  
echo '#SBATCH -n 1' >>$o/$sn.map.sh 
echo '#SBATCH --error='map'.%J.err' >>$o/$sn.map.sh 
echo '#SBATCH --output='map'.%J.out' >>$o/$sn.map.sh
echo '#SBATCH --cpus-per-task='$t >>$o/$sn.map.sh
echo '#SBATCH --workdir='$o >>$o/$sn.map.sh
echo '#SBATCH --dependency=afterok:'$j >>$o/$sn.map.sh
echo '#SBATCH --mem='$m >>$o/$sn.map.sh

echo "srun $starexc --genomeDir $g --runThreadN $t --readFilesCommand zcat --sjdbGTFfile $gtf --outFileNamePrefix $o/$sn. --outSAMtype BAM Unsorted --outSAMmultNmax 1 --outFilterMultimapNmax 50 --outSAMunmapped Within --sjdbOverhang $rl --twopassMode Basic --readFilesIn $o/$sn.cdnaread.filtered.fastq.gz $x" >>$o/$sn.map.sh
echo "srun $samtoolsexc sort -n -O bam -T $o/temp.$sn -@ $t -m 2G -o $o/$sn.aligned.sorted.bam $o/$sn.Aligned.out.bam" >>$o/$sn.map.sh

sbatch $o/$sn.map.sh > $o/$sn.mapjobid.txt
