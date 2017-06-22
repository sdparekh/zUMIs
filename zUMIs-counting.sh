#!/bin/bash
sn=$1
o=$2
n=$3
t=$4
g=$5
s=$6
xc=$7
xm=$8
subs=$9
d="${10}"
stats="${11}"
whichStage="${12}"
isstrt="${13}"
f3="${14}"
xc2="${15}"

j=`cat $o/$sn.mapjobid.txt | cut -f4 -d' '`
u=`cat $o/$sn.unmapjobid.txt | cut -f4 -d' '`

if [[ "$isstrt" == "no" ]] ; then
	xcst=`echo $xc | cut -f1 -d '-'`
	xcend=`echo $xc | cut -f2 -d '-'`
	xmst=`echo $xm | cut -f1 -d '-'`
	xmend=`echo $xm | cut -f2 -d '-'`
else
	xcst=1
	a=`echo $xc2 | cut -f2 -d '-'`
	b=`echo $xc | cut -f2 -d '-'`
	xcend=`expr $a + $b`
	xmst=`expr $xcend + 1`
	c=`echo $xm | cut -f2 -d '-'`
	xmend=`expr $c + $xmst - 1`
fi

re='^[0-9]+$'

# MAKING THE HEADER
echo '#!/bin/bash' >$o/$sn.dge.sh
echo '#SBATCH -n 1' >>$o/$sn.dge.sh
echo '#SBATCH --error='dge'.%J.err' >>$o/$sn.dge.sh
echo '#SBATCH --output='dge'.%J.out' >>$o/$sn.dge.sh
echo '#SBATCH --workdir='$o >>$o/$sn.dge.sh
echo '#SBATCH --dependency=afterok:'$j':'$u >>$o/$sn.dge.sh
echo '#SBATCH --mem=50000' >>$o/$sn.dge.sh
echo '#SBATCH --cpus-per-task='$t >>$o/$sn.dge.sh

echo "ln -s -f $o/$sn.aligned.sorted.bam $o/$sn.aligned.sorted.bam.in" >>$o/$sn.dge.sh
echo "ln -s -f $o/$sn.aligned.sorted.bam $o/$sn.aligned.sorted.bam.ex" >>$o/$sn.dge.sh

if [[ "$whichStage" != "summarising" ]] ; then
	if [[ $n =~ $re ]] ; then
		echo "srun --chdir=$o Rscript $d/zUMIs-dge.R --gtf $g --abam $o/$sn.aligned.sorted.bam --ubam $o/$sn.barcodelist.filtered.sort.sam --barcodenumber $n --out $o --sn $sn --cores $t --strandedness $s --bcstart $xcst --bcend $xcend --umistart $xmst --umiend $xmend --subsamp $subs" >>$o/$sn.dge.sh
	else
		echo "srun --chdir=$o Rscript $d/zUMIs-dge.R --gtf $g --abam $o/$sn.aligned.sorted.bam --ubam $o/$sn.barcodelist.filtered.sort.sam --barcodefile $n --out $o --sn $sn --cores $t --strandedness $s --bcstart $xcst --bcend $xcend --umistart $xmst --umiend $xmend --subsamp $subs" >>$o/$sn.dge.sh
	fi

	if [[ $stats == "yes" ]] ; then
		echo "srun --chdir=$o Rscript $d/zUMIs-stats.R --out $o --sn $sn  --bcstart $xcst --bcend $xcend --umistart $xmst --umiend $xmend" >>$o/$sn.dge.sh
	fi
else
	if [[ $stats == "yes" ]] ; then
		echo "srun --chdir=$o Rscript $d/zUMIs-stats.R --out $o --sn $sn  --bcstart $xcst --bcend $xcend --umistart $xmst --umiend $xmend" >>$o/$sn.dge.sh
	else
		echo "You need to switch on -S <isStats> option."
	fi
fi


sbatch $o/$sn.dge.sh > $o/$sn.dgejobids.txt
