#!/bin/bash

f1=$1
f2=$2
sn=$3
o=$4
xc=$5
xm=$6
cbq=$7
mbq=$8
mq=$9
cq="${10}"
t="${11}"
d="${12}"

# MAKING THE HEADER
echo '#!/bin/bash' >$o/$sn.prep.sh
echo '#SBATCH -n 1' >>$o/$sn.prep.sh
echo '#SBATCH --error='prep'.%J.err' >>$o/$sn.prep.sh
echo '#SBATCH --output='prep'.%J.out' >>$o/$sn.prep.sh
echo '#SBATCH --cpus-per-task='$t >>$o/$sn.prep.sh
echo '#SBATCH --workdir='$o >>$o/$sn.prep.sh

echo "srun perl $d/fqfilter.pl $f1 $f2 $cq $cbq $mq $mbq $xc $xm $t $sn $o" >>$o/$sn.prep.sh

sbatch $o/$sn.prep.sh > $o/$sn.preparejobid.txt
