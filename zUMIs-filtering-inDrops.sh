#!/bin/bash

f1=$1
f2=$2
f3=$3
f4=$4
sn=$5
o=$6
xm=$7
cbq=$8
mbq="${9}"
mq="${10}"
cq="${11}"
t="${12}"
d="${13}"

# MAKING THE HEADER
echo '#!/bin/bash' >$o/$sn.prep.sh
echo '#SBATCH -n 1' >>$o/$sn.prep.sh
echo '#SBATCH --error='prep'.%J.err' >>$o/$sn.prep.sh
echo '#SBATCH --output='prep'.%J.out' >>$o/$sn.prep.sh
echo '#SBATCH --cpus-per-task='$t >>$o/$sn.prep.sh
echo '#SBATCH --workdir='$o >>$o/$sn.prep.sh

echo "srun perl $d/fqfilter-inDrops.pl $f1 $f2 $f3 $f4 $cq $cbq $mq $mbq $xm $t $sn $o" >>$o/$sn.prep.sh

sbatch $o/$sn.prep.sh > $o/$sn.preparejobid.txt
