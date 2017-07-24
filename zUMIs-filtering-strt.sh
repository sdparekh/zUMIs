#!/bin/bash

f1=$1
f2=$2
sn=$3
o=$4
xm=$5
cbq=$6
mbq=$7
mq=$8
cq="${9}"
t="${10}"
d="${11}"
f3="${12}"
bt="${13}"

# MAKING THE HEADER
echo '#!/bin/bash' >$o/$sn.prep.sh
echo '#SBATCH -n 1' >>$o/$sn.prep.sh
echo '#SBATCH --error='prep'.%J.err' >>$o/$sn.prep.sh
echo '#SBATCH --output='prep'.%J.out' >>$o/$sn.prep.sh
echo '#SBATCH --cpus-per-task='$t >>$o/$sn.prep.sh
echo '#SBATCH --workdir='$o >>$o/$sn.prep.sh

echo "srun perl $d/fqfilter-strt.pl $f1 $f2 $f3 $cq $cbq $mq $mbq $xm $bt $t $sn $o" >>$o/$sn.prep.sh

sbatch $o/$sn.prep.sh > $o/$sn.preparejobid.txt
