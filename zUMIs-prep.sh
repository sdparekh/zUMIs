#!/bin/bash

f1=$1
f2=$2
sn=$3
o=$4
d=$5

# MAKING THE HEADER
echo '#!/bin/bash' >$o/$sn.prep.sh
echo '#SBATCH -n 1' >>$o/$sn.prep.sh
echo '#SBATCH --error='prep'.%J.err' >>$o/$sn.prep.sh
echo '#SBATCH --output='prep'.%J.out' >>$o/$sn.prep.sh
echo '#SBATCH --workdir='$o >>$o/$sn.prep.sh

echo "srun perl $d/fqcheck.pl $f1 $f2 $sn $o" >>$o/$sn.prep.sh

sbatch $o/$sn.prep.sh > $o/$sn.preparejobid.txt
