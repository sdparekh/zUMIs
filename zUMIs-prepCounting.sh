#!/bin/bash

sn=$1
o=$2
t=$3
m=`expr $t \* 2`
j=`cat $o/$sn.preparejobid.txt | cut -f4 -d' '`
samtoolsexc=$4

# MAKING THE HEADER
echo '#!/bin/bash' >$o/$sn.unmap.sh  
echo '#SBATCH -n 1' >>$o/$sn.unmap.sh 
echo '#SBATCH --error='unmap'.%J.err' >>$o/$sn.unmap.sh 
echo '#SBATCH --output='unmap'.%J.out' >>$o/$sn.unmap.sh
echo '#SBATCH --cpus-per-task='$t >>$o/$sn.unmap.sh
echo '#SBATCH --workdir='$o >>$o/$sn.unmap.sh
echo '#SBATCH --dependency=afterok:'$j >>$o/$sn.unmap.sh
echo '#SBATCH --mem='$m >>$o/$sn.unmap.sh

echo "srun $samtoolsexc sort -n -O sam -T $o/tmp.$sn -@ $t -m 2G -o $o/$sn.barcodelist.filtered.sort.sam $o/$sn.barcodelist.filtered.sam" >>$o/$sn.unmap.sh

sbatch $o/$sn.unmap.sh > $o/$sn.unmapjobid.txt
