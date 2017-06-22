#!/bin/bash

sn=$1
o=$2
j=`cat $o/$sn.dgejobids.txt | cut -f4 -d' '`

# MAKING THE HEADER
echo '#!/bin/bash' >$o/$sn.clean.sh
echo '#SBATCH -n 1' >>$o/$sn.clean.sh
echo '#SBATCH --error='clean'.%J.err' >>$o/$sn.clean.sh
echo '#SBATCH --output='clean'.%J.out' >>$o/$sn.clean.sh
echo '#SBATCH --workdir='$o >>$o/$sn.clean.sh
echo '#SBATCH --dependency=afterok:'$j >>$o/$sn.clean.sh

echo "srun rm $o/$sn.Aligned.out.bam $o/$sn.aligned.sorted.bam.in $o/$sn.aligned.sorted.bam.ex $o/$sn.barcodelist.filtered.sam" >>$o/$sn.clean.sh
echo "srun mv $o/$sn.barcoderead.filtered.fastq.gz $o/zUMIs_output/filtered_fastq/" >>$o/$sn.clean.sh
echo "srun mv $o/$sn.cdnaread.filtered.fastq.gz $o/zUMIs_output/filtered_fastq/" >>$o/$sn.clean.sh

sbatch $o/$sn.clean.sh
