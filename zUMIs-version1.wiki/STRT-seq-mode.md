STRT-seq data can be processed by zUMIs by switching on the STRT-seq mode with the -y flag. 
Trimming of the TSO-derived G homopolymers is handled within the filtering step (eg -j 3).
Please set the barcode and UMI base-ranges from 1 (eg. -m 1-6).

One needs to carefully set the following keys in STRT-seq mode.

### STRT-seq mode parameters setting
-	-c  <XC baserange>       : Base range for cell/sample barcode in -f Barcode read(e.g. 1-6).  Required.
                              For STRT-seq give this as 1-n where n is your first cell barcode(-f) length.
-	-m  <XM baserange>       : Base range for UMI barcode in -f Barcode read(e.g. 7-16).  Required.
                              For STRT-seq give this as 1-n where n is your UMI barcode length.
-	-l  <readlength>         : Read length of -r cDNA reads (e.g. 50).  Required.
                              For STRT-seq give this as a total length of your umicdna reads.
-	-y  <STRT-seq>		 : Do you have STRT-seq data? yes/no Default: no.
-	-F  <BarcodeRead2 fastq> : In case dual index is used, provide the second cell barcode index read <-F> here. Default: NA.
-	-C  <XC2 baserange> 	 : Base range for cell/sample barcode in -F Barcode read(e.g. 1-5).  Required.
-	-j  <BaseTrim>		 : <-j INT> fixed number of bases(G) will be trimmed between UMI and cDNA read for STRT-seq. Default: 3.

Here is a STRT-2i example:

```
bash zUMIs-master.sh -f barcode1read.fastq -F barcode2read.fastq -r umicdnaread.fastq -n test2i -g hg38_STAR5idx_noGTF/ -o ./ -a Homo_sapiens.GRCh38.84.gtf -p 8 -c 1-8 -C 1-5 -m 1-6 -l 51 -b 9600 -q 17 -Q 17 -y yes -j 3

```
