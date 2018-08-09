zUMIs is also compatible with combinatorial indexing protocols, such as sci-RNA-seq (Cao et al., 2017) and SPLiT-seq (Rosenberg et al., 2017).

However, because of their structure, these protocols need a preprocessing step before they may be used in zUMIs.

We provide perl scripts for preprocessing this type of data.

### SPLiT-seq
SPLiT-seq contains cell barcodes that are ligated after each split/pool step during the library preparation.
Thus, the final barcode read contains the actual barcode bases interspersed by fixed ligation linkers that should be removed prior to invoking zUMIs.

Use the `preprocess_splitseq.pl` script provided and use as follows.

Example:
`preprocess_splitseq.pl read2.fq.gz 1-18 49-56 87-94 16 read2 /your/output/dir pigz`

The input arguments are defined as follows:
- Read2 fasta file
- Range of UMI sequence + first barcode segment (Round 1 RT barcode)
- Range of second barcode segment (Round 2 Ligation Barcodes)
- Range of third barcode segment (Round 3 Ligation Barcodes)
- Number of threads to use for zipping the output fastq file
- Output file prefix, `.barcoderead.preprocess.fastq.gz` will be added
- Output directory path
- path to the pigz executable

This example call will generate `read2.barcoderead.preprocess.fastq.gz` as output.
Input that into zUMIs with `-f read2.barcoderead.preprocess.fastq.gz` and define the barcode ranges depending on the ranges selected in the preprocessing: eg. `-m 1-10 -c 11-34`


### sci-RNA-seq
sci-RNA-seq contains three levels of barcodes: RT barcodes and Illumina i7 and i5 indices.
Thus, the barcode sequence is split over several reads/fastq files that should be combined prior to invoking zUMIs.

Use the `cat3fq.pl` script provided and use as follows.

Example:
`cat3fq.pl R1.fastq.gz I1.fastq.gz I2.fastq.gz Combined_barcodeRead.fq 16`
The input arguments are defined as follows:
- Read1 fasta file (contains RT barcode)
- Index 1 fastq file (contains i7 barcode)
- Index 2 fastq file (contains i5 barcode)
- Output file prefix, `.gz` will be added
- Number of threads to use for zipping the output fastq file

This example call will generate `Combined_barcodeRead.fq.gz` as output.
Input that into zUMIs with `-f Combined_barcodeRead.fq.gz` and define the barcode ranges depending on the length of the input files selected in the preprocessing: eg. `-m 1-8 -c 9-48`
