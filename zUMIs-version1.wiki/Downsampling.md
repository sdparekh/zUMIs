zUMIs has powerful downsampling capabilites. Independent of downsampling mode, the full data is always exported aswell.

- Adaptive downsampling: According to the recommendation of the Scater package (McCarthy et al., 2017) reads are downsampled to be within 3 times median absolute deviation. This is the default setting (eg. -d 0).
- downsampling to a fixed depth: Reads are downsampled to a user-specified depth. Any barcodes that do not reach the requested depth are omitted. Example: -d 10000
- downsampling to a depth range: Barcodes with read depth above the maximum of the range are downsampled to this value. All barcodes within the range are reported without downsampling and barcodes below the minimum specified read depth are ommited. Example: -d 10000-20000
- downsampling to several depths: Several depths can be requested by comma separation. Combinations of fixed depth and depth ranges may be given. Example: -d 10000,10000-20000,30000

```
bash zUMIs-master.sh -f barcoderead.fastq -r cdnaread.fastq -n test -g hg38_STAR5idx_noGTF/ -o ./ -a Homo_sapiens.GRCh38.84.gtf -p 8 -s 0 -d 10000,10000-20000,30000 -c 1-6 -m 7-16 -l 50 -b 384

```

### Output
Downsampled count tables are reported in <StudyName>.dgecounts.rds for each feature type (exons, introns, intron.exon). It is a list of all the downsamplings requested. Each of the downsampling list contains "readcounts" & "umicounts".

These tables can be saved as a Tab delimited text file using the code below.
For example:
```
AllCounts <- readRDS("zUMIs_output/expression/example.dgecounts.rds")

downsamp <- unlist(x = AllCounts$exons$downsampled,recursive = F,use.names = T)
lapply(names(downsamp),function(x) write.table(downsamp[[x]],file=paste("zUMIs_output/expression/",x,".txt",sep=""),sep = "\t",row.names = T,col.names = T))

```

