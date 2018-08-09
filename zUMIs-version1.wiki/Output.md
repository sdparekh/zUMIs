zUMIs' output is structured in three subdirectories:

```
zUMIs_output/filtered_fastq
zUMIs_output/expression
zUMIs_output/stats
```

- "filtered_fastq" contains the filtered, gzipped fastq files for cDNA and barcode reads.
- "expression" contains .rds files
   - generated reference annotation(*.annotationSAF.rds)
   - list of count matrices as sparseMatrix (*.dgecounts.rds)
   - intermediate tbl with the data in long format (*.tbl.rds)
- "stats" contains plots and data files with descriptive statistics
- STAR output files and feautreCounts per reads files are stored in the parent directory defined by the user (-o)
- If the jobs were run using slurm, every step creates .err and .out files with the prefix of stage name.
   e.g. map.???.err and map.???.out files contain error messages and standard output from mapping, respectively.

## Structure of the output dgecounts object in <sn>.dgecounts.rds 

zUMIs produces dge output in .rds format that can be read in R with the following command.

```
AllCounts <- readRDS("zUMIs_output/expression/example.dgecounts.rds")
names(AllCounts)
[1] "introns"     "exons"       "intron.exon"

names(AllCounts$exons)
[1] "readcounts"  "umicounts"   "downsampled"

names(AllCounts$exons$downsampled)
[1] "downsampled_7358"

```
AllCounts is a list of lists with all the count matrices as sparseMatrix. The parent list is three feature types (introns,exons and intron+exon) and each of them contain three subtypes with "readcounts" -- (without removing duplicates), "umicounts" -- (removed duplicates) and "downsampled" -- a list of all the downsampling sizes requested. Each of the downsampling list also contains "readcounts" & "umicounts".

The sparseMatrix can be converted into a conventional count table using "as.matrix" function in R and saved as a text file using the code below.
```
#Feature exons umicounts table
dge <- as.matrix(AllCounts$exons$umicounts)
write.table(dge,"exons.umicounts.txt",quote=F,sep="\t")
```
To convert all the count tables recursively and save them as text files. 
```
a <- lapply(AllCounts,function(x) {
  tmp <- unlist(x)
  lapply(tmp,as.matrix)
  })

#To save each count table as a text file
lapply(names(a),function(x) 
  lapply(names(a[[x]]),function(xx) 
    write.table(a[[x]][[xx]],file=paste(x,xx,"txt",sep = "."),quote = F,sep="\t")))

```
NOTE: It is highly recommended to retain sparseMatrix format for the downstream analysis, especially when you have >20,000 cells to save time and space.