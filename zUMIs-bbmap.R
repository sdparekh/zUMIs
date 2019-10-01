#!/usr/bin/env Rscript
suppressMessages(require(yaml))

y <- commandArgs(trailingOnly = T)

inp<-yaml::read_yaml(y)
samtools <- inp$samtools_exec
pigz <- inp$pigz_exec
cores <- inp$num_threads
mem <- inp$mem_limit

###### setup files
outbam <- paste0(inp$out_dir,"/",inp$project,".filtered.tagged.Aligned.out.bam")
starbam <- paste0(inp$out_dir,"/",inp$project,".filtered.tagged.Aligned.STAR.bam")
bbbam <- paste0(inp$out_dir,"/",inp$project,".filtered.tagged.Aligned.bbmap.sam.gz")
fastqfile <- paste0(inp$out_dir,"/",inp$project,".filtered.tagged.cDNA.fastq.gz")
#fasta <- paste0(inp$reference$STAR_index,"/Genome")
#fasta <- "/data/ngs/genomes/Human/hg38/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
fasta <- inp$reference$bbmap_fasta

###### make Fastq
system(paste("mv",outbam,starbam))
#system(paste0("picard SamToFastq I=",starbam," O=",fastqfile))
#
###### run BBmap
bb_cmd <- paste0(samtools," fastq ",starbam,
                 " | bbmap local=t maxindel=200000 nhtag=t nmtag=t ambiguous=best ordered=t int=f out=stdout.sam in=stdin.fq",
                " out=",bbbam,
                " threads=",cores,
                " ref=",fasta,
                " -Xmx",mem,"g ",
                "nodisk"
              )
system(bb_cmd)
#~/programs/bbmap/bbmap.sh in=../run_test/Example.filtered.cDNA.fastq.gz unpigz=t pigz=t threads=32 local=t ambiguous=best ordered=t out=Example.filtered.bbmap.sam.gz -Xmx50g ref=hg38.primary_assembly.fa nodisk


###### combine outputs
merge_cmd <- paste("perl",paste0(inp$zUMIs_directory,"/merge_bbmap_alignment.pl"),starbam,bbbam,outbam,samtools,pigz)
system(merge_cmd)
#perl merge_bbmap_alignment.pl Example.filtered.tagged.Aligned.STAR.bam Example.filtered.bbmap.sam.gz Example.filtered.tagged.Aligned.out.bam samtools pigz
