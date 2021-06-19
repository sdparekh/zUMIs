#!/usr/bin/env python3

import pysam
import argparse

def read_cellBCs (bcfile):
    f = open(bcfile,"r")
    bclines = f.readlines()
    bclist = []
    for x in bclines:
        bclist.append(x.split(',')[0])
    bclist.pop(0)
    return(bclist)

def demultiplex_bam (bamfile, bcdict, outpath, pout, pin, chr):
    inbam = pysam.AlignmentFile(bamfile, 'rb', threads = pin)
    
    if chr == 'zunmapped':
      label = 'zunmapped'
      chr = '*'
    elif chr == 'allreads':
      chr = None
    else:
      label = chr

    bcfilehandles = {}
    for bc in bcdict:
        if chr is None:
          path = outpath + bc + '.demx.bam'
        else: 
          path = outpath + bc + "." + label + '.demx.bam'
        bcfilehandles[bc] = pysam.AlignmentFile(path, "wb", template=inbam, threads = pout)

    for read in inbam.fetch(contig=chr, until_eof=True):
        thisBC = read.get_tag("BC")
        if thisBC in bcfilehandles:
            bcfilehandles[thisBC].write(read)

    inbam.close()
    for outfile in bcfilehandles.values():
        outfile.close()


def main():
    parser = argparse.ArgumentParser(add_help=True)
    parser.add_argument('--bam', type=str, metavar='FILENAME',
                        help='Path to BAM file')
    parser.add_argument('--out', type=str, metavar='DIRNAME',
                        help='Path where demultiplexed files go')
    parser.add_argument('--bc', type=str, metavar='FILENAME',
                        help='Path to barcode dictionary.')
    parser.add_argument('--pout', type=int,
                        help='Number of processes for output bams')
    parser.add_argument('--pin', type=int,
                        help='Number of processes for input bam'),
    parser.add_argument('--chr', type=str,
                        help='Chromosome to use')
    args = parser.parse_args()

    #print("Demultiplexing zUMIs bam file...")
    bc_whitelist = read_cellBCs(args.bc)
    demultiplex_bam(
        bamfile = args.bam,
        bcdict = bc_whitelist,
        outpath = args.out,
        pout = args.pout,
        pin = args.pin,
        chr = args.chr)
    #print("Demultiplexing complete.")

if __name__ == "__main__":
    main()
