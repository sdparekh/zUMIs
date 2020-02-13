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

def demultiplex_bam (bamfile, bcdict, outpath, pout, pin):
    inbam = pysam.AlignmentFile(bamfile, 'rb', threads = pin)

    bcfilehandles = {}
    for bc in bcdict:
        path = outpath + bc + '.demx.bam'
        bcfilehandles[bc] = pysam.AlignmentFile(path, "wb", template=inbam, threads = pout)

    for read in inbam:
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
                        help='Number of processes for input bam')
    args = parser.parse_args()

    print("Demultiplexing zUMIs bam file...")
    bc_whitelist = read_cellBCs(args.bc)
    demultiplex_bam(
        bamfile = args.bam,
        bcdict = bc_whitelist,
        outpath = args.out,
        pout = args.pout,
        pin = args.pin)
    print("Demultiplexing complete.")

if __name__ == "__main__":
    main()
