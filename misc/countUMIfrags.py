#!/usr/bin/env python3
import pysam
import argparse

def load_bcs(bcpath):
    with open(bcpath) as f:
        x = f.readline() # remove header
        y = f.readlines()
        bc = []
        for l in y:
            l = l.split(',')
            bc.append(l[0])
    return(bc)

def count_UMItags(inpath, bcs, threads, outpath):
    bccounts = {}
    for b in bcs:
        bccounts[b] = {}
        bccounts[b]['umi'] = 0
        bccounts[b]['int'] = 0

    inp = pysam.AlignmentFile(inpath, 'rb', threads = threads)
    for read in inp:
        bc = read.get_tag('BC')
        if bc in bcs:
            ub = read.get_tag('UB')
            if ub is '':
                bccounts[bc]['int'] += 1
            else:
                 bccounts[bc]['umi'] += 1

    inp.close()
    with open(outpath, 'w') as out:
        out.write('XC\tnNontagged\tnUMItag\n')
        for bc in bccounts:
            out.write(bc+'\t'+str(bccounts[bc]['int'])+'\t'+str(bccounts[bc]['umi'])+'\n')



def main():
    parser = argparse.ArgumentParser(add_help=True)
    parser.add_argument('--bam', type=str, metavar='FILENAME',
                        help='Path to input BAM file')
    parser.add_argument('--p', type=int, default = 10,
                        help='Number of processes for bams')
    parser.add_argument('--bcs', type=str, metavar='FILENAME',
                        help='Path to kept barcodes')

    args = parser.parse_args()

    bcs = load_bcs(args.bcs)

    count_UMItags(inpath = args.bam, bcs = bcs, threads = args.p, outpath = args.bcs+".BCUMIstats.txt")

if __name__ == "__main__":
    main()
