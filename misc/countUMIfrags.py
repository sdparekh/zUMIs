#!/usr/bin/env python3
import pysam
import argparse
import multiprocessing as mp

def load_bcs(bcpath):
    with open(bcpath) as f:
        x = f.readline() # remove header
        y = f.readlines()
        bc = []
        for l in y:
            l = l.split(',')
            bc.append(l[0])
    return(bc)

def count_UMItags(inpath, bcs, chr):
    bccounts = {}
    for b in bcs:
        bccounts[b] = {}
        bccounts[b]['umi'] = 0
        bccounts[b]['int'] = 0
    inp = pysam.AlignmentFile(inpath, 'rb')
    for read in inp.fetch(chr):
        bc = read.get_tag('BC')
        if bc in bcs:
            ub = read.get_tag('UB')
            if ub is '':
                bccounts[bc]['int'] += 1
            else:
                 bccounts[bc]['umi'] += 1
    inp.close()
    return(bccounts)


def collect_write_stats(chrcounts, outpath):
    bccounts = chrcounts.pop(0) #get first dict
    for b in bccounts: #for every cell collect counts
        bccounts[b]['umi'] += sum( [chrcounts[i][b]['umi'] for i in range(len(chrcounts))] )
        bccounts[b]['int'] += sum( [chrcounts[i][b]['int'] for i in range(len(chrcounts))] )
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

    inp = pysam.AlignmentFile(args.bam, 'rb')
    chrs = list(inp.references)
    chrs.append('*') #don't forget unmapped reads
    inp.close()

    if(args.p > len(chrs)):
        num_threads = len(chrs)
    else:
        num_threads = args.p

    pool = mp.Pool(num_threads)
    results = [pool.apply_async(count_UMItags, (args.bam, bcs, chr, )) for chr in chrs]
    x = [r.get() for r in results]

    collect_write_stats(x, args.bcs+".BCUMIstats.txt")

if __name__ == "__main__":
    main()
