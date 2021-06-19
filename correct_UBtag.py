#!/usr/bin/env python3
import os
import pysam
import argparse
import multiprocessing as mp

def collect_bam_chunks(inpath, chrs, outpath):
    allpaths = [inpath+".tmp."+c+".bam" for c in chrs[:-1]]
    allpaths.append(inpath+".tmp."+"unmapped"+".bam")
    cat_args = ['-o', outpath]+allpaths
    pysam.cat(*cat_args)
    x = [os.remove(f) for f in allpaths]
    #pysam.index(outpath)

def load_bcs(bcpath):
    with open(bcpath) as f:
        x = f.readline() # remove header
        y = f.readlines()
        bc = []
        for l in y:
            l = l.split(',')
            bc.append(l[0])
    return(bc)

def load_dict(stub, bcs):
    molecules_dict = {}
    for i in bcs:
        fp = stub+i+".txt"
        if os.path.exists(fp):
            molecules_dict[i] = {}
            with open(fp) as f:
              x = f.readline() # remove header
              y = f.readlines()
              for l in y:
                l = l.strip().split('\t')
                if l[3] not in molecules_dict[i]:
                  molecules_dict[i][l[3]] = {}
                if l[0] not in molecules_dict[i][l[3]]:
                  molecules_dict[i][l[3]][l[0]] = {}
                molecules_dict[i][l[3]][l[0]] = l[1]
    return(molecules_dict)

# def return_UB(moldict, BC, GE, UX):
#     UB = UX
#     if BC in moldict:
#         if GE in moldict[BC]:
#             if UX in moldict[BC][GE]:
#                 UB = moldict[BC][GE][UX]
#     return(UB)

def return_UB(moldict, BC, GE, UX):
    try:
        UB = moldict[BC][GE][UX]
    except KeyError:
        UB = UX
    return(UB)

def correct_tags(inpath, threads, chr):
    global mols
    #nreads = 0
    if chr == '*':
        chrlabel = 'unmapped'
    else:
        chrlabel = chr
    outpath = inpath+".tmp."+chrlabel+".bam"
    inp = pysam.AlignmentFile(inpath, 'rb', threads = threads)
    out = pysam.AlignmentFile(outpath, 'wb', template = inp, threads = threads)
    for read in inp.fetch(chr):
        #nreads += 1
        umi = read.get_tag('UB')
        cell = read.get_tag('BC')
        if read.has_tag('GE'):
            gene = read.get_tag('GE')
        else:
            gene = 'NA'
        read.set_tag(tag = 'UX', value = umi, value_type = 'Z')
        umi_new = return_UB(moldict = mols, BC = cell, GE = gene, UX = umi)
        read.set_tag(tag = 'UB', value = umi_new, value_type = 'Z')
        out.write(read)
    inp.close()
    out.close()
    #print("Number of reads processed: "+nreads)

def main():
    parser = argparse.ArgumentParser(add_help=True)
    parser.add_argument('--bam', type=str, metavar='FILENAME',
                        help='Path to input BAM file')
    parser.add_argument('--out', type=str, metavar='FILENAME',
                        help='Path to output bam file')
    parser.add_argument('--p', type=int, default = 10,
                        help='Number of processes for bams')
    parser.add_argument('--bcs', type=str, metavar='FILENAME',
                        help='Path to kept barcodes')
    parser.add_argument('--stub', type=str, metavar='FILENAME',
                        help='Molecule table path stub')

    args = parser.parse_args()

    bcs = load_bcs(args.bcs)
    print("Loading molecule correction dictionary...")
    global mols
    mols = load_dict(args.stub, bcs)
    print("Correcting UB tags...")

    chrs = pysam.idxstats(args.bam).split('\n')
    chrs = [c.split('\t')[0] for c in chrs[:-1]]

    if args.p > 8:
        pysam_workers = 2
        n_jobs = int(args.p/2)
    else:
        pysam_workers = 1
        n_jobs = args.p

    pool = mp.Pool(n_jobs)
    results = [pool.apply_async(correct_tags, (args.bam,pysam_workers,chr, )) for chr in chrs]
    x = [r.get() for r in results]


    collect_bam_chunks(inpath = args.bam, chrs = chrs, outpath = args.out)

if __name__ == "__main__":
    main()
