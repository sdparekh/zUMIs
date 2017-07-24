# Welcome to zUMIs :red_car::dash:

zUMIs is a fast and flexible pipeline to process RNA-seq data with UMIs.

The input to this pipeline is paired-end fastq files, where one read contains the cDNA sequence and the other read contains UMI and Cell Barcode information. Furthermore, you will need a STAR index for your genome (see below).

![zUMIs Workflow](https://github.com/sdparekh/zUMIs/blob/master/zUMIs.png?raw=true)

You can read more about zUMIs in our [biorxiv preprint](http://www.biorxiv.org/content/early/2017/06/22/153940)!

## Installation and Usage

Please find information on [installation](https://github.com/sdparekh/zUMIs/wiki/Installation) and [usage](https://github.com/sdparekh/zUMIs/wiki/Usage) in the [zUMIs wiki](https://github.com/sdparekh/zUMIs/wiki/).

## Compatibility

zUMIs is compatible with these single-cell UMI protocols:

- CEL-seq with UMI (Grün et al., 2014)
- SCRB-seq (Soumillon et al., 2014)
- MARS-seq (Jaitin et al., 2014)
- STRT-C1 (Islam et al., 2014)
- Drop-seq (Macosko et al., 2015)
- CEL-seq2 (Hashimshony et al., 2016)
- SORT-seq (Muraro et al., 2016)
- DroNc-seq (Habib et al., 2017)
- SPLiT-seq (Rosenberg et al., 2017)
- STRT-2i (Hochgerner et al., 2017)
- Quartz-seq2 (Sasagawa et al., 2017)
- 10x Genomics Chromium (Zheng et al., 2017)

For InDrops compatibility, users need to preprocess the barcode and UMI read because of variable length cell barcodes.

## Getting help

Refer to [zUMIs Github wiki](https://github.com/sdparekh/zUMIs/wiki) for help.

Please report bugs :beetle::bug: to the [zUMIs Github issue page](https://github.com/sdparekh/zUMIs/issues)

If you encounter issues when using zUMIs for the first time, please try to [run the example data set](https://github.com/sdparekh/zUMIs/wiki/Usage) included in this repository.
