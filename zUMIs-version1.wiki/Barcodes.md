zUMIs provides three main options for selecting relevant barcodes controlled by the `-b` switch:

* automatic detection
* number of barcodes with most reads
* barcode list annotation

Here is more information on each of the modes:
## Automatic barcode detection
zUMIs infers which barcodes mark good cells from the observed sequences. To this end, we fit a k-dimensional multivariate normal distribution using the R-package mclust for the number of reads/BC, where k is empirically determined by mclust via the Bayesian Information Criterion (BIC). We reason that only the kth normal distribution with the largest mean contains barcodes that identify reads originating from intact cells. We exclude all barcodes that fall in the lower 1% tail of this kth normal-distribution to exclude spurious barcodes.

## Number of barcodes with most reads
zUMIs will make a summary statistic over all observed barcode sequences and their frequency. The user-specified number of barcodes will be selected in descending order.

e.g. `-b 1000`

## Barcode annotation
If expected barcodes are known a priori, it is usually advisable to provide these.
The format should be a plain text file without headers, where each line contains the exact barcode sequence.

For instance:
> GGGGCA
>
> TATTGT
>
> GCACGG
>
> CAATAA
>
> CGCGTG

Attention: If you have specified a 6-mer in the barcode range (eg. `-c 1-6`), this annotation should also contain 6-mer reference barcodes!
 
In case you are using an additional plate barcode in zUMIs, the expected barcodelist should contain the concatenated string of all possible expected plate+cell barcode combinations!

`[plateBC][cellBC]`

For instance, take the above cell barcodes that should all have the same plate barcode:
> CGTACTAGGGGGCA
>
> CGTACTAGTATTGT
>
> CGTACTAGGCACGG
>
> CGTACTAGCAATAA
>
> CGTACTAGCGCGTG

Attention: Make sure the annotation always contains reference barcodes with correct length (sum of plate and cell barcode lengths)!
