## Cell Barcodes

In order to be compatible with well-based and droplet-based scRNA-seq methods, zUMIs is flexible with handling of cell barcodes.
As default behavior, zUMIs tries to guesstimate the relevant barcodes from the data by finding a cluster with the highest number of reads using model based clustering.


![](https://github.com/sdparekh/zUMIs/blob/master/ExampleData/zUMIs_output/stats/detected_BC.png)


To override automatic detection of barcodes, users can either give a fixed number of barcodes to consider (e.g. "-b 100") or refer to a plain text file containing known expected barcodes (e.g. "-b barcodefile.txt").
The text file should contain just one column with a list of barcodes without headers and without sample names:

```
ATGAAT
ATCAAA
GGAGCC
TAAGAT
AAAACT
GCGCTG
CCAACC
CTTTAA
TCATAT
TACTAT
```