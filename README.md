### Summary and License
This script was used to obtain the coverage plots depicted in Figure 2 of the Ponsard _et al_ paper ; "Conserved and distinct expression of circular RNAs in commercially used Marek’s disease vaccine viruses."

If you use parts of this script, cite the paper, as this script is protected under a GPLv3 license.

## Remarks
### !!! This pipeline only works after modifications of lines 4-5 ; 150-151 in the Python script !!!
In these lines, you should indicate the files location of the vCircTrappist output.

### !!! Same remark for lines 1-3 in the R script !!!

Coverages were obtained using the command
`samtools depth -a input.sam > output.csv`
on all the alignment files.
