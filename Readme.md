This is an R package that compares a user's single nucleotide variant list to calls derived from ultra-deep whole-genome and validation data generated from the AML31 tumor, relapse, and matched normal.

##Installation

    #install dependencies
    install.packages(c("venneuler","devtools"))
    library(devtools)
    install_github("genome/aml31Benchmarking")


## Usage
    library(aml31Benchmarking)
    benchmark("snvlist","outTable","out.pdf")



Descriptions of how the 'gold' and 'platinum' lists were generated can be found in this publication:
(insert publication info here).

Sequence data for those who wish to make their own variant calls is available through dbGaP (insert accession numbers here)
