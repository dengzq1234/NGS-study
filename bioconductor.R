# install bioconductor
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
# BiocManager::install()
# BiocManager::available()
BiocManager::valid()
library(BiocManager)


BiocManager::install("BSgenome")
packageVersion("BSgenome")
library(BSgenome)
BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
library(BSgenome.Hsapiens.UCSC.hg38)
BiocManager::install("BSgenome.Scerevisiae.UCSC.sacCer3")
library(BSgenome.Scerevisiae.UCSC.sacCer3)
help("BSgenome")


setwd("/home/deng/Educations/rna-seq/")

A <- "try1"
isS4(A)
str(A)

showClass("BSgenome")
metadata(Scerevisiae)
a_genome <- Scerevisiae

# What is a_genome's main class?
class(a_genome)  # "BSgenome"

# What is a_genome's other classes?
is(a_genome)  # "BSgenome", "GenomeDescription"

# Is a_genome an S4 representation?
isS4(a_genome)  # TRUE

show(a_genome)
organism(a_genome)
provider(a_genome)
seqinfo(a_genome)

yeastGenome <- BSgenome.Scerevisiae.UCSC.sacCer3

available.genomes()

getSeq(yeastGenome, 'chrM')

getSeq(yeastGenome, end=10)


# Get the head of seqnames and tail of seqlengths for yeastGenome
head(seqnames(yeastGenome))
tail(seqlengths(yeastGenome))

yeastGenome$chrM

# Count characters of the chrM sequence
nchar(yeastGenome$chrM)


# Get the first 30 bases of chrM
getSeq(yeastGenome, names = "chrM", start=1, end = 30)

library(Biostrings)

dna_seq <- DNAString("ATGATCTCGTAA")

zikaVirus <- readDNAStringSet('data/zikavirus.fa')
# Create zikv with one collated sequence using `zikaVirus`
zikv <- unlist(zikaVirus)
zikv

# Check the length of zikaVirus and zikv
length(zikaVirus)
length(zikv)

# Check the width of zikaVirus
width(zikaVirus)
# Subset zikv to only the first 30 bases
subZikv <- subseq(zikv, end = 30)
subZikv

# Reverse the zikv sequence
reverse(zikv)

# Complement the zikv sequence
complement(zikv)

# Reverse complement the zikv sequence
reverseComplement(zikv)

# Translate the zikv sequence
translate(zikv)
