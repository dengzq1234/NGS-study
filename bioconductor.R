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


# Print rnaframesZikaSet
rnaframesZikaSet

# Translate rnaframesZikaSet
AAzika6F <- translate(rnaframesZikaSet)
AAzika6F

# Count NS5 protein matches in AAzika6F, allowing 15 mismatches
vcountPattern(pattern = NS5, subject = AAzika6F, max.mismatch = 15)

# Subset the frame that contains the match from AAzika6F
selectedSet <- AAzika6F[3] 

# Convert selectedSet into a single sequence
selectedSeq <- unlist(selectedSet)

# Use vmatchPattern() with the set
vmatchPattern(pattern = ns5, subject = selectedSet, max.mismatch = 15)
"MIndex object of length 1
$`Pos 2`
IRanges object with 1 range and 0 metadata columns:
start       end     width
<integer> <integer> <integer>
[1]      3023      3347       325
# Take your time to see the similarities/differences in the result."


matchPattern(pattern=ns5, subject=selectedSeq, max.mismatch=15)
"matchPattern(pattern=ns5, subject=selectedSeq, max.mismatch=15)
Views on a 3597-letter AAString subject
subject: VVDLCESDCDSSSLKRELTTVSTGLIWIWKREFL...IDVGKTRDSMSFHHAGRQAQIAELRRPVWGNPWF
views:
start  end width
[1]  3023 3347   325 [SRAIWYMWLGARFLEFEALGFLNEDHW...HRRDLRLMANAICSAVPVDWVPTGRTT]"


# Load IRanges package
library(IRanges)

# IRnum1: start - vector 1 through 5, end - 100  
IRnum1 <- IRanges(1:5, end=100)

# IRnum2: end - 100, width - 89 and 10
IRnum2 <- IRanges(width=c(89,10), end=100)

# IRlog1: start = Rle(c(F, T, T, T, F, T, T, T)))
IRlog1 <- IRanges(start = Rle(c(F, T, T, T, F, T, T, T)))

# Print objects in a list
print(list(IRnum1 = IRnum1, IRnum2 = IRnum2, IRlog1 = IRlog1))

# Create the first sequence seq_1
seq_1 <- IRanges(start = 10, end = 37)

# Create the second sequence seq_2
seq_2 <- IRanges(start = c(5, 35, 50),
                 end = c(12, 39, 61),
                 names = LETTERS[1:3])

# Check the width of seq_1 and seq_2
width(seq_1)
width(seq_2)

# Check the width of seq_1 and seq_2
lengths(seq_1)
lengths(seq_2)

# Load GenomicRanges package
library(GenomicRanges)

# Print seq_intervals
seq_intervals

# Create myGR
myGR <- as(seq_intervals, "GRanges")

# Print myGR
myGR

# Print the seqinfo of myGR
seqinfo(myGR)

# Check the metadata
mcols(myGR)

# Load human reference genome hg38
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

# Assign hg38 to hg, then print it
hg <- TxDb.Hsapiens.UCSC.hg38.knownGene
hg

# Extract all the genes in chromosome X as hg_chrXg, then print it
hg_chrXg <- genes(hg, filter = list(tx_chrom = c("chrX")))
hg_chrXg

# Extract all positive stranded genes in chromosome X, assign to hg_chrXgp, then sort it
hg_chrXgp <- genes(hg, filter = list(tx_chrom = c("chrX"), tx_strand = "+"))
sort(hg_chrXgp)

# Store the overlapping range in rangefound
rangefound <- subsetByOverlaps(hg_chrX, ABCD1)

# Print names of rangefound
names(rangefound)

# Print the gene of interest 
GRanges(rangefound)

# Print rangefound
rangefound


# Load the human transcripts DB to hg
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
hg <- TxDb.Hsapiens.UCSC.hg38.knownGene

# Prefilter chromosome X "chrX" using seqlevels()
seqlevels(hg) <- c("chrX")

# Get all transcripts by gene and print it
hg_chrXt <- transcriptsBy(hg, by = "gene")
hg_chrXt

# Select gene `215` from the hg_chrXt
hg_chrXt$'215'

library(ShortRead)

# Print fqsample
print(fqsample)

# Print fqsample
fqsample

# Check class of fqsample
class(fqsample)

# Check class sread fqsample
class(sread(fqsample))

# Check ids of fqsample
id(fqsample)

f<-"/usr/local/share/datasets/small_SRR1971253.fastq"

# Set a seed for sampling
set.seed(1234)

# Use FastqSampler with f and select 100 reads
fs <- FastqSampler(con = f, n = 100)

# Generate new sample yield
my_sample <- yield(fs)

# Print my_sample
my_sample

# Reverse the zikv sequence
reverse(zikv)

# Complement the zikv sequence
complement(zikv)

# Reverse complement the zikv sequence
reverseComplement(zikv)

# Translate the zikv sequence
translate(zikv)
