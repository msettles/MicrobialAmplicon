###########################################################################################
###########################################################################################
# PARAMTER SECTION FOR PREPROCESSING 16S MICROBIAL AMPLICONS (CURRENTLY 454 ONLY)
###########################################################################################
###########################################################################################
### TODO: need to split out basefilename required options
### VERSION OF THIS CODE
my_ver <- "2.0"
### cross_match parameters
cross_match_ver <- "1.090518"
cross_match_minmatch <- 8
cross_match_minscore <- 12
cross_match_screenfile <- file.path(microbe.amplicon.home,"ext.data","screen_27f-534r.combined.fa")
### for multicore capabilities
cross_match_call <- function(x) paste("cross_match ",basefilename,".",x,".raw.fasta ", cross_match_screenfile, " -minmatch ",cross_match_minmatch," -minscore ", cross_match_minscore," -tags > ",basefilename,".",x,".cmout",sep="")

### screenfile properties
tagkey <- "^FP_"
tag_nucs <- 8
primerkey <- "^RP_"

### lucy parameters
lucy_ver <- "1.20p"
lucy_max_avg_error <- 0.002 # Qscore 27
lucy_max_error_at_ends <- 0.002
lucy_call <- paste("lucy -xtra ", nproc," -minimum 0 -debug ",basefilename,".lucy_clip.txt -error ",lucy_max_avg_error, " ", lucy_max_error_at_ends, " -output ",basefilename,".lucy.fasta ",basefilename,".lucy.fasta.qual ",basefilename,".adapterClip.fasta ",basefilename,".adapterClip.fasta.qual",sep="")

### mothur parameters

mothur_ver <- "1.27.0"
mothur_alignment_db <- "silva.bacteria.fasta"
mothur.template="/mnt/home/msettles/projects/Amplicon_Preprocessing/Alignment_db/silva.bacteria.fasta"
mothur_align_call <- paste("mothur \"#align.seqs(candidate=",basefilename,".rdp.fasta, template=", mothur.template ,", flip=T, processors=",nproc,")\"",sep="")
#mothur_filter_call <- paste("mothur \"#filter.seqs(fasta=",basefilename,".rdp.align, processors=",nproc,");\"",sep="")
TE_exp <- 534
align_length_max_error <- 5
TE_max_dist <- 75

### Ribosomal database project
rdp_ver <- 2.5
rdp_path <- "/mnt/home/msettles/opt/rdp_classifier_2.5/rdp_classifier-2.5.jar"
rdp_call <- function(x) paste("java -Xmx1g -jar ",rdp_path," -q ",basefilename,".",x,".rdp.fasta -o ",basefilename,".",x,".lucy.rdpV6.fix -f fixrank",sep="")
reverseSeq <- TRUE

## filter parameters
minlength <- 350  ## minimum length of acceptable sequence
maxlength <- 600  ## maximum length of sequence allowed
maxhammingdisttag <- 1 ## max hamming distance allowed for barcode
maxforwardprimererrors <- 2 ## max number of errors allowed in the forward primer
maxNs <- 2 ## max number of Ns post lucy filtering
maxhomopol <- 10  ## maximum homopolymer allowed

version = paste("Rcode:",my_ver,";rdp:",rdp_ver,";mothur:",mothur_ver,";alignment_db:",mothur_alignment_db,collapse="")
