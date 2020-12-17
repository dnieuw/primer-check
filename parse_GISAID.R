library(ape)
library(Biostrings)
library(data.table)
library(stringr)
library(lubridate)

setwd("D://OneDrive/Rprojects/primer-check-v2/")

#May need some changes if GISAID decides to change the format of their tar file
ALIGNMENT_TAR <- list.files(pattern = ".tar.xz")
ALIGNMENT_FILE <- paste(str_remove(ALIGNMENT_TAR, ".tar.xz"),str_replace(ALIGNMENT_TAR,".tar.xz",".fasta"), sep="/")

#Untar file (this xzfile stuff is needed on windows, otherwise just untar)
zz <- xzfile(ALIGNMENT_TAR, open = "rb")
untar(zz, files = ALIGNMENT_FILE)
close(zz)

#Read the complete file into memory (works for ~200k sequences, may need a better solution in the future)
aln <- readDNAStringSet(ALIGNMENT_FILE, format = "fasta")
#Define and capture the reference sequences (Wuhan WH01)
ref <- aln[grep("EPI_ISL_406798",names(aln)),]

#Remove sequences that are too short
aln_width <- median(width(aln))
aln <- aln[which(width(aln)>=aln_width),]

#Get metadata from sequence headers by splitting on | and /
metadata <- as.data.table(cbind(names(aln),
                                str_split(str_split(names(aln),'\\|', simplify = T)[,1], "/", simplify = T),
                                str_split(names(aln),'\\|', simplify = T)[,-1]))
colnames(metadata) <- c("seq_id","hcov","country","note","year","epi","date","continent")
setkey(metadata, "seq_id")

#Fix some dates
metadata[,date:=str_remove_all(date, '-00')]
metadata[,date:=ymd(date,truncated = 2)]

#Initialize reduction routine
progress <- 0
total <- length(metadata[,.N,date][,date])
days <- metadata[,.N,date][,date]

#Function for clustering and reducing sequences by day/country at a 5nt difference threshold
#This take time (an hour or so)
reduce_daily_sequences <- function(day) {
  daily <- metadata[date==day]
  #remove sequences which are too short or long....
  daily <- daily[width(aln[daily$seq_id,])==aln_width,]
  
  #Return if only 1 sequence was produced on that day
  if (nrow(daily)==1) {
    progress <<- progress + 1
    cat('\r',paste(progress,total, sep = "/"))
    flush.console() 
    return(daily)
  }
  
  #Transform sequence selection to DNAbin class for further analysis
  daily_subset <- as.DNAbin(as.matrix(aln[daily$seq_id,]))
  
  #extract variable sites from the sequence
  variable_sites <- which(sapply(1:ncol(daily_subset), function(x) !any(base.freq(daily_subset[,x])>0.99)))
  
  #Return if all sequences on that day have no variable sites
  if (length(variable_sites)==0) {
    progress <<- progress + 1
    cat('\r',paste(progress,total, sep = "/"))
    flush.console() 
    return(daily)
  }
  
  #Calculate all vs all nucleotide distance based on variable sites
  variable_dist <- dist.dna(daily_subset[,variable_sites], model = "N", pairwise.deletion = T)
  
  #Cluster based on 5 nucleotide difference
  seq_clusters <- cutree(hclust(variable_dist),h = 5)
  
  daily[,cluster:=seq_clusters]
  
  return_most_common <- function(seq_id) {
    seq_dist <- as.matrix(variable_dist)[seq_id,,drop=F]
    sum_seq_dist <- rowSums(seq_dist)
    most_common <- order(sum_seq_dist)[1]
    return(most_common)
  }
  
  #Find most common (least total mutation difference) sequence per cluster and per country
  daily <- daily[,.SD[return_most_common(seq_id)],.(cluster,country)]
  daily <- daily[,.(seq_id,hcov,country,note,year,epi,date, continent)]
  
  progress <<- progress + 1
  cat('\r',paste(progress,total, sep = "/"))
  flush.console() 
  
  return(daily)
}

#Combine all results
result <- rbindlist(lapply(days, reduce_daily_sequences))

#Save (for testing)
# save(result, file="reduction_data.Rdata")
# load("reduction_data.Rdata")

#Keep only the sequences in the daily reduction result and only from the past 2 months
aln <- aln[result[date>today()-months(2),seq_id]]

#Add the reference sequence tot the start of the alignment
aln <- c(ref,aln)

#Remove gaps from the alignment if >90% are gaps at that position
aln <- DNAStringSet(maskGaps(DNAMultipleAlignment(aln), min.fraction=0.9, min.block.width=1))

#Write reduced data to file
writeXStringSet(aln, "example_data.fasta", format="fasta")

#Remove tar file and msa file
# unlink(c(ALIGNMENT_FILE,ALIGNMENT_TAR))

#Reupload and update the app 
# rsconnect::deployApp('.', account = "viroscience-emc", appName = "primer-check", appTitle = "SARS-CoV-2 primer check visualization")
