#! /usr/bin/Rscript

suppressMessages(library(optparse))

option_list = list(
  make_option(c("-o", "--outdir"), type="character", default=NULL, help="Path to output dir", metavar="character"),
  make_option(c("-i", "--indir"), type="character", default=NULL, help="Path to the input dir (Should have the counts file)", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$i) || is.null(opt$o)){
  print_help(opt_parser)
  stop("Input Parameters Not Accurate", call.=FALSE)
}

message("Input Dir is ", opt$i)
message("Ouput Dir is ", opt$o)

basedir <- opt$i
setwd(basedir)
outdir <- opt$o

pat <- ".count"

myfiles <- list.files(path=basedir, pattern="^[ACGT].*.count")
print(myfiles)

# read each file as array element of DT and rename the last 2 cols
# we created a list of single sample tables
DT <- list()
for (i in 1:length(myfiles) ) {
  infile = paste(basedir, myfiles[i], sep = "/")
  DT[[myfiles[i]]] <- read.table(infile, header = F, stringsAsFactors = FALSE)
  cnts <- gsub("(.*).count", "\\1", myfiles[i])
  colnames(DT[[myfiles[i]]]) <- c("ID", cnts)
}

# merge all elements based on first ID columns
data <- DT[[myfiles[1]]]

# we now add each other table with the ID column as key
for (i in 2:length(myfiles)) {
  y <- DT[[myfiles[i]]]
  z <- merge(data, y, by = c("ID"), all=TRUE)
  data <- z
}

# ID column becomes rownames
rownames(data) <- data$ID
data <- data[,-1]

# write data to file
write.csv(data, file = paste(outdir, "CRISPR_counts_all_Raw1.csv", sep ="/"))

