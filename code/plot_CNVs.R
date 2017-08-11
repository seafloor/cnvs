#!/usr/bin/env Rscript
# header ------------------------------------------------------------------
# July 2017
# Matt Smith
#
# Outline:
# read in the list of cnvs to visualise, columns = c("id", "chromosome", "start", "end")
# read in the IDs file, search for individual to get 
# sort list of cnvs by chromosome then by start position
# for each unique chromosome, call function that:
#   - read in the chromosomes snp file
#   - search for start and end locations to get row numbers
#   - also store a vector of positions from the start to the end
#   - search for the individual in the individual IDs file to get column number
#   - read in only the needed rows from that column, using skip and n_max etc. for both baf and l2r
#   - make two plots - the vector of positions on x, baf or l2r on y
#   - try to append plots to a pdf
#
# scripts (try to) follow Google's R style guide: https://google.github.io/styleguide/Rguide.xml

# setup -------------------------------------------------------------------
options(warn = 2)
#begin <- Sys.time()
# set default values
plot_region <- "file"
all_chr <- vector()
all_start <- vector(mode = "integer")
all_end <- vector(mode = "integer")
# out <- ""

# header to be print if incorrent number of arguments given
PrintDescription <- function() {
  write(paste("## ----------------------------------- plot_CNVs ----------------------------------- ##",
    "##\n##\n## plot_CNVs creates a pdf of LRR and BAF plots for CNV regions you supply",
    "## you have 3 options for plotting CNVs:",
    "##\n## 1. a basic call to plot_CNVs looks like '.plot_CNVs cnvlist.txt directory_to_check', where:",
    "##      * cnvlist.txt has the 4 headers id, chromosome, start and end, and the id doesn't include the split3. section",
    "##      * directory_to_check is any directory that contains split files (they can be in subfolders). Don't include a slash at the end.",
    "##      * e.g. './plot_CNVs ./cnvstocheck.txt /c8000xd3/uk-biobank-kirov/BATCHES'",
    "##\n## 2. you can override cnv locations for IDs by calling './plot_CNVs cnvlist.txt dir_to_check override chr start end'",
    "##      * e.g. './plot_CNVs ./testlist_override.txt /c8000xd3/uk-biobank-kirov/BATCHES override 22 44000000 46000000'",
    "##      * this will ignore any cnv locations in the testlist_override.txt file",
    "##      * if you're overriding values, your testlist_override.txt file can just have the one header, 'id', with just a single column of IDs",
    "##\n## 3. you can plot the region you supply in the command line, AND have the actual CNV region from the cnvlist.txt file annotated",
    "##      * e.g. './plot_CNVs ./testlist.txt /c8000xd3/uk-biobank-kirov/BATCHES both 22 44000000 46000000'",
    "##      * make sure you type 'both', and not 'override'",
    "##      * you obviously need to supply a full cnvlist.txt file as described in 1",
    "##\n## the CNV_visualisations.pdf output is saved in the directory you supply\n##\n##",
    "## ----------------------------------- enjoy :) ----------------------------------- ##\n", sep = "\n"), stdout())
}

# read in command line arguments, catch errors in input, return PrintDescription if inputs are incorrect
args <- commandArgs(TRUE)
if(length(args) < 2 | length(args) > 6) {
  PrintDescription()
  stop("Not enough arguments. Please follow instructions above")
} else if(length(args) == 2) {
  # set command-line args
  cnv.list.path <- args[1]
  if(!file_test("-f", cnv.list.path)) {
    stop("cnv.list.path is not a file")
  }
  batch.dir <- args[2]
  if(!file_test("-d", batch.dir)) {
    stop("batch.dir is not a directory")
  }
} else if(length(args) > 2 & length(args) < 6) {
  PrintDescription()
  stop("Incorrect number of arguments. Please follow instructions above")
} else if(length(args) == 6) {
  # setup for overriding cnv locations or plotting both
  cnv.list.path <- args[1]
  if(!file_test("-f", cnv.list.path)) {
    stop("cnv.list.path is not a file")
  }
  batch.dir <- args[2]
  if(!file_test("-d", batch.dir)) {
    stop("batch.dir is not a directory")
  }
  plot_region <- args[3]
  if(!plot_region %in% c("override", "both")) {
    stop("plot_region argument not recognised. Please use either 'override' or 'both' without quotation marks")
  }
  all_chr <- args[4]
  if(!all_chr %in% c(1:22, "X", "Y")) {
    stop("all_chr is not a recognised chromosome")
  }
  tryCatch ({
    all_start <- as.integer(args[5])
    all_end <- as.integer(args[6])
    },
    error = function(c) {
      c$message <- ("Error: start and end positions must be integers")
      stop(c)
    }
  ) 
  if(plot_region == "override") {
    write(paste("Overriding CNV locations in csv file with chromosome ", all_chr, ", starting location ", all_start, ", ending location ", all_end, sep = ""), stdout())
  } else if(plot_region == "both") {
    write(paste("Plotting BOTH CNV locations in csv file AND chromosome ", all_chr, ", starting location ", all_start, ", ending location ", all_end, sep = ""), stdout())
  }
}

# import files
if(plot_region %in% c("file", "both")) {
  cnv.list <- read.table(cnv.list.path,
                       header = TRUE,
                       stringsAsFactors = FALSE,
                       colClasses = c("character", "character", "integer", "integer"))
  # checking start and end positions are in the right order
  apply(cnv.list, 1, function(x) {
    if(x[3] >= x[4]) {
      stop("The start of a CNV must come before the end!")
    }
  })
} else if(plot_region == "override") {
  cnv.list <- read.table(cnv.list.path,
			 header = TRUE,
			 stringsAsFactors = FALSE,
			 colClasses = c("character"))
}

# if plot_region = both, check that the outer regions are more extreme than the inner regions in the file
if(plot_region == "both") {
  if(sum(sapply(cnv.list$start, function(x) x < all_start)) > 0) {
    stop("start region passed to command line must be lower than all start regions in the cnvlist file")
  }
  if(sum(sapply(cnv.list$end, function(x) x > all_end)) > 0) {
    stop("end region passed to command line must be lower than all start regions in the cnvlist file")
  }
}

# create new variables
people.to.look.for <- unique(cnv.list$id)

# define plotting functions --------------------------------------------------------
PlotLrr <- function(id, start.position, end.position, chromosome, range, inner.chromosome, inner.start.position, inner.end.position, plot_region) {
  # PlotCNV:
  #   - make two plots - the vector of positions on x, baf or lrr on y
  #
  # Arguments:
  #   - positions
  #   - baf
  #   - lrr
  #
  # Returns:
  #   - plot?
  #print("plotlrr called")
  xlabs <- round(seq(from = start.position, to = end.position, length.out = 5), digits = -2)
  plot(range[, "Position"], range[, "LRR"], type = "p", col = "blue", pch = 20, xaxt = "n", yaxt = "n", xlab = "", ylab = "", ylim = c(-1, 1))
  title(main = paste("Chromosome ", 
                     chromosome,
                     ", ID ", id,
                     ", position ",
                     formatC(start.position, format="d", big.mark=","),
                     " to ", formatC(end.position, format="d", big.mark=","),
                     sep = ""),
        xlab = "Location (bp)",
        ylab = "LRR")
  abline(h = 0, col = "grey", lty = 2)
  if(plot_region == "both") {
    abline(v = inner.start.position, col = "plum3", lty = 1)
    abline(v = inner.end.position, col = "plum3", lty = 1)
  }
  axis(1, at = xlabs, labels = formatC(xlabs, format="d", big.mark=","))
  axis(2, at = c(-1, 1), labels = c("-1", "1"))
}

PlotBaf <- function(id, start.position, end.position, chromosome, range, inner.chromosome, inner.start.position, inner.end.position, plot_region) {
  #print("plotbaf called")
  xlabs <- round(seq(from = start.position, to = end.position, length.out = 5), digits = -2)
  plot(range[, "Position"], range[, "BAF"], type = "p", col = "blue", pch = 20, xaxt = "n", yaxt = "n", xlab = "", ylab = "", ylim = c(0, 1))
  title(main = paste("Chromosome ", 
                     chromosome,
                     ", ID ", id,
                     ", position ",
                     formatC(start.position, format="d", big.mark=","),
                     " to ", formatC(end.position, format="d", big.mark=","),
                     sep = ""),
        xlab = "Location (bp)",
        ylab = "BAF")
  abline(h = 0.5, col = "grey", lty = 2)
  if(plot_region == "both") {
    abline(v = inner.start.position, col = "plum3", lty = 1)
    abline(v = inner.end.position, col = "plum3", lty = 1)
  }
  axis(1, at = xlabs, labels = formatC(xlabs, format="d", big.mark=","))
  axis(2, at = c(0, 1), labels = c("0", "1"))
}

# define functions for individuals ----------------------------------------

FindCNV <- function(file.path, cnv.list, people.to.look.for, id.present, all_chr, all_start, all_end, plot_region) {
  # FindCNV:
  #   - read in the individual's split file
  #   - iterate through regions to check, calling the plot functions for each
  #
  # Arguments:
  #   - 
  #
  # Returns:
  #   - 
  
  # define empty vectors for holding inner regions if "both" is called
  inner.chromosome <- vector(mode = "character")
  inner.start.position <- vector(mode = "integer")
  inner.end.position <- vector(mode = "integer")

  chromosome.file <- tryCatch(read.table(file.path,
                              header = FALSE,
                              skip = 1,
                              col.names = c("Name", "Chr", "Position", "LRR", "BAF"),
                              stringsAsFactors = FALSE,
                              colClasses = c("character", "character", "integer", "double", "double")), error = function(c) {
	  c$message <- paste("Error when trying to open file", file.path, sep = " ")
	  stop(c)
  })
  chromosome.file <- chromosome.file[order(chromosome.file$Chr, chromosome.file$Position), ]
  
  if(sum(id.present) > 1) {
    stop("File name matches more than one ID")
  }
  
  id <- people.to.look.for[id.present]
  short.list <- cnv.list[cnv.list$id == id, ]
  if(!is.data.frame(short.list)) {
    looplen <- 1:length(short.list)
  } else {
    looplen <- 1:nrow(short.list)
  }

  if(plot_region %in% c("both", "override")) {
    # if user wants to override the csv file locations for CNVs
    chromosome <- all_chr
    # add start and end positions
    # get closest positions if SNP locations are not present in the file
    if(all_start %in% chromosome.file[, "Position"]) {
      start.position <- all_start
    } else {
      write("start snp not present: choosing closest match", stdout())
      start.position <- chromosome.file[, "Position"][which.min(abs(chromosome.file[, "Position"] - all_start))]
    }
    
    if(all_end %in% chromosome.file[, "Position"]) {
      end.position <- all_end
    } else {
      write("end snp not present: choosing closest match", stdout())
      end.position <- chromosome.file[, "Position"][which.min(abs(chromosome.file[, "Position"] - all_end))]
    }
  }
  #print(short.list)
  for(line in looplen) {
    if(plot_region == "both") {
      inner.chromosome <- short.list[line, "chromosome"]
      
      # get closest positions if SNP locations are not present in the file
      if(short.list[line, "start"] %in% chromosome.file[, "Position"]) {
        inner.start.position <- short.list[line, "start"]
      } else {
        write("start snp not present: choosing closest match", stdout())
        inner.start.position <- chromosome.file[, "Position"][which.min(abs(chromosome.file[, "Position"] - short.list[line, "start"]))]
      }
      
      if(short.list[line, "end"] %in% chromosome.file[, "Position"]) {
        inner.end.position <- short.list[line, "end"]
      } else {
        write("end snp not present: choosing closest match", stdout())
        inner.end.position <- chromosome.file[, "Position"][which.min(abs(chromosome.file[, "Position"] - short.list[line, "end"]))]
      }
    } else if(plot_region == "file") {
      # if user just supplies csv file locations
      chromosome <- short.list[line, "chromosome"]
      
      # get closest positions if SNP locations are not present in the file
      if(short.list[line, "start"] %in% chromosome.file[, "Position"]) {
        start.position <- short.list[line, "start"]
      } else {
        write("start snp not present: choosing closest match", stdout())
        start.position <- chromosome.file[, "Position"][which.min(abs(chromosome.file[, "Position"] - short.list[line, "start"]))]
      }
      
      if(short.list[line, "end"] %in% chromosome.file[, "Position"]) {
        end.position <- short.list[line, "end"]
      } else {
        write("end snp not present: choosing closest match", stdout())
        end.position <- chromosome.file[, "Position"][which.min(abs(chromosome.file[, "Position"] - short.list[line, "end"]))]
      }
    }
    
    # filter for our cnv
    #range <- chromosome.file %>%
    #  filter(Chr == chromosome, Position >= start.position, Position <= end.position)
    range <- chromosome.file[which(chromosome.file$Chr == chromosome), ]
    range <- range[which(range$Position >= start.position), ]
    range <- range[which(range$Position <= end.position), ]
    
    #head(range)
    
    # call plotting functions
    write("calling plotting functions", stdout())
    PlotBaf(id, start.position, end.position, chromosome, range, inner.chromosome, inner.start.position, inner.end.position, plot_region)
    PlotLrr(id, start.position, end.position, chromosome, range, inner.chromosome, inner.start.position, inner.end.position, plot_region)
  }
}
 

# apply function ----------------------------------------------------------
#batch.num <- 1
#batch.parent <- "~/Desktop/testing cnv plot/BATCHES"
#batch.location <- paste("~/Desktop/testing cnv plot/BATCHES/Batch", batch.num, "/Batch", batch.num, "_split", sep = "")
batch.pattern <- "^split[a-zA-Z0-9]+\\..+"

pdf.name = paste(batch.dir, "CNV_visualisations.pdf", sep = "/") # choose pdf name
pdf(file = pdf.name, width = 8.27, height = 11.69) 
par(mfrow = c(6, 1),
    mar = c(2, 5, 2, 1),
    oma = c(3, 0, 4, 0))

# for each batch number, get the folder location and call FindCNV
write("listing files in subdirectories", stdout())
tryCatch(files <- list.files(path = batch.dir, pattern = batch.pattern, full.names=T, recursive=TRUE), error = function(c) {
	 c$message <- paste("failed to list files")
	 stop(c)
})
if(length(files) > 0) {
  write(paste(length(files), "files listed for checking", sep = " "), stdout())
} else {
  stop("no .split files found in this folder")
}

#print(files)
invisible(
lapply(files, function(file.path) {
  #print(file.path)
  id.present <- unlist(lapply(people.to.look.for, function(x) grepl(x, file.path)))
  
  #print(id.present)
  if(sum(id.present) > 0) {
    # if the file name contains one of the IDs we want
    # call FindCNV to read in the file and plot all their CNVs we want
    write(paste("checking file:", file.path, sep = " "), stdout())
    FindCNV(file.path, cnv.list, people.to.look.for, id.present, all_chr, all_start, all_end, plot_region)
  } else {
    #write(paste("ignoring file:", file.path, sep = " "), stdout())
  }
}))

mtext(paste("\nCNV traces generated using", batch.dir, "on", format(Sys.time(), "%d-%m-%Y"), "using R version", getRversion(), sep = " "), outer = TRUE, side = 1, cex = 0.6)
dev.off() 
#finish <- Sys.time()
#print(paste("Time taken to plot CNVs:", (finish - begin), sep = " "))
