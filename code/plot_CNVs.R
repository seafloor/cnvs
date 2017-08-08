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
# import files
args <- commandArgs(TRUE)
cnv.list.path <- args[1]
batch.dir <- args[2]
#cnv.list.path <- "~/Dropbox/PhD/Year 1/Rotation 3/code/cnvlist.txt"
#batch.dir <- "~/Desktop/testing cnv plot/BATCHES"


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

# create new variables
people.to.look.for <- unique(cnv.list$id)

# define plotting functions --------------------------------------------------------
PlotLrr <- function(id, start.position, end.position, chromosome, range) {
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
  axis(1, at = xlabs, labels = formatC(xlabs, format="d", big.mark=","))
  axis(2, at = c(-1, 1), labels = c("-1", "1"))
}

PlotBaf <- function(id, start.position, end.position, chromosome, range) {
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
  axis(1, at = xlabs, labels = formatC(xlabs, format="d", big.mark=","))
  axis(2, at = c(0, 1), labels = c("0", "1"))
}

# define functions for individuals ----------------------------------------

FindCNV <- function(file.path, cnv.list, people.to.look.for, id.present) {
  # FindCNV:
  #   - read in the individual's split file
  #   - iterate through regions to check, calling the plot functions for each
  #
  # Arguments:
  #   - 
  #
  # Returns:
  #   - 
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
  #print(short.list)
  for(line in 1:nrow(short.list)) {
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
    
    # filter for our cnv
    #range <- chromosome.file %>%
    #  filter(Chr == chromosome, Position >= start.position, Position <= end.position)
    range <- chromosome.file[which(chromosome.file$Chr == chromosome), ]
    range <- range[which(range$Position >= start.position), ]
    range <- range[which(range$Position <= end.position), ]
    
    #head(range)
    
    # call plotting functions
    write("calling plotting functions", stdout())
    PlotBaf(id, start.position, end.position, chromosome, range)
    PlotLrr(id, start.position, end.position, chromosome, range)
  }
}


# apply function ----------------------------------------------------------
#batch.num <- 1
#batch.parent <- "~/Desktop/testing cnv plot/BATCHES"
#batch.location <- paste("~/Desktop/testing cnv plot/BATCHES/Batch", batch.num, "/Batch", batch.num, "_split", sep = "")
batch.pattern <- "split\\..+"

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
    FindCNV(file.path, cnv.list, people.to.look.for, id.present)
  } else {
    #write(paste("ignoring file:", file.path, sep = " "), stdout())
  }
}))

mtext(paste("CNV traces generated using", batch.dir, "\non", format(Sys.time(), "%d-%m-%Y"), "using R version", getRversion(), sep = " "), outer = TRUE, side = 1)
dev.off() 
