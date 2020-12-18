library(stringr)
library(plyr)
library(openxlsx)

data.dir <- "~/Dropbox (Partners HealthCare)/Melanoma P2 Analysis/patient15/tcr_data/"
setwd(data.dir)
samples <- list.files(data.dir)
print(samples)
prefix <- samples[3]
print(prefix)
setwd(prefix)

# reading data
clonotypes <- read.csv("clonotypes.csv", header = TRUE, sep = ",")
contigs <- read.csv("filtered_contig_annotations.csv", header = TRUE, sep = ",")

# categorizing clonotypes (Ken's scheme)
clonotypes$category <- NA
for (i in 1:nrow(clonotypes)) {
  tra <- str_count(clonotypes[i, ]$cdr3s_aa, "TRA:") # count the number of alpha chains
  trb <- str_count(clonotypes[i, ]$cdr3s_aa, "TRB:") # count the number of beta chains
  if (tra == 1 & trb == 1) {
    clonotypes[i, ]$category <- "1.1A1B"
  }
  else if (tra == 2 & trb == 1) {
    clonotypes[i, ]$category <- "2.2A1B"
  }
  else if (tra == 1 & trb == 2) {
    clonotypes[i, ]$category <- "3.1A2B"
  }
  else if ((tra == 1 | tra == 2) & trb == 0) {
    clonotypes[i, ]$category <- "4.A0B"
  }
  else if (tra == 0 & (trb == 1 | trb == 2)) {
    clonotypes[i, ]$category <- "5.B0A"
  }
  else if (tra >= 2 & trb >= 2) {
    clonotypes[i, ]$category <- "6.multiple"
  }
  else {
    cat('no tcrs are detected\n')
    print(clonotypes[i, ]$cdr3s_aa)
  }
}

# count for each category
table(factor(clonotypes$category))

# getting unique cell barcodes
cells <- unique(contigs$barcode)

# getting clonotype per barcode
clonotypes.per.cell.aa <- c() # amino acids?
clonotypes.per.cell.nt <- c() # nucleotides?
clonotypes.per.cell.genes <- c() # the genes?
for (cell in cells) {
  ## pull out the contigs that poss quality thresholds
  data <- contigs[contigs$barcode == cell & contigs$high_confidence == "True" & 
                    contigs$full_length == "True" & contigs$productive == "True" & 
                    contigs$chain %in% c("TRA", "TRB") & contigs$cdr3 != "None",]
  
  if (nrow(data) != 0) {
    # sorting by chain and cdr3_nt
    data$cdr3_nt <- as.character(data$cdr3_nt)
    data <- data[order(data$cdr3_nt, decreasing = FALSE, method = "radix"),]
    data <- data[order(data$chain, decreasing = FALSE),]
    
    # getting chains in correct format
    aa <- c()
    nt <- c()
    genes <- c()
    for (i in 1:nrow(data)) {
      aa[i] <- paste(data[i,]$chain, data[i,]$cdr3, sep = ":") 
      nt[i] <- paste(data[i,]$chain, data[i,]$cdr3_nt, sep = ":")
      if (data[i,]$chain == "TRA") {
        genes[i] <- paste(data[i,]$chain, paste(data[i,]$v_gene, data[i,]$j_gene, sep = ","), sep = ":") 
      } else if (data[i,]$chain == "TRB") {
        genes[i] <- paste(data[i,]$chain, paste(data[i,]$v_gene, data[i,]$d_gene, data[i,]$j_gene, sep = ","), sep = ":") 
      }
    }
    if (length(aa) > 1) {
      tcr.aa <- paste(aa, collapse = ";")
      tcr.nt <- paste(nt, collapse = ";")
      tcr.genes <- paste(genes, collapse = ";")
    }
    else {
      tcr.aa <- aa[1]
      tcr.nt <- nt[1]
      tcr.genes <- genes[1]
    }
    
    # adding tcr for that cell
    clonotypes.per.cell.aa[cell] <- tcr.aa
    clonotypes.per.cell.nt[cell] <- tcr.nt
    clonotypes.per.cell.genes[cell] <- tcr.genes
    
  }
  else {
    tcr.aa <- NA
    tcr.nt <- NA
    tcr.genes <- NA
    
    # adding tcr for that cell
    clonotypes.per.cell.aa[cell] <- tcr.aa
    clonotypes.per.cell.nt[cell] <- tcr.nt
    clonotypes.per.cell.genes[cell] <- tcr.genes
  }
}

# remove cells with no tcr
clonotypes.per.cell.aa <- clonotypes.per.cell.aa[!is.na(clonotypes.per.cell.aa)]
clonotypes.per.cell.nt <- clonotypes.per.cell.nt[!is.na(clonotypes.per.cell.nt)]
clonotypes.per.cell.genes <- clonotypes.per.cell.genes[!is.na(clonotypes.per.cell.genes)]

# generating table
df <- data.frame(names(clonotypes.per.cell.aa), clonotypes.per.cell.aa, clonotypes.per.cell.nt, clonotypes.per.cell.genes, 
                 row.names = NULL)
colnames(df) <- c("cell barcode", "cdr3_aa", "cdr3_nt", "V(D)J genes")

# adding frequency and category
df$frequency <- NA
df$category <- NA
df$cdr3_aa <- as.character(df$cdr3_aa)
df$cdr3_nt <- as.character(df$cdr3_nt)
for (i in 1:nrow(df)) {
  df[i,]$frequency <- clonotypes[clonotypes$cdr3s_aa == df[i,]$cdr3_aa & clonotypes$cdr3s_nt == df[i,]$cdr3_nt,]$frequency
  df[i,]$category <- clonotypes[clonotypes$cdr3s_aa == df[i,]$cdr3_aa & clonotypes$cdr3s_nt == df[i,]$cdr3_nt,]$category
}

# category stats
categories <- as.data.frame(table(factor(df$category)))
colnames(categories) <- c("category", "frequency")

# writing output file
setwd("../")
wb <- createWorkbook()
name <- "TCR clones"
addWorksheet(wb, name)
writeData(wb, sheet = name, df)

name <- "categories"
addWorksheet(wb, name)
writeData(wb, sheet = name, categories)
saveWorkbook(wb, file = paste(prefix, "TCRs.xlsx", sep = "-"), overwrite = TRUE)

