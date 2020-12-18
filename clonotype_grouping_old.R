library(tidyverse)
library(openxlsx)
library(Biostrings)


LoadData <- function(filelabels, samplelabels, patient) {
  # p2 = readRDS('R_output/p2.seurat.01282020.rds')
  show(filelabels[1])
  data <- read.xlsx(filelabels[1])
  data$origin <- samplelabels[1]
  data$patient <- patient
  
  if (length(filelabels) > 1) {
    for (i in 2:length(filelabels)) {
      show(filelabels[i])
      temp <- read.xlsx(filelabels[i])
      sub = paste0('-', i)
      temp$cell.barcode <- gsub("-1", sub, temp$cell.barcode)
      temp$origin <- samplelabels[i]
      temp$patient<- patient
      show('dimensions of excel file:')
      show(nrow(temp))
      data = rbind(data, temp)
    }
  }
  data <- data[!data$category %in% c("6.multiple", NA), ]
  
  return(data)
}


SplitData <- function(data) {
  show(nrow(data))
  show('in split data')
  data = data[order(-data$frequency),]
  
  data$alpha.1 <- 0
  data$alpha.2 <- 0
  data$beta.1 <- 0
  data$beta.2 <- 0
  
  for (i in 1:nrow(data)) {
    if (i %% 1000 == 0) {
      show(i)
    }
    tcr <- strsplit(data[i, ]$cdr3_aa, ";")[[1]]
    alphas <- grep("^TRA", tcr, value = TRUE)
    betas <- grep("^TRB", tcr, value = TRUE)
    
    genes <- strsplit(data[i, ]$`V(D)J.genes`, ";")[[1]]
    alpha.genes <- grep("TRA", genes, value = TRUE)
    beta.genes <- grep("TRB", genes, value = TRUE)
    
    # getting alpha chain(s)
    if (length(alphas) == 1) {
      alpha <- paste(strsplit(alpha.genes, ":")[[1]][2], strsplit(alphas, ":")[[1]][2], sep = ",")
      data[i, ]$alpha.1 <- alpha
    } else if (length(alphas) == 2) {
      m <- matchPattern(alphas[1], alphas[2], max.mismatch = 2, min.mismatch = 0, with.indels = FALSE, fixed = TRUE, algorithm = "auto")
      if (length(m@ranges@start) != 0) {
        alpha <- paste(strsplit(alpha.genes[1], ":")[[1]][2], strsplit(alphas[1], ":")[[1]][2], sep = ",")
        data[i, ]$alpha.1 <- alpha
      }
      else {
        data[i, ]$alpha.1 <- paste(strsplit(alpha.genes[1], ":")[[1]][2], strsplit(alphas[1], ":")[[1]][2], sep = ",")
        data[i, ]$alpha.2 <- paste(strsplit(alpha.genes[2], ":")[[1]][2], strsplit(alphas[2], ":")[[1]][2], sep = ",")
      }
    }
    
    # getting beta chain(s)
    if (length(betas) == 1) {
      beta <- paste(strsplit(beta.genes, ":")[[1]][2], strsplit(betas, ":")[[1]][2], sep = ",")
      data[i, ]$beta.1 <- beta
    } else if (length(betas) == 2) {
      m <- matchPattern(betas[1], betas[2], max.mismatch = 2, min.mismatch = 0, with.indels = FALSE, fixed = TRUE, algorithm = "auto")
      if (length(m@ranges@start) != 0) {
        beta <- paste(strsplit(beta.genes[1], ":")[[1]][2], strsplit(betas[1], ":")[[1]][2], sep = ",")
        data[i, ]$beta.1 <- beta
      }
      else {
        data[i, ]$beta.1 <- paste(strsplit(beta.genes[2], ":")[[1]][2], strsplit(betas[2], ":")[[1]][2], sep = ",")
        data[i, ]$beta.2 <- paste(strsplit(beta.genes[2], ":")[[1]][2], strsplit(betas[2], ":")[[1]][2], sep = ",")
      }
    }
  }
  return(data)
}


ReorderAlphasBetas <- function(data) {
  #### re-classify any with an alpha or beta that has been repeated.
  a1 = 0
  b1 = 0
  for (i in 1:nrow(data)) {
    if (i %% 1000 == 0) {
      show(i)
    }
    alpha.1 = data[i,]$alpha.1
    alpha.2 = data[i,]$alpha.2
    beta.1 = data[i,]$beta.1
    beta.2 = data[i,]$beta.2
    
    if (alpha.1 != '0' & alpha.1==alpha.2) { data[i,]$alpha.2<-'0'; a1 = a1+1; data[i,]$category<-'1.1A1B'; }
    if (beta.1 != '0' & beta.1==beta.2) { data[i,]$beta.2<-'0'; b1 = b1+1; data[i,]$category<-'1.1A1B'; }
    
    ### MAKING SURE THE SEQUENCES ARE IN THE CORRECT ORDER
    if (alpha.1!='0' & alpha.2!='0' & alpha.1 < alpha.2) { 
      tmp = alpha.1
      data[i,]$alpha.1<-alpha.2
      data[i,]$alpha.2<-tmp
    }
    if (beta.1!='0' & beta.2!='0' & beta.1 < beta.2) {
      tmp = beta.1
      data[i,]$beta.1<-beta.2
      data[i,]$beta.2<-tmp
    }
    
  }
  return(data)
}



PrepTcrData <- function(data) {
  require(tidyverse)
  data = data %>% separate(alpha.1, ',', into=c('TRAV', 'TRAJ', 'alpha.aa'), remove=F)
  data = data %>% separate(beta.1, ',', into=c('TRBV', 'TRBD', 'TRBJ', 'beta.aa'), remove=F)
  data = data %>% separate(alpha.2, ',', into=c('TRAV.2', 'TRAJ.2', 'alpha.aa.2'), remove=F)
  data = data %>% separate(beta.2, ',', into=c('TRBV.2', 'TRBD.2', 'TRBJ.2', 'beta.aa.2'), remove=F)
  
  to_remove = c('TRAV', 'TRAJ', 'alpha.aa', 'TRBV', 'TRBD', 'TRBJ', 'beta.aa', 'TRAV.2', 'TRAJ.2', 'alpha.aa.2', 'TRBV.2', 'TRBD.2', 'TRBJ.2', 'beta.aa.2')
  
  data$clone <- paste(data$cdr3_aa, data$`V(D)J.genes`, sep = "-")
  clones <- unique(data$clone)
  
  data = data[order(-data$frequency), ]
  data = data[order(data$category), ]
  
  clones <- data[,c('clone', 'category')];  clones$found<-FALSE
  clones <- unique(clones)
  clones = rbind(clones[clones$category=='2.2A1B',], clones[clones$category!='2.2A1B', ])
  
  
  all.tcrgroups = data.frame(matrix(nrow=0, ncol=(ncol(data)+1)))
  colnames(all.tcrgroups) = c(colnames(data), 'clonotype')
  return(list(data, clones, all.tcrgroups))
}



GroupClonotypeFamilies <- function(allsplit, clones, all.tcrgroups) {
  rest = allsplit
  rest$found <-FALSE
  clono.type = 1
  small = 0
  match_rules = c(1,0,0,1)
  show(nrow(clones))
  cont = TRUE
  
  for (i in 1:nrow(clones)) {
    if (i %% 500 == 0){
      cat('row', i, '\n')
    }
    
    cl = clones[i,]$clone
    
    if (cont == TRUE) {
      
      if ( clones[i,]$found==FALSE ) {
        ## PULL OUT ANY IDENTICAL MATCHES, REMOVE FROM LARGER TABLE
        rest[rest$clone==cl, ]$found <- TRUE
        exact.match = rest[rest$clone==cl, ]
        
        if (nrow(exact.match)>0){
          v.alpha <- exact.match[1,]$TRAV
          j.alpha <- exact.match[1,]$TRAJ
          
          v.alpha.2 <- exact.match[1,]$TRAV.2
          j.alpha.2 <- exact.match[1,]$TRAJ.2
          
          v.beta <- exact.match[1,]$TRBV
          d.beta <- exact.match[1,]$TRBD
          j.beta <- exact.match[1,]$TRBJ
          
          v.beta.2 <- exact.match[1,]$TRBV.2
          d.beta.2 <- exact.match[1,]$TRBD.2
          j.beta.2 <- exact.match[1,]$TRBJ.2
          
          cdr3.a <- exact.match[1,]$alpha.aa
          cdr3.a.2 <- exact.match[1,]$alpha.aa.2
          cdr3.b <- exact.match[1,]$beta.aa
          cdr3.b.2 <- exact.match[1,]$beta.aa.2
          
          ## LOOK THROUGH THE REST OF THE CLONES TO FIND ANY IN THE SAME FAMILY
          ## category 1
          if (exact.match[1,]$category=='1.1A1B') {
            ## category 1: match exactly
            vdj.match = rest %>% subset(category=='1.1A1B' & TRAV==v.alpha & TRAJ==j.alpha & TRBV==v.beta & TRBJ==j.beta & found==FALSE)
            if (nrow(vdj.match)>0) {
              cdr3.match = vdj.match %>% subset(agrepl(cdr3.a, alpha.aa, max.distance = match_rules) & agrepl(cdr3.b, beta.aa, max.distance = match_rules))
              exact.match = rbind(exact.match, cdr3.match)
            }
            
            ## category 2: either alpha matches and the beta matches
            vdj.match = rest %>% subset(category=='2.2A1B' & ((TRAV==v.alpha & TRAJ==j.alpha) | (TRAV.2==v.alpha & TRAJ.2==j.alpha)) & TRBV==v.beta & TRBJ==j.beta & found==FALSE)
            if (nrow(vdj.match)>0) {
              cdr3.match = vdj.match %>% subset((agrepl(cdr3.a, alpha.aa, max.distance = match_rules) | agrepl(cdr3.a, alpha.aa.2, max.distance = match_rules)) & agrepl(cdr3.b, beta.aa, max.distance = match_rules))
              exact.match = rbind(exact.match, cdr3.match)
            }
            
            ## category 3: alpha matches and one of betas matches
            vdj.match = rest %>% subset(category=='3.1A2B' & (TRAV==v.alpha & TRAJ==j.alpha) & ((TRBV.2==v.beta & TRBJ.2==j.beta) | (TRBV==v.beta & TRBJ==j.beta)) & found==FALSE)
            if (nrow(vdj.match)>0) {
              cdr3.match = vdj.match %>% subset(agrepl(cdr3.a, alpha.aa, max.distance = match_rules) & (agrepl(cdr3.b, beta.aa.2, max.distance = match_rules) | agrepl(cdr3.b, beta.aa, max.distance = match_rules)))
              exact.match = rbind(exact.match, cdr3.match)
            }
            ## category 4: match alpha but no beta
            vdj.match = rest %>% subset(category=='4.A0B' & TRAV==v.alpha & TRAJ==j.alpha & found==FALSE)
            if (nrow(vdj.match)>0) {
              cdr3.match = vdj.match %>% subset(agrepl(cdr3.a, alpha.aa, max.distance = match_rules))
              exact.match = rbind(exact.match, cdr3.match)
            }
            
            ## category 5: match beta but no alpha
            vdj.match = rest %>% subset(category=='5.B0A' & TRBV==v.beta & TRBJ==j.beta & found==FALSE)
            if (nrow(vdj.match)>0) {
              cdr3.match = vdj.match %>% subset(agrepl(cdr3.b, beta.aa, max.distance = match_rules))
              exact.match = rbind(exact.match, cdr3.match)
            }
            rest[rest$clone %in% exact.match$clone, ]$found<-TRUE
            clones[clones$clone %in% exact.match$clone, ]$found<-TRUE
            
          } else if (exact.match[1,]$category=='2.2A1B') {
            ## category 1: match exactly
            vdj.match = rest %>% subset(category=='2.2A1B' & TRAV==v.alpha & TRAJ==j.alpha & TRAV.2==v.alpha.2 & TRAJ.2==j.alpha.2 & TRBV==v.beta & TRBJ==j.beta & found==FALSE)
            if (nrow(vdj.match)>0) {
              cdr3.match = vdj.match %>% subset(agrepl(cdr3.a, alpha.aa, max.distance = match_rules) & agrepl(cdr3.a.2, alpha.aa.2, max.distance = match_rules) & agrepl(cdr3.b, beta.aa, max.distance = match_rules))
              exact.match = rbind(exact.match, cdr3.match)
            }
            
            ## category 2: first alpha and beta matches
            vdj.match = rest %>% subset(category=='1.1A1B' & TRAV==v.alpha & TRAJ==j.alpha & TRBV==v.beta & TRBJ==j.beta & found==FALSE)
            if (nrow(vdj.match)>0) {
              cdr3.match = vdj.match %>% subset(agrepl(cdr3.a, alpha.aa, max.distance = match_rules) & agrepl(cdr3.b, beta.aa, max.distance = match_rules))
              exact.match = rbind(exact.match, cdr3.match)
            }
            
            ## category 3: second alpha and beta matches
            vdj.match = rest %>% subset(category=='1.1A1B' & TRAV==v.alpha.2 & TRAJ==j.alpha.2 & TRBV==v.beta & TRBJ==j.beta & found==FALSE)
            if (nrow(vdj.match)>0) {
              cdr3.match = vdj.match %>% subset(agrepl(cdr3.a.2, alpha.aa, max.distance = match_rules) & agrepl(cdr3.b, beta.aa, max.distance = match_rules))
              exact.match = rbind(exact.match, cdr3.match)
            }
            
            ## category 4: an alpha matches but no beta
            vdj.match = rest %>% subset(category=='4.A0B' & ((TRAV==v.alpha & TRAJ==j.alpha) | (TRAV==v.alpha.2 & TRAJ==j.alpha.2)) & found==FALSE)
            if (nrow(vdj.match)>0) {
              cdr3.match = vdj.match %>% subset(agrepl(cdr3.a, alpha.aa, max.distance = match_rules) | agrepl(cdr3.a.2, alpha.aa, max.distance = match_rules))
              exact.match = rbind(exact.match, cdr3.match)
            }
            
            ## category 5: match beta but no alpha
            vdj.match = rest %>% subset(category=='5.B0A' & TRBV==v.beta & TRBJ==j.beta & found==FALSE)
            if (nrow(vdj.match)>0) {
              cdr3.match = vdj.match %>% subset(agrepl(cdr3.b, beta.aa, max.distance = match_rules))
              exact.match = rbind(exact.match, cdr3.match)
              # rest = rest[!(rest$clone %in% cdr3.match$clone), ]
            }
            
            # show(nrow(rest[rest$clone %in% clonotype.group$clone, ]))
            rest[rest$clone %in% exact.match$clone, ]$found<-TRUE
            clones[clones$clone %in% exact.match$clone, ]$found<-TRUE
            # clonotype.group$clonotype <- clono.type
            
          } else if (exact.match[1,]$category=='3.1A2B') {
            # show('3.1A2B')
            ## category 1: match exactly
            vdj.match = rest %>% subset(category=='3.1A2B' & TRAV==v.alpha & TRAJ==j.alpha & TRBV.2==v.beta.2 & TRBJ.2==j.beta.2 & TRBV==v.beta & TRBJ==j.beta & found==FALSE)
            cdr3.match = vdj.match %>% subset(agrepl(cdr3.a, alpha.aa, max.distance = match_rules) & agrepl(cdr3.b.2, beta.aa.2, max.distance = match_rules) & agrepl(cdr3.b, beta.aa, max.distance = c(1,1,1,1)))
            exact.match = rbind(exact.match, cdr3.match)
            
            ## category 2: first alpha and beta matches
            vdj.match = rest %>% subset(category=='1.1A1B' & TRAV==v.alpha & TRAJ==j.alpha & TRBV==v.beta & TRBJ==j.beta & found==FALSE)
            cdr3.match = vdj.match %>% subset(agrepl(cdr3.a, alpha.aa, max.distance = match_rules) & agrepl(cdr3.b, beta.aa, max.distance = match_rules))
            exact.match = rbind(exact.match, cdr3.match)
            
            ## category 3: first alpha and second beta matches
            vdj.match = rest %>% subset(category=='1.1A1B' & TRAV==v.alpha & TRAJ==j.alpha & TRBV==v.beta.2 & TRBJ==j.beta.2 & found==FALSE)
            cdr3.match = vdj.match %>% subset(agrepl(cdr3.a, alpha.aa, max.distance = match_rules) & agrepl(cdr3.b.2, beta.aa, max.distance = match_rules))
            exact.match = rbind(exact.match, cdr3.match)
            
            ## category 4: alpha matches but no betas
            vdj.match = rest %>% subset(category=='4.A0B' & TRAV==v.alpha & TRAJ==j.alpha & found==FALSE)
            cdr3.match = vdj.match %>% subset(agrepl(cdr3.a, alpha.aa, max.distance = match_rules))
            exact.match = rbind(exact.match, cdr3.match)
            
            ## category 5: either beta matches but no alpha
            vdj.match = rest %>% subset(category=='5.B0A' & ((TRBV==v.beta & TRBJ==j.beta) | (TRBV.2==v.beta.2 & TRBJ.2==j.beta.2)) & found==FALSE)
            cdr3.match = vdj.match %>% subset(agrepl(cdr3.b, beta.aa, max.distance = match_rules) | agrepl(cdr3.b.2, beta.aa.2, max.distance = match_rules))
            exact.match = rbind(exact.match, cdr3.match)
            rest[rest$clone %in% exact.match$clone, ]$found<-TRUE
            clones[clones$clone %in% exact.match$clone, ]$found<-TRUE
            
          } else if (exact.match[1,]$category=='4.A0B') {
            ## category 4: alpha matches but no betas
            vdj.match = rest %>% subset(category=='4.A0B' & TRAV==v.alpha & TRAJ==j.alpha & found==FALSE)
            if (nrow(vdj.match)>0) {
              cdr3.match = vdj.match %>% subset(agrepl(cdr3.a, alpha.aa, max.distance = match_rules))
              exact.match = rbind(exact.match, cdr3.match)
            }
            rest[rest$clone %in% exact.match$clone, ]$found<-TRUE
            clones[clones$clone %in% exact.match$clone, ]$found<-TRUE
            
          } else if (exact.match[1,]$category=='5.B0A') {
            ## category 4: alpha matches but no betas
            vdj.match = rest %>% subset(category=='5.B0A' & TRBV==v.beta & TRBJ==j.beta & found==FALSE)
            if (nrow(vdj.match)>0) {
              cdr3.match = vdj.match %>% subset(agrepl(cdr3.b, beta.aa, max.distance = match_rules))
              exact.match = rbind(exact.match, cdr3.match)
            }
            rest[rest$clone %in% exact.match$clone, ]$found<-TRUE
            clones[clones$clone %in% exact.match$clone, ]$found<-TRUE
          } # end category checking
        } 
        
        if (nrow(exact.match) < 3) {
          small = small+1
        }
        if (nrow(exact.match) > 0) {
          exact.match$clonotype <- clono.type
          clono.type = clono.type + 1
        } else {
          exact.match$clonotype <- 0
        }
        all.tcrgroups = rbind(all.tcrgroups, exact.match)
        clones[i,]$found<-TRUE
      }
    }
    if (nrow(clones)-i == nrow(subset(rest, rest$found==TRUE))) {
      cont = FALSE
      print('only singletons left')
      r = rest %>% subset(clone %in% clones[i:nrow(clones),]$clone)
      r$clonotype = clono.type:(clono.type+nrow(r)-1) 
      all.tcrgroups = rbind(all.tcrgroups, r)
    }
  }
  
  return(all.tcrgroups)
}

