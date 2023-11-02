library(dplyr)
# https://subtyping.geno2pheno.org

path <- "./"
alignm <- function(ref,seq) {
  seq0 <- gapidx <- score <- list()
  for (i in 1:length(seq)) {
    al <- msa::msa(c(ref,tolower(gsub(whitespace,"",seq[[i]]))),method="Muscle",type="dna")
    al <- t(apply(matrix(as.character(al@unmasked),ncol=1),1,
                  function(x) return(unlist(strsplit(x,"")))))
    ref <- tolower(al[1,])
    seq0[[i]] <- tolower(al[2,])
    gapidx[[i]] <- which(ref!="-")
    score[[i]] <- sum(ref==seq0[[i]])/sum(seq0[[i]]!="-")
  }
  return(list(seq=seq0,gapidx=gapidx,score=score))
}
# load los alamos hiv refs and convert to consensus
ref <- read.fasta(paste0(path,"ref.fasta"))
refs <- read.fasta(paste0(path,"HIV1_REF_2020_genome_DNA.fasta"))
mafft <- function(m,thread=1) {
  require(seqinr)
  names <- names(m)
  fname <- paste0("mafft_",runif(1),"_",as.numeric(Sys.time()),".fasta")
  write.fasta(m,names=names,file.out=fname)
  fname2 <- paste0("out_",fname)
  system(paste0("mafft --thread ",thread," ",fname," > ",fname2))
  n <- read.fasta(fname2)
  system(paste0("rm ",fname))
  system(paste0("rm ",fname2))
  return(n)
}
refs <- mafft(refs,thread=46)
refs <- do.call("rbind",refs)
ref <- apply(refs,2,function(x) {
  y <- table(x)
  y <- names(y)[which.max(y)]
  return(y)
})
refs <- refs[,ref!="-"]
ref <- ref[ref!="-"]
refs <- rbind(reference=ref,refs)
rownames(refs) <- gsub("^Ref\\.","",rownames(refs))
rownames(refs)[grep("A6_pol_reference",rownames(refs))] <- "A6.pol_reference"
refnames <- unlist(lapply(rownames(refs),function(x) {
  y <- unlist(strsplit(x,"\\."))[1]
}))
refnames[is.na(refnames)] <- "A6"
refnames[1] <- "reference"
refnames <- gsub(".*_","",refnames)
A6idx <- grep("A6",rownames(refs))
CURFidx <- grep("^[0-9]",rownames(refs))
CURFidx <- CURFidx[which(!CURFidx %in% grep("AE",rownames(refs)))]
exoten <- grep("^[N-P]$",refnames)
apes <- grep("GOR|CPZ",refnames)
refidx <- (2:nrow(refs))[which(!(2:nrow(refs)) %in% c(CURFidx,apes))]
bpidx <- 1:ncol(refs)
refspart <- refs[refidx,bpidx]
reference <- refs[1,bpidx]
refnames0 <- refnames[refidx]
# reduce
newrefs <- NULL
for (i in unique(refnames)) {
  print(i)
  idx <- which(refnames==i)
  newref <- unlist(apply(refs[idx,,drop=FALSE],2,function(x) {
    y <- x[x!="-"]
    y <- table(y)
    y <- names(y)[rev(order(y))[1]]
    if (all(x=="-")) y <- "-"
    return(y)
  }))
  newrefs <- rbind(newrefs,newref)
}
rownames(newrefs) <- unique(refnames)
newrefs <- newrefs[rownames(newrefs) %in% refnames0,]
# here we manually check if all references are assigned to their consensus reference
# now we start building the list with hiv and later add hev
newrefidx <- grep(paste("^",unique(refnames[refidx]),"$",sep="",collapse="|"),rownames(newrefs))
refspart <- newrefs[newrefidx,]
refnames0 <- rownames(refspart)
allrefs <- allref <- allrefnames <- allrefids <- list()
allrefs[["HIV"]] <- refspart
allref[["HIV"]] <- paste(reference,collapse="")
allrefnames[["HIV"]] <- refnames0
allrefids[["HIV"]] <- "msa consensus"
if (file.exists(paste0(path,"HIVrefs"))) {
  files <- list.files(paste0(path,"HIVrefs"))
  for (file in files) {
    tmp <- read.fasta(paste0(path,"HIVrefs/",file))
    for (i in 1:length(tmp)) {
      tmptype <- gsub("\\..*","",names(tmp)[i])
      tmpal <- alignm(allref[["HIV"]],paste(tmp[[i]],collapse=""))
      allrefs[["HIV"]] <- rbind(allrefs[["HIV"]],tmpal$seq[[1]][tmpal$gapidx[[1]]])
      rownames(allrefs[["HIV"]])[nrow(allrefs[["HIV"]])] <- tmptype
      allrefnames[["HIV"]] <- c(allrefnames[["HIV"]],tmptype)
    }
  }
}
# now hev
filename <- "hevrefs.rds"
if (file.exists(filename)) {
  refs <- readRDS(filename)
} else {
  refs <- NULL
  files <- list.files("hevtypes/") # this is the folder with references from hev-glue
  for (file in files) {
    subtype <- unlist(strsplit(file,"_"))[2]
    tmp <- read.fasta(paste0("hevtypes/",file))
    names(tmp) <- paste0(subtype,".",names(tmp))
    refs <- c(refs,tmp)
  }
  mafft <- function(m,thread=1) {
    require(seqinr)
    names <- names(m)
    fname <- paste0("mafft_",runif(1),"_",as.numeric(Sys.time()),".fasta")
    write.fasta(m,names=names,file.out=fname)
    fname2 <- paste0("out_",fname)
    system(paste0("mafft --thread ",thread," ",fname," > ",fname2))
    n <- read.fasta(fname2)
    system(paste0("rm ",fname))
    system(paste0("rm ",fname2))
    return(n)
  }
  refs <- mafft(refs,thread=46)
  refs <- do.call("rbind",refs)
  refs <- apply(refs,c(1,2),tolower)
  refs <- refs[-grep("4\\.",rownames(refs)),]
  ref <- apply(refs,2,function(x) {
    y <- table(x)
    y <- names(y)[which.max(y)]
    return(y)
  })
  refs <- refs[,ref!="-"]
  ref <- ref[ref!="-"]
  refs <- rbind(reference=ref,refs)
  saveRDS(refs,filename)
}
rownames(refs) <- gsub("\\..*","",rownames(refs))
newrefs <- NULL
for (i in unique(rownames(refs))) {
  print(i)
  idx <- which(rownames(refs)==i)
  newref <- unlist(apply(refs[idx,,drop=FALSE],2,function(x) {
    y <- x[x!="-"]
    y <- table(y)
    y <- names(y)[rev(order(y))[1]]
    if (all(x=="-")) y <- "-"
    return(y)
  }))
  newrefs <- rbind(newrefs,newref)
}
rownames(newrefs) <- unique(rownames(refs))
# again we check the quality of the consensus refs like with hiv
# last we add hev to the list
refs <- newrefs
bpidx <- 1:ncol(refs)
allrefs[["HEV"]] <- refs[,bpidx]
allref[["HEV"]] <- paste(allrefs[["HEV"]][1,],collapse="")
allrefs[["HEV"]] <- allrefs[["HEV"]][-1,]
allrefnames[["HEV"]] <- rownames(allrefs[["HEV"]])
allrefids[["HEV"]] <- "msa consensus based on http://hev.glue.cvr.ac.uk/"

# now we start the hiv analysis
refspart <- allrefs$HIV

# we convert reference names into a nicer format and fully resolve CRFs which include numbers
fixed <- lapply(1:nrow(refs),function(i) {
  y <- unlist(strsplit(rownames(refs)[i],"\\."))[2]
  y <- gsub(".*_|/","",y)
  z <- gsub("[A-Z].*|[a-z].*","",y)
  y <- gsub("[0-9]","",y)
  if (is.na(z)) z <- ""
  if (is.na(y)) y <- ""
  if (nchar(z)>1) {
    z <- apply(matrix(unlist(strsplit(z,"")),nrow=2),2,paste,collapse="")
    z <- unlist(lapply(z,function(u) {
      v <- unique(refnames[grep(paste0("\\.",u,"_"),rownames(refs))])
    }))
    z <- unlist(strsplit(gsub("[0-9]","",z),""))
  } else {
    z <- ""
  }
  y <- unlist(strsplit(gsub("[0-9]","",y),""))
  y <- unique(c(y,z))
  y <- y[y!=""]
  return(y)
})
# we performe prediction for a wide range of alphas
alphas <- seq(0.01,0.5,0.01)
for (alpha in alphas) {
  print(alpha)
  df0[[as.character(alpha)]] # <- data.frame with annotated c(subtype,predicted subtype,fit,breaks) in each row
}
# we compute the accuracies for each alpha
received <- lapply(1:length(fixed),function(i) { # annotated subtype
  x <- fixed[[i]]
  y <- sort(x)
  if (is.na(x[1])) {
    x <- refnames[i]
    y <- unlist(strsplit(x,""))
    y <- gsub("[0-9]","",y)
  }
  y <- y[y!=" " & y!="" & y!="/" & y!="U" & y!="c" & y!="p" & y!="x"]
  y <- unique(y)
  return(y)
})
noidx <- "cpx|GOR|CPZ"
auc <- matrix(NA,length(alphas),2)
colnames(auc) <- c("reference in prediction","prediction in reference")
rownames(auc) <- alphas
for (alpha in alphas) {
  df <- do.call("rbind",df0[[as.character(alpha)]])
  called <- lapply(df[,2],function(x) { # predicted subtype
    y <- unlist(strsplit(x,""))
    y <- gsub("[0-9]","",y)
    y <- sort(y)
    y <- y[y!=" " & y!="" & y!="/"]
    y <- unique(y)
    return(y)
  })
  acc <- matrix(NA,nrow(refs),2)
  for (i in 2:nrow(df)) {
    if (length(grep(noidx,rownames(refs)[i]))>0) next()
    acc[i,] <- c(sum(received[[i]] %in% called[[i]])/length(received[[i]]),
                 sum(called[[i]] %in% received[[i]])/length(called[[i]]))
  }
  colnames(acc) <- c("reference in prediction","prediction in reference")
  auc[as.character(alpha),] <- apply(acc,2,mean,na.rm=TRUE)
}

par(mfrow=c(1,1))
plot(alphas,auc[,1],type="l",col=2,lwd=3,ylim=c(0,1),ylab="accuracy",xlab=expression(alpha))
lines(alphas,auc[,2],type="l",col=4,lwd=3)
legend("bottomright",fill=c(2,4),legend=colnames(auc))
idx <- which(round(auc[,1],2)==max(round(auc[,1],2)))[1]
abline(v=alphas[idx],lty=2)
axis(3,alphas[idx],round(alphas[idx],4))
idx <- which(round(auc[,2],2)==max(round(auc[,2],2)))
idx <- idx[length(idx)]
abline(v=alphas[idx],lty=2)
axis(3,alphas[idx],round(alphas[idx],4))
idx <- which.min(abs(round(auc[,1],2)-round(auc[,2],2)))
abline(v=alphas[idx],lty=2)
axis(3,alphas[idx],round(alphas[idx],4))

# fix alpha and compute accuracy
alpha <- alpha0 <- round(alphas[idx],2)
keep <- 0
df <- do.call("rbind",df0[[as.character(alpha)]])
head(df)
called <- lapply(df[,2],function(x) {
  y <- unlist(strsplit(x,""))
  y <- gsub("[0-9]","",y)
  y <- sort(y)
  y <- y[y!=" " & y!="" & y!="/"]
  y <- unique(y)
  return(y)
})
acc <- matrix(NA,nrow(refs),2)
for (i in 2:nrow(df)) {
  if (length(grep(noidx,rownames(refs)[i]))>0) next()
  acc[i,] <- c(sum(received[[i]] %in% called[[i]])/length(received[[i]]),
               sum(called[[i]] %in% received[[i]])/length(called[[i]]))
}
colnames(acc) <- c("reference in prediction","prediction in reference")
summary(acc)

par(mfrow=c(1,1))
mnem::moreboxplot(acc[CURFidx,],xaxt="n",dcol=rev(rgb(c(0,1),0,c(1,0),0.5)),
                  main="accuracy",gmin=0.1,ylim=c(0,1))
axis(1,1:ncol(acc),colnames(acc),las=1,cex.axis=0.75)

refspart <- rbind(allrefs$HIV,refs[apes,]) # add apes

cometSet <- read.fasta("training_HIV-1_v6.fasta") # comet training set
names(cometSet) <- gsub("\\..*","",names(cometSet))

comets # <- g2p-recco result with list(subtype1=list(predicted subtype,fit,break points),subtype2...)
cometfits <- as.numeric(unlist(lapply(comets,function(x) return(x$fit))))
comets <- unlist(lapply(comets,function(x) return(paste(sort(names(x$subtype)),collapse=", "))))

# accuracy the same as before
fixedC <- lapply(1:length(cometSet),function(i) {
  y <- unlist(strsplit(names(cometSet)[i],"\\."))[1]
  y <- gsub(".*_|/","",y)
  z <- gsub("[A-Z].*|[a-z].*","",y)
  y <- gsub("[0-9]","",y)
  if (is.na(z)) z <- ""
  if (is.na(y)) y <- ""
  if (nchar(z)>1) {
    z <- apply(matrix(unlist(strsplit(z,"")),nrow=2),2,paste,collapse="")
    z <- unlist(lapply(z,function(u) {
      v <- unique(names(cometSet)[grep(paste0("^",u,"_"),names(cometSet))])
    }))
    z <- unlist(strsplit(gsub("[0-9]","",z),""))
  } else {
    z <- ""
  }
  y <- unlist(strsplit(gsub("[0-9]","",y),""))
  y <- unique(c(y,z))
  y <- y[y!="" & y!="_"]
  return(y)
})
receivedC <- lapply(1:length(fixedC),function(i) {
  x <- fixedC[[i]]
  y <- sort(x)
  if (is.na(x[1])) {
    x <- names(cometSet)[i]
    y <- unlist(strsplit(x,""))
    y <- gsub("[0-9]","",y)
  }
  y <- y[y!=" " & y!="" & y!="/" & y!="U" & y!="c" & y!="p" & y!="x"]
  y <- unique(y)
  return(y)
})
calledC <- lapply(comets,function(x) {
  y <- unlist(strsplit(x,""))
  y <- gsub("[0-9]","",y)
  y <- sort(y)
  y <- y[y!=" " & y!="" & y!="/" & y!=","]
  y <- unique(y)
  return(y)
})
acc <- matrix(NA,length(calledC),2)
for (i in 1:length(comets)) {
  acc[i,] <- c(sum(receivedC[[i]] %in% calledC[[i]])/length(receivedC[[i]]),
               sum(calledC[[i]] %in% receivedC[[i]])/length(calledC[[i]]))
}
colnames(acc) <- c("reference in prediction","prediction in reference")
summary(acc[-grep(noidx,names(cometSet)),])
CURFidxC <- (1:nrow(acc))[-grep(noidx,names(cometSet))]

par(mfrow=c(1,1))
mnem::moreboxplot(acc[CURFidxC,],xaxt="n",dcol=rev(rgb(c(0,1),0,c(1,0),0.5)),
                  main="accuracy",gmin=0.1,ylim=c(0,1))
axis(1,1:ncol(acc),colnames(acc),las=1,cex.axis=0.75)

# synthetic analysis
# we construct a set of synthetic recombinations
acc2 <- fit <- dfsyn <- breaks <- list()

noise <- 0.5
keep <- 0
synrefs <- synreceived <- synidx <- synbreaks <- list()
for (i in 1:100) { # 100 sequence on subtyping.g2p.org takes approx 256 seconds, comet takes 18 seconds
  n <- sample(2:5,1)
  ns <- c(1,sort(sample(1:ncol(refs),n-1)),ncol(refs))
  while(any(diff(ns)/ncol(refs)<=keep)) {
    ns <- c(1,sort(sample(1:ncol(refs),n-1)),ncol(refs))
  }
  # while(all(diff(ns)/ncol(refs)>0.05) & keep==0) {
  #   ns <- c(1,sort(sample(1:ncol(refs),n-1)),ncol(refs))
  # } # this does not really influence the result much
  idx <- sample(refidx[!refidx %in% apes],n)
  newseq <- refs[1,]
  for (j in 1:length(idx)) {
    newseq[ns[j]:ns[j+1]] <- refs[idx[j],ns[j]:ns[j+1]]
  }
  newseq[sample(1:length(newseq),
                floor(noise*length(newseq)))] <- sample(c("a","c","g","t"),
                                                        floor(noise*length(newseq)),
                                                        replace=TRUE)
  synrefs[[i]] <- newseq
  synreceived[[i]] <- sort(unique(unlist(received[idx])))
  synidx[[i]] <- idx
  synbreaks[[i]] <- ns[-c(1,length(ns))]
}
write.fasta(synrefs,1:length(synrefs),file.out="syntest.fasta")

dfsyn[[as.character(noise)]] # <- g2-recco with c(synthetic subtype annotation,recco subtype,fit,annotated break points, predicted break points)

# as before we compute accuracy
df <- do.call("rbind",dfsyn[[as.character(noise)]])
called <- lapply(df[,2],function(x) {
  y <- unlist(strsplit(x,""))
  y <- gsub("[0-9]","",y)
  y <- sort(y)
  y <- y[y!=" " & y!="" & y!="/"]
  y <- unique(y)
  return(y)
})

acc <- matrix(NA,length(synrefs),2)
for (i in 1:nrow(df)) {
  acc[i,] <- c(sum(synreceived[[i]] %in% called[[i]])/length(synreceived[[i]]),
               sum(called[[i]] %in% synreceived[[i]])/length(called[[i]]))
}
colnames(acc) <- c("reference in prediction","prediction in reference")
mnem::moreboxplot(acc,xaxt="n",dcol=rev(rgb(c(0,1),0,c(1,0),0.5)),ylim=c(0,1))
summary(acc)
summary(as.numeric(df[,3]))
if (keep==0) {
  fit[[as.character(noise)]] <- as.numeric(df[,3])
} else {
  fit[[paste0(as.character(noise),"_",keep)]] <- as.numeric(df[,3])
}
breaks[[as.character(noise)]] <- df[,4:5]
head(df)

comet <- read.table("results.csv",fill=TRUE,header=TRUE,sep="\t")
comet <- comet[naturalsort::naturalorder(as.numeric(comet$name)),]
cometcalled <- lapply(comet$subtype,function(x) {
  x <- gsub("CPZ|GOR","",x) # we filter GOR and CPZ subtype annotation == accuracy of 0
  y <- unlist(strsplit(x,""))
  y <- gsub("[0-9]","",y)
  y <- sort(y)
  y <- y[y!=" " & y!="" & y!="/"]
  y <- y[!y %in% letters & !y %in% c("-",",","_","U","(",")","\"")]
  y <- unique(y)
  return(y)
})
accomet <- matrix(NA,length(synrefs),2)
for (i in 1:nrow(df)) {
  accomet[i,] <- c(sum(synreceived[[i]] %in% cometcalled[[i]])/length(synreceived[[i]]),
                   sum(cometcalled[[i]] %in% synreceived[[i]])/length(cometcalled[[i]]))
}
colnames(accomet) <- c("reference in prediction","prediction in reference")
mnem::moreboxplot(accomet,xaxt="n",dcol=rev(rgb(c(0,1),0,c(1,0),0.5)))
summary(accomet)

stan <- read.csv("sequenceSummaries.csv")
stan <- stan$Subtype....
stancalled <- lapply(stan,function(x) {
  x <- gsub("CRF|CPZ|GOR|Unknown","",x)
  x <- gsub("\\(.*","",x)
  y <- unlist(strsplit(x,""))
  y <- gsub("[0-9]","",y)
  y <- sort(y)
  y <- y[y!=" " & y!="" & y!="/"]
  y <- y[!y %in% letters & !y %in% c("-",",","_","U","(",")","\"")]
  y <- unique(y)
  return(y)
})
acstan <- matrix(NA,length(synrefs),2)
for (i in 1:nrow(df)) {
  acstan[i,] <- c(sum(synreceived[[i]] %in% stancalled[[i]])/length(synreceived[[i]]),
                  sum(stancalled[[i]] %in% synreceived[[i]])/length(stancalled[[i]]))
}
colnames(acstan) <- c("reference in prediction","prediction in reference")
mnem::moreboxplot(acstan,xaxt="n",dcol=rev(rgb(c(0,1),0,c(1,0),0.5)))
summary(acstan)

# break point accuracy
breaks2 <- list()
for (i in 1:length(breaks)) {
  breaks2[[i]] <- numeric()#numeric(nrow(breaks[[i]]))
  for (j in 1:nrow(breaks[[i]])) {
    a <- as.numeric(unlist(strsplit(breaks[[i]][j,1],",")))
    b <- as.numeric(unlist(strsplit(breaks[[i]][j,2],",")))
    breakmat <- as.matrix(dist(c(a,b)))
    diag(breakmat) <- Inf
    tmp <- apply(breakmat[,1:length(a),drop=FALSE],2,min)
    #if (any(is.infinite(tmp))) tmp <- unlist(lapply(a,function(x) return(min(x,ncol(refspart)-x))))
    breaks2[[i]] <- c(breaks2[[i]],tmp)
  }
  breaks2[[i]][is.infinite(breaks2[[i]])] <- ncol(refspart)/2
  breaks2[[i]] <- log10(breaks2[[i]]+1)
}

dcol <- c(rev(rgb(c(0,1),0,c(1,0),0.5)),rev(rgb(c(0,1),c(1,0),c(0,1),0.5)))
laymat <- matrix(c(1,4,2,5,3,5),2,3)
layout(laymat)
par(mfrow=c(2,3))
for (i in order(names(acc2))) {
  tmpname <- ifelse(length(grep("_",names(acc2)[i]))>0,paste0(gsub("_.*","",names(acc2)[i])," (<10%)"),names(acc2)[i])
  mnem::moreboxplot(acc2[[i]],xaxt="n",dcol=rep(dcol[1:2],ncol(acc2[[i]])/2),
                    main=paste0("noise: ",tmpname),
                    gcol=rgb(0,0,0,0),ylim=c(0,1),boxwex=0.5,ylab="accuracy")
  if (i==4) {
    legend("topright",fill=dcol[1:2],c("reference in prediction","prediction in reference"))
  }
  axis(1,1:ncol(acc2[[i]]),rep(c("G2P","COMET","STA"),each=2),las=2)
  lines(1:ncol(acc2[[i]]),apply(acc2[[i]],2,mean,na.rm=TRUE),col="white",type="p",pch=4)
}
cols <- rgb(0,1,0,0.5)
hist(fit[[1]],xlim=c(0.45,1),col=rgb(0,1,0,0.5),main="fit of query to prediction",
     xlab="fit",breaks=seq(0,1,0.01),ylim=c(0,70))
for (i in 2:(length(fit)-0)) {
  hist(fit[[i]],col=rgb(as.numeric(names(fit)[i])*2,1-as.numeric(names(fit)[i])*2,0,0.5),
       add=TRUE,breaks=seq(0,1,0.01))
  cols <- c(cols,rgb(as.numeric(names(fit)[i])*2,1-as.numeric(names(fit)[i])*2,0,0.5))
}
legend("topleft",fill=c(rgb(0,0,0,0),cols),c("noise:",sort(names(fit)[-5])),border=0)
mnem::moreboxplot(breaks2,dcol=cols,gcol=rgb(0,0,0,0),
                  xlab="noise",ylab="mean distance in log to base 10",yaxt="none")
print(10^unlist(lapply(breaks2,median)))
print(10^unlist(lapply(breaks2,mean)))
axis(2,0:100,10^(0:100))
lines(1:length(breaks2),unlist(lapply(breaks2,mean,na.rm=TRUE)),col="white",type="p",pch=4)

files <- paste0("A6seqs/",list.files("A6seqs/")) # folder with fastas from genbank
A6gb <- NULL
for (file in files) {
  tmp <- read.fasta(file)
  tmp <- tmp[grep("pol",unlist(lapply(tmp,attr,"Annot")))]
  A6gb <- c(A6gb,tmp)
}
A6 <- c(A6gb,A6isr,A6h) # we add not published A6 from Tel Aviv and Hamburg
isridx <- (length(A6gb)+1):(length(A6gb)+length(A6isr))
hidx <- (length(A6gb)+length(A6isr)+1):(length(A6gb)+length(A6isr)+length(A6h))
A6care # <- list(subtypeA6_1=list(predicted subtype,fit,aligned sequence))
A6care2 <- unlist(lapply(A6care,function(x) return(paste(sort(names(x$subtype)),collapse=", "))))
names(A6care2) <- names(A6)

par(mfrow=c(1,2))
bars <- log10(sort(table(A6care2)))
bars[bars==0] <- 0.01
barplot(bars,border=rgb(0,0,0,0),col=rgb(0,0,1,0.5),las=2,horiz=TRUE,xaxt="n",
        xlab="number of sequences (log10)",cex.names=0.8)
ticks <- c(1,10,100,700)
axis(1,log10(ticks),ticks)
A6fit <- as.numeric(unlist(lapply(A6care,function(x) return(x$fit))))
hist(A6fit,border=rgb(0,0,0,0),col=rgb(0,0,1,0.5),main="",xlab="fit to reference(s)")

par(mfrow=c(1,2))
bars <- sort(table(A6care2[isridx]))
bars[bars==0] <- 0.01
barplot(bars,border=rgb(0,0,0,0),col=rgb(0,0,1,0.5),las=2,horiz=TRUE,xaxt="n",
        xlab="number of sequences",cex.names=0.8)
ticks <- 1:10
axis(1,ticks,ticks)
A6fit <- as.numeric(unlist(lapply(A6care[isridx],function(x) return(x$fit))))
hist(A6fit,border=rgb(0,0,0,0),col=rgb(0,0,1,0.5),main="",xlab="fit to reference(s)")

par(mfrow=c(1,2))
bars <- sort(table(A6care2[hidx]))
bars[bars==0] <- 0.01
barplot(bars,border=rgb(0,0,0,0),col=rgb(0,0,1,0.5),las=2,horiz=TRUE,xaxt="n",
        xlab="number of sequences",cex.names=0.8)
ticks <- 1:100
axis(1,ticks,ticks)
A6fit <- as.numeric(unlist(lapply(A6care[hidx],function(x) return(x$fit))))
hist(A6fit,border=rgb(0,0,0,0),col=rgb(0,0,1,0.5),main="",xlab="fit to reference(s)")

alphas <- seq(0.01,1,0.01)

breakmat # this is matrix with row i and column j the number of recombination for alpha i and A6 sequence j
recIdx <- grep(",",A6care2)

par(mfrow=c(1,2))
plot(alphas,apply(breakmat[-recIdx,],2,median),type="l",
     ylim=c(0,max(quantile(breakmat[-recIdx,],0.975),quantile(breakmat[recIdx,],0.975))+10),col=rgb(0,0,1,1),
     ylab="breaks",xlab=expression(alpha),lwd=3)
polygon(c(alphas,rev(alphas)),c(apply(breakmat[-recIdx,],2,quantile,0.975),rev(apply(breakmat[-recIdx,],2,quantile,0.025))),
        col=rgb(0,0,1,0.5),border=rgb(0,0,0,0))
lines(alphas,apply(breakmat[recIdx,],2,median),col=rgb(1,0,0,1),lwd=3)
polygon(c(alphas,rev(alphas)),c(apply(breakmat[recIdx,],2,quantile,0.975),rev(apply(breakmat[recIdx,],2,quantile,0.025))),
        col=rgb(1,0,0,0.5),border=rgb(0,0,0,0))
legend("topleft",legend=c("recombination","no recombination"),
       border=rgb(0,0,0,0),fill=c(rgb(1,0,0,1),rgb(0,0,1,1)))
ymax <- 8
xmax <- 0.25
lwd <- 1
lty <- 2
segments(0,-2,0,ymax,lwd=lwd,lty=lty)
segments(0,ymax,xmax,ymax,lwd=lwd,lty=lty)
segments(xmax,ymax,xmax,-2,lwd=lwd,lty=lty)
segments(xmax,-2,0,-2,lwd=lwd,lty=lty)
plot(alphas,apply(breakmat[-recIdx,],2,median),type="l",ylim=c(0,ymax),col=rgb(0,0,1,1),
     ylab="breaks",xlab=expression(alpha),lwd=3,xlim=c(0,xmax))
polygon(c(alphas,rev(alphas)),c(apply(breakmat[-recIdx,],2,quantile,0.975),rev(apply(breakmat[-recIdx,],2,quantile,0.025))),
        col=rgb(0,0,1,0.5),border=rgb(0,0,0,0))
lines(alphas,apply(breakmat[recIdx,],2,median),col=rgb(1,0,0,1),lwd=3)
polygon(c(alphas,rev(alphas)),c(apply(breakmat[recIdx,],2,quantile,0.975),rev(apply(breakmat[recIdx,],2,quantile,0.025))),
        col=rgb(1,0,0,0.5),border=rgb(0,0,0,0))
blues <- apply(breakmat[-recIdx,],2,min)
reds <- apply(breakmat[recIdx,],2,min)
niceidx <- alphas[which(reds>0 & blues<1)[1]]
blues <- apply(breakmat[-recIdx,],2,quantile,0.975)
niceidx <- c(niceidx,alphas[which(reds==blues & reds>0)[1]])
reds <- apply(breakmat[recIdx,],2,median)
niceidx <- c(niceidx,alphas[which(reds<blues & reds>0)[1]])
abline(v=niceidx[1:2],lty=2)
axis(3,niceidx[1:2],niceidx[1:2])
abline(v=niceidx[3],lty=2)
axis(3,niceidx[3],"")
axis(3,niceidx[3],niceidx[3],line=1,tick=FALSE)

# HEV analysis

filename <- "hevrefs.rds"
refs <- readRDS(filename)
rownames(refs) <- gsub("\\..*","",rownames(refs))
unident <- read.fasta("AL_4_AB369688_sequences.fasta")
unlist(lapply(unident,function(x)return(length(x))))

refs <- allrefs$HEV
ref <- allref$HEV

Hev4 # <- list of lists, one for each sequence in unident with list(subtype="vector of subtypes",fit=fit))
Hev4types <- unlist(lapply(Hev4,function(x) return(paste(sort(names(x$subtype)),collapse="_"))))
Hev4fit <- unlist(lapply(Hev4,function(x) return(x$fit)))

alphas <- seq(0.01,0.99,0.01)
costs <- recombs <- breaks <- numeric(length(alphas))
for (i in 1:length(alphas)) {
  tmp # contains internal g2p-recco analysis
  costs[i] <- tmp$cost # cost value
  recombs[i] <- length(unique(tmp$z)) # number of subtypes
  breaks[i] <- sum(diff(tmp$z)!=0) # number of breaks
}

par(mfrow=c(1,1))
plot(alphas,costs,type="l",xaxt="n",xlab=expression(alpha),ylab="cost",col=rgb(0,0,1,0.5),lwd=3)
axis(1,seq(0,1,0.05),seq(0,1,0.05))
recidx <- which(c(0,diff(breaks))!=0)
recn <- breaks[which(c(0,diff(breaks))!=0)]
abline(v=alphas[recidx],col=rgb(recn/max(recn),1-recn/max(recn),0,0.5))
axis(3,alphas[recidx],breaks[which(c(0,diff(breaks))!=0)])
axisidx <- c(2,4:6,8,9,22)
axis(3,alphas[which(c(0,diff(breaks))!=0)],recombs[which(c(0,diff(breaks))!=0)],line=1,tick=FALSE,lwd=0)

tmp1 <- breaks
tmp2 <- recombs
tmp2 <- tmp2/max(tmp2)*max(tmp1)
lines(alphas,tmp2,col=rgb(0,0,1,0.5),lwd=3)
axis(4,seq(0,max(tmp2),length.out=10),round(seq(0,max(recombs),length.out=10)))

library(countrycode)
place <- as.data.frame(countrycode::codelist)[,c("country.name.en","eurocontrol_pru","continent","region","region23")]
hevglue <- read.fasta("HEV_all.fasta") # all hev-glue sequences
hevgluemeta <- read.delim("sequences/HEV_all.txt") # meta data
hevtypes # <- same type of list object as for the unidentified hev sequences

table(place$continent)
Hevcont <- list()
for (cont in unique(place$continent)) {
  if (is.na(cont)) next()
  cidx <- which(hevgluemeta$gb_country_short %in% place$country.name.en[place$continent==cont])
  metas <- hevgluemeta[cidx,]
  types <- do.call("rbind",lapply(cidx,function(i) {
    y <- hevtypes[[i]]
    y <- c(paste(names(y$subtype),collapse=","),y$fit)
    return(y)
  }))
  Hevcont[[cont]] <- c(length(grep(",",types[,1]))/nrow(types),
                       median(as.numeric(types[grep(",",types[,1]),2])),
                       median(as.numeric(types[,2])),
                       length(grep(",",types[,1])),length(cidx))
  if (is.na(Hevcont[[cont]][1])) {
    Hevcont[[cont]] <- NULL
  } else {
    cat(paste(c(cont,round(Hevcont[[cont]],2)),collapse=" & "))
    cat(" \\\\\n")
  }
}
Hevcont <- do.call("rbind",Hevcont)
Ftest <- matrix(NA,nrow(Hevcont),nrow(Hevcont))
colnames(Ftest) <- rownames(Ftest) <- rownames(Hevcont)
for (i in 1:nrow(Hevcont)) {
  for (j in 1:nrow(Hevcont)) {
    Ftest[i,j] <- fisher.test(matrix(c(Hevcont[i,4],Hevcont[i,5]-Hevcont[i,4],
                                       Hevcont[j,4],Hevcont[j,5]-Hevcont[j,4]),2,2))$p.value
  }
}
Ftest
