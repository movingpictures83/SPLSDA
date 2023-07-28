dyn.load(paste("RPluMA", .Platform$dynlib.ext, sep=""))
source("RPluMA.R")

library(mixOmics)
library(RCurl)
library(bitops)
library(DiscriMiner)

input <- function(inputfile) {
  pfix = prefix()
  if (length(pfix) != 0) {
     pfix <- paste(pfix, "/", sep="")
  }
  print(pfix)
  parameters <- read.table(inputfile, as.is=T);
  rownames(parameters) <- parameters[,1]
  print(rownames(parameters))
  
  metabolime <- read.csv(paste(pfix, toString(parameters["samples",2]), sep=""), header = F)
  obs_names=as.matrix(read.table(paste(pfix, toString(parameters["categories", 2]), sep="")))
  var_ID<<-as.matrix(read.table(paste(pfix, toString(parameters["observables", 2]), sep="")))
  #colnames(metabolime) <- obs_names
  #print(nrow(metabolime))
  #print(nrow(var_ID))
  rownames(metabolime) <- var_ID
  metabolime <- as.data.frame(t(metabolime))
  metabolime$obs_names<-as.factor(obs_names)


  #names=levels(metabolime$obs_names)
  #ind_total=c()
  #for (j in 1:(length(names)-1)) {
  #
  #  sub_met=metabolime[which(metabolime$obs_names==names[j]),]
  #  print(sub_met)
  #  x <- readline("Any key to continue")
  #  col_var=c()
  #  for(i in 1:(dim(sub_met)[2]-1)){
  #    col_var[i]=var(sub_met[,i])
  #  }
  #  ind=which(col_var<=0.01)
  #  ind_total=union(ind_total,ind)
  #}
  #metabolime=metabolime[,-ind_total]
  testers = read.delim(paste(pfix, toString(parameters["targets", 2]), sep=""), header=F, sep='\n', as.is=T)
  testindex=c()
  for (i in 1:length(testers[,1])) {
     testindex = c(testindex, which(metabolime$obs_names==testers[i,1]))
  }
  print(testindex)
  test = metabolime[testindex,]
  #print(test)
  X <<- subset(test,select = -obs_names)
  #print(X)
  Y <<- as.character(test$obs_names)
}

run <- function() {
   plsda.metabolite <<- splsda(X, Y, ncomp = 2)
   #print(plsda.metabolite)
   #my_pls1 <<- plsDA(X, Y, autosel=FALSE, comps=2)
   #print(levels(my_pls1$classification))
   #VIP <<- names(which(my_pls1$VIP[,2]>1))
}

output <- function(outputfile) {
   #write.csv(VIP, paste(outputfile, ".VIP.csv", sep=""))
   #x <- my_pls1$functions
   #print(nrow(x))
   #print(length(c("INTERCEPT", var_ID)))
   #rownames(x) <- c("INTERCEPT", colnames(X))
   #colnames(x) <- levels(my_pls1$classification)
   #write.csv(x, paste(outputfile, ".functions.csv", sep=""))
   #print(my_pls1$scores)
   #print(str(my_pls1$scores))
   #y <- my_pls1$scores
   #colnames(y) <- levels(my_pls1$classification)
   #write.csv(y, paste(outputfile, ".scores.csv", sep=""))
   #write.csv(my_pls1$scores, paste(outputfile, ".scores.csv", sep=""))
   y <- plotIndiv(plsda.metabolite, ellipse=TRUE, legend=TRUE)
   write.csv(y$df, outputfile)
   perf.plsda <- perf(plsda.metabolite, validation="Mfold", folds=5, progressBar=FALSE, auc=TRUE, nrepeat=10)
   #print(str(perf.plsda))
   print(perf.plsda$error.rate$overall)
   plot(perf.plsda, col=color.mixo(1:3), sd=TRUE, ylim=c(0,1), legend.position="horizontal")
  #plotIndiv(plsda.metabolite, ind.names = F,
  #        add.legend =TRUE, plot.ellipse = TRUE,
  #        ellipse.level = 0.5, blocks = "lipid", main = 'PLSDA',
  #        plot.star = TRUE, plot.centroid = TRUE,style='3d')
}
