checkmembership <- function(kclust,kclusttrue,truelabels,estmemb) {
  relabelmatrix<-matrix(nrow=kclust,ncol=kclusttrue) #i row_num label, j col_num estmemb
  estlabels<-list()
  for (i in 1:kclust) {
    estlabels[[i]]<-which(estmemb==i)
  }
  for (i in 1:kclust) {
    for (j in 1:kclusttrue) {
      relabelmatrix[i,j]<-length(which(estlabels[[i]]%in%truelabels[[j]]))
    }
  }
  rowcol<-solve_LSAP(relabelmatrix,TRUE)
  res<-list()
  res$relabel<-as.vector(rowcol)
  res$ncorr<-0
  for (j in 1:min(kclust,kclusttrue)) {
    res$ncorr<-res$ncorr+relabelmatrix[j,res$relabel[j]]
  }
  res$relabelmatrix<-relabelmatrix
  return(res)
}
generatetriple<-function(n) {
  resmat<-matrix(nrow=n,ncol=3)
  for (i in 1:n) {
    ordr<-sample.int(3,3)
    res<-vector(length=3)
    res[1]<-runif(1,min=0,max=1)
    res[2]<-runif(1,min=0,max=1-res[1])
    res[3]<-1-res[1]-res[2]
    resmat[i,ordr]<-res
  }
  return(resmat)
}
generatefour<-function(n) {
  resmat<-matrix(nrow=n,ncol=4)
  for (i in 1:n) {
    ordr<-sample.int(4,4)
    res<-vector(length=4)
    res[1]<-runif(1,min=0,max=1)
    res[2]<-runif(1,min=0,max=1-res[1])
    res[3]<-runif(1,min=0,max=1-res[1]-res[2])
    res[4]<-1-res[1]-res[2]-res[3]
    resmat[i,ordr]<-res
  }
  return(resmat)
}

propersample <- function(x){if(length(x)==1) x else sample(x,1)}

calcloglike <- function(samplescores,tau) {
 # samplescores<-t(samplescores)
  maxscorey<-apply(samplescores,1,max) # find the max of each column
  loglike<-sum(log(colSums(t(exp(samplescores-maxscorey))*tau))+maxscorey) # remove max for numerical stability and exponentiate
  return(loglike)
}

reassignsamples <- function(samplescores,numsamps){
  newclustermembership <-rep(0,numsamps) # to store the new cluster
  for(s in 1:numsamps){ # run through the samples
    clusterscores<-samplescores[s,]
    maxscorey<-max(clusterscores) # take the maximum
    maxscoreelem<-which(clusterscores==maxscorey)
    newclustermembership[s]<-propersample(maxscoreelem)
  }  
  return(newclustermembership)
}
relativeprobs <- function(samplescores,numsamps){
  relativeprobabs <-rep(0,numsamps) # to store the relative probabilities
  for(s in 1:numsamps){ # run through the samples
    clusterscores<-samplescores[s,]
    maxscorey<-max(clusterscores) # take the maximum
    shifty<-exp(clusterscores-maxscorey)
    rescaley<-shifty/sum(shifty)
    relativeprobabs[s]<-max(rescaley) # relative probabilities
  }  
  return(relativeprobabs)
}

allrelativeprobs <- function(samplescores,numsamps){
  relativeprobabs <-samplescores # to store the relative probabilities
  for(s in 1:numsamps){ # run through the samples
    clusterscores<-samplescores[s,]
    maxscorey<-max(clusterscores) # take the maximum
    shifty<-exp(clusterscores-maxscorey)
    rescaley<-shifty/sum(shifty)
    relativeprobabs[s,]<-rescaley # relative probabilities
  }  
  return(relativeprobabs)
}

relativeprobswithtau <- function(sampleprobs,tau){
  temp<-tau*t(sampleprobs)
  relativeprobabswithtau<-1/colSums(temp)*t(temp) # ugly code
  return(relativeprobabswithtau)
}

avescore <- function(samplescores,numsamps){
  averagescores <-rep(0,numsamps) # to store the relative probabilities
  for(s in 1:numsamps){ # run through the samples
    clusterscores<-samplescores[s,]
    maxscorey<-max(clusterscores) # take the maximum
    shifty<-exp(clusterscores-maxscorey) # exponentiate
    rescaley<-log(mean(shifty))+maxscorey # mean and turn back to log
    averagescores[s]<-rescaley # averagescore
  }  
  return(averagescores)
}
reassignsamplesprop <- function(samplescores,numsamps,gamma){
  newclustermembership <-rep(0,numsamps) # to store the new cluster
  for(s in 1:numsamps){ # run through the samples
    clusterscores<-samplescores[s,]*gamma
    maxscorey<-max(clusterscores) # take the maximum
    shifty<-exp(clusterscores-maxscorey)
    rescaley<-shifty/sum(shifty)
    scorelength<-length(rescaley)
    newclustermembership[s]<-sample.int(scorelength,1,prob=rescaley) # sample according to scores
  }  
  return(newclustermembership)
}
