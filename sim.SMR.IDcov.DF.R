e2dist<- function (x, y){
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}

sim.SMR.IDcov.DF<-
  function(N=50,n.marked=10,lam0=NA,theta.d=NA,sigma=0.50,K=10,X=X,buff=3,
           theta.marked=c(1,0,0),theta.unmarked=1,K1D=NA,tlocs=0,marktype="natural",obstype="poisson",
           n.cat=NA,IDcovs=NA,gamma=NA,theta.cat=NA){
    if(!marktype%in%c("natural","premarked")){
      stop("marktype must be 'natural' or 'premarked'")
    }
    if(sum(theta.marked)!=1)stop("theta.marked must sum to 1.")
    if(theta.unmarked<0|theta.unmarked>1)stop("theta.unmarked must be between 0 and 1.")
    if(length(lam0)!=length(sigma))stop("Need the same number of values for each detection function parameter")
    if(length(IDcovs[[1]])!=length(lam0))stop("IDcovs[[1]] should have the same number of values as detection function parameter values")
    library(abind)
    # simulate a population of activity centers
    X=as.matrix(X)
    xlim=c(min(X[,1]),max(X[,1]))+c(-buff,buff)
    ylim=c(min(X[,2]),max(X[,2]))+c(-buff,buff)
    s<- cbind(runif(N, xlim[1],xlim[2]), runif(N,ylim[1],ylim[2]))
    D<- e2dist(s,X)
    J=nrow(X)
    
    #trap operation
    if(!any(is.na(K1D))){
      if(any(K1D>K)){
        stop("Some entries in K1D are greater than K.")
      }
      if(is.null(dim(K1D))){
        if(length(K1D)!=J){
          stop("K1D vector must be of length J.")
        }
      }
    }else{
      K1D=rep(K,J)
    }
    
    #simulate IDcovs
    G.true=matrix(NA,nrow=N,ncol=n.cat) #all IDcovs in population.
    for(i in 1:N){
      for(j in 1:n.cat){
        G.true[i,j]=sample(IDcovs[[j]],1,prob=gamma[[j]])
      }
    }
    
    # Capture and mark individuals
    y <-array(0,dim=c(N,J,K))
    if(obstype=="poisson"){
      if(any(is.na(lam0)))stop("must provide lam0 for poisson obstype")
      lamd <-lam0[1] * exp(-D * D/(2 * sigma[1] * sigma[1]))
      for(i in 2:length(sigma)){
        lamd[which(G.true[,1]==i), ]<-lam0[i]*exp(-D[which(G.true[,1]==i), ]^2/(2*sigma[i]^2))
      }
      for(i in 1:N){
        for(j in 1:J){
          y[i,j,1:K1D[j]]=rpois(K1D[j],lamd[i,j])
        }
      } 
    }else if(obstype=="negbin"){
      if(any(is.na(lam0)))stop("must provide lam0 for negbin obstype")
      if(is.na(theta.d))stop("Must provide theta.d for negbin obstype")
      lamd <-lam0[1] * exp(-D * D/(2 * sigma[1] * sigma[1]))
      for(i in 2:length(sigma)){
        lamd[which(G.true[,1]==i), ]<-lam0[i]*exp(-D[which(G.true[,1]==i), ]^2/(2*sigma[i]^2))
      }
      for(i in 1:N){
        for(j in 1:J){
          for(k in 1:K){
            y[i,j,1:K1D[j]]=rnbinom(K1D[j],mu=lamd[i,j],size=theta.d)
          }
        }
      } 
    }else{
      stop("obstype not recognized")
    }

    if(marktype=="natural"){
      #reorder data so enough marked guys are at the top
      #must be random sample of marked guys, e.g. not ordered by # of captures
      idx=which(rowSums(y)>0)
      if(length(idx)<n.marked){
        stop("Fewer than n.marked individuals captured. Cannot naturally mark uncaptured individuals.")
      }
      move=sample(idx,n.marked)
      idx=setdiff(1:N,move)
      y=y[c(move,idx),,]
      s=s[c(move,idx),]
      G.true[c(move,idx),]
    }
    
    IDmarked=1:n.marked
    umguys=setdiff(1:N,IDmarked)
    
    #split sightings into marked and unmarked histories, considering occasion of marking
    y.marked=y[IDmarked,,]
    G.marked=G.true[IDmarked,]
    n.samples=sum(y[umguys,,])
    G.unmarked=matrix(NA,nrow=n.samples,ncol=n.cat)
    y.unmarked=array(0,dim=c(n.samples,J,K))
    IDum=rep(NA,n.samples)
    idx=1
    for(i in 1:length(umguys)){
      for(j in 1:J){ #then traps
        for(k in 1:K){ #then occasions
          if(y[umguys[i],j,k]>0){ #is there at least one sample here?
            for(l in 1:y[umguys[i],j,k]){ #then samples
              G.unmarked[idx,]=G.true[umguys[i],]
              y.unmarked[idx,j,k]=1
              IDum[idx]=umguys[i]
              idx=idx+1
            }
          }
        }
      }
    }
    if(!is.matrix(G.marked)){
      G.marked=matrix(G.marked,ncol=1)
    }
    if(!is.matrix(G.unmarked)){
      G.unmarked=matrix(G.unmarked,ncol=1)
    }
    
    #ID/marked status observation model for marked individuals
    if(theta.marked[1]!=0&sum(y.marked)>0){#bug fix from Glenn Stauffer
      idx1=which(y.marked>0)#used to extract counts
      count=y.marked[idx1]
      idx2=which(y.marked>0,arr.ind=TRUE)#used to move counts
      idx3=rep(1:nrow(idx2),count)#repeat array indices for counts>1
      idx2=idx2[idx3,]
      if(!is.matrix(idx2)){
        idx2=matrix(idx2,ncol=3)
      }
      n.sample.mark=sum(count)
      outcome=rmultinom(n.sample.mark,1,theta.marked)
      mnoID.idx=which(outcome[2,]==1)
      munk.idx=which(outcome[3,]==1)
      
      #remove marked but no ID and unk marked status from y.marked
      #create y.marked.noID and y.unk
      nmnoID=length(mnoID.idx)
      nmunk=length(munk.idx)
      if(nmnoID>0){
        y.marked.noID=array(0,dim=c(nmnoID,J,K))
        G.marked.noID=matrix(0,nrow=nmnoID,ncol=n.cat)
        for(i in 1:length(mnoID.idx)){
          #delete this sighting
          y.marked[idx2[mnoID.idx[i],1],idx2[mnoID.idx[i],2],idx2[mnoID.idx[i],3]]=
            y.marked[idx2[mnoID.idx[i],1],idx2[mnoID.idx[i],2],idx2[mnoID.idx[i],3]]-1
          #add this sighting
          y.marked.noID[i,idx2[mnoID.idx[i],2],idx2[mnoID.idx[i],3]]=1
          #extract G
          G.marked.noID[i,]=G.marked[idx2[mnoID.idx[i],1],]
        }
        IDmnoID=idx2[mnoID.idx,1]
      }else{
        y.marked.noID=G.marked.noID=IDnoID=NA
        IDmnoID=NA
      }
      if(nmunk>0){
        y.unk=array(0,dim=c(nmunk,J,K))
        G.unk=matrix(0,nrow=nmunk,ncol=n.cat)
        for(i in 1:length(munk.idx)){
          #delete this sighting
          y.marked[idx2[munk.idx[i],1],idx2[munk.idx[i],2],idx2[munk.idx[i],3]]=
            y.marked[idx2[munk.idx[i],1],idx2[munk.idx[i],2],idx2[munk.idx[i],3]]-1
          #add this sighting
          y.unk[i,idx2[munk.idx[i],2],idx2[munk.idx[i],3]]=1
          #extract G
          G.unk[i,]=G.marked[idx2[munk.idx[i],1],]
        }
        IDunk=idx2[munk.idx,1]
        IDunkType=rep("marked",length(munk.idx))
      }else{
        y.unk=G.unk=IDunk=IDunkType=NA
      }
    }else{
      y.marked.noID=G.marked.noID=IDmnoID=NA
      y.unk=G.unk=IDunk=IDunkType=NA
    }
    
    #marked status observation model for unmarked individuals
    if(theta.unmarked!=1){
      outcome=rbinom(n.samples,1,theta.unmarked)
      unk.idx=which(outcome==0)
      nunk=length(unk.idx)
      if(nunk>0){
        #extract um history to unk
        y.unk2=y.unmarked[unk.idx,,]
        G.unk2=G.unmarked[unk.idx,]
        IDunk2=IDum[unk.idx]
        IDunkType2=rep("unmarked",length(IDunk2))
        #remove unk from um history
        y.unmarked=y.unmarked[-unk.idx,,]
        G.unmarked=G.unmarked[-unk.idx,]
        IDum=IDum[-unk.idx]
      }else{
        IDunk2=IDunkType2=NA
      }
    }else{
      y.unk2=NA
      G.unk2=NA
      IDunk2=IDunkType2=NA
    }
    #If only 1 G cat, convert back to matrix
    if(n.cat==1){
      G.unmarked=matrix(G.unmarked,ncol=1)
      if(!all(is.na(G.unk))){
        G.unk=matrix(G.unk,ncol=1)
      }
      if(!all(is.na(G.unk2))){
        G.unk2=matrix(G.unk2,ncol=1)
      }
      if(!all(is.na(G.marked.noID))){
        G.marked.noID=matrix(G.marked.noID,ncol=1)
      }
    }
    
    #combine y.unk and G.sight.unk if there are members from both marked and unmarked
    if(!any(is.na(IDunk))&!any(is.na(IDunk2))){
      IDunk=c(IDunk,IDunk2)
      IDunkType=c(IDunkType,IDunkType2)
      y.unk=abind(y.unk,y.unk2,along=1)
      G.unk=rbind(G.unk,G.unk2)
    }else if(!any(is.na(IDunk2))){#or if just unk from unmarked, rename them
      IDunk=IDunk2
      IDunkType=IDunkType2
      y.unk=y.unk2
      G.unk=G.unk2
    }
    
    #Observation failure for category levels in unmarked, unknown, and marked no ID if present
    for(l in 1:n.cat){
      drop=which(rbinom(nrow(G.unmarked),1,theta.cat[l])==0)
      if(length(drop)>0){
        G.unmarked[drop,l]=0 #coding missing data as 0
      }
    }
    if(!is.na(G.unk[1])){
      for(l in 1:n.cat){
        drop=which(rbinom(nrow(G.unk),1,theta.cat[l])==0)
        if(length(drop)>0){
          G.unk[drop,l]=0 #coding missing data as 0
        }
      }
    }
    if(!is.na(G.marked.noID[1])){
      for(l in 1:n.cat){
        drop=which(rbinom(nrow(G.marked.noID),1,theta.cat[l])==0)
        if(length(drop)>0){
          G.marked.noID[drop,l]=0 #coding missing data as 0
        }
      }
    }
    
    #check data
    y.check=y*0
    for(i in 1:length(IDmarked)){
      y.check[IDmarked[i],,]=y.marked[i,,]
    }
    if(length(IDum)>0){
      for(i in 1:length(IDum)){
        y.check[IDum[i],,]=y.check[IDum[i],,]+y.unmarked[i,,]
      }
    }
    if(all(!is.na(IDunk))){
      for(i in 1:length(IDunk)){
        y.check[IDunk[i],,]=y.check[IDunk[i],,]+y.unk[i,,]
      }
    }
    if(theta.marked[2]>0){
      if(all(is.finite(IDmnoID))){
        for(i in 1:length(IDmnoID)){
          y.check[IDmnoID[i],,]=y.check[IDmnoID[i],,]+y.marked.noID[i,,]
        }
      }
    }
    if(!all(y==y.check)){
      stop("Error rebuilding data. Report bug.")
    }
    dimnames(y.unk)=NULL
    if(all(!is.na(G.unk))){
      for(i in 1:nrow(G.unk)){
        for(m in 1:n.cat){
          if(G.unk[i,m]!=0){
            if(G.true[IDunk[i],m]!=G.unk[i,m])stop("G.unk not in correct order. Report bug.")
          }
        }
      }
    }
    if(all(!is.na(G.unmarked))){
      for(i in 1:nrow(G.unmarked)){
        for(m in 1:n.cat){
          if(G.unmarked[i,m]!=0){
            if(G.true[IDum[i],m]!=G.unmarked[i,m])stop("G.unmarked not in correct order. Report bug.")
          }
        }
      }
    }
    if(all(!is.na(G.marked.noID))){
      for(i in 1:nrow(G.marked.noID)){
        for(m in 1:n.cat){
          if(G.marked.noID[i,m]!=0){
            if(G.true[IDmnoID[i],m]!=G.marked.noID[i,m])stop("G.marked.noID not in correct order. Report bug.")
          }
        }
      }
    }
    
    #Telemetry observations
    if(tlocs>0){
      if(marktype=="natural")warning("Simulating telemetry for naturally marked individuals, but this probably does not make sense.")
      locs=array(NA,dim=c(n.marked,tlocs,2))
      for(i in 1:n.marked){
        for(j in 1:tlocs){
          locs[i,j,]=c(rnorm(1,s[IDmarked[i],1],sigma[G.true[IDmarked[i],1]]),rnorm(1,s[IDmarked[i],2],sigma[G.true[IDmarked[i],1]]))
        }
      }
    }else{
      locs=NA
    }
    
    #convert unknown ID observations to "this.j" and "this.k"
    #Are there unknown marked status guys?
    useUnk=FALSE
    if(!all(is.na(y.unk))){
        useUnk=TRUE
    }else{
      y.unk=array(0,dim=c(0,J,K))
      ID.unk=c()
    }
    #Are there marked no ID guys?
    useMarkednoID=FALSE
    if(!all(is.na(y.marked.noID))){
      useMarkednoID=TRUE
    }else{
      y.marked.noID=array(0,dim=c(0,J,K))
      ID.mnoID=c()
    }
    
    #disassemble y.marked
    G.marked.ID=matrix(NA,nrow=sum(y.marked),ncol=n.cat)
    y.marked.ID=array(0,dim=c(sum(y.marked),J,K))
    ID.marked=rep(NA,sum(y.marked))
    idx=1
    for(i in 1:n.marked){
      for(j in 1:J){
        for(k in 1:K){
          if(y.marked[i,j,k]>0){
            for(l in 1:y.marked[i,j,k]){
              y.marked.ID[idx,j,k]=1
              ID.marked[idx]=i
              G.marked.ID[idx,]=G.marked[i,]
              idx=idx+1
            }
          }
        }
      }
    }
    
    n.samp1=nrow(y.marked.ID) #1
    if(useMarkednoID){
      n.samp2=nrow(y.marked.noID) #2
    }else{
      n.samp2=0
    }
    n.samp3=nrow(y.unmarked) #3
    if(useUnk){
      n.samp4=nrow(y.unk) #4
    }else{
      n.samp4=0
    }
    n.samples=n.samp1+n.samp2+n.samp3+n.samp4
    
    y.obs=abind(y.marked.ID,y.marked.noID,y.unmarked,y.unk,along=1)
    y.obs2D=apply(y.obs,c(1,2),sum)
    this.j=apply(y.obs2D,1,function(x){which(x>0)})
    y.obs2Dk=apply(y.obs,c(1,3),sum)
    this.k=apply(y.obs2Dk,1,function(x){which(x>0)})
    
    samp.type=c(rep("markedID",n.samp1),
                rep("markednoID",n.samp2),
                rep("unmarked",n.samp3),
                rep("unk",n.samp4))
    ID=c(ID.marked,IDmnoID,IDum,IDunk)
    n.M=sum(rowSums(y[1:n.marked,,])>0)
    n.UM=sum(rowSums(y[(n.marked+1):N,,])>0)
    
    #actually, let's combine all types into G.obs here
    if(!any(is.na(G.unk))&!any(is.na(G.marked.noID))){
      G.obs=rbind(G.marked.ID,G.marked.noID,G.unmarked,G.unk)
    }else if(!any(is.na(G.unk))&any(is.na(G.marked.noID))){
      G.obs=rbind(G.marked.ID,G.unmarked,G.unk)
    }else if(any(is.na(G.unk))&!any(is.na(G.marked.noID))){
      G.obs=rbind(G.marked.ID,G.marked.noID,G.unmarked)
    }else{
      G.obs=rbind(G.marked.ID,G.unmarked)
    }
    

    out<-list(this.j=this.j,this.k=this.k,samp.type=samp.type,ID.marked=ID.marked, #observed data
              G.marked=G.marked,G.obs=G.obs, #observed IDcov data
              n.marked=n.marked,locs=locs,n.M=n.M,n.UM=n.UM,IDlist=list(n.cat=n.cat,IDcovs=IDcovs),
              y=y,s=s, ID=ID,#true data
              X=X,K=K,K1D=K1D,buff=buff,xlim=xlim,ylim=ylim)
    return(out)
  }