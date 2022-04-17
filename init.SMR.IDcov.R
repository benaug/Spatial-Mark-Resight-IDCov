e2dist<- function (x, y){
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}

init.SMR.IDcov=function(data,inits=NA,M1=NA,M2=NA,marktype="premarked",obstype="poisson"){
  library(abind)
  #extract observed data
  this.j=data$this.j
  # this.k=data$this.k #not used in this 2D sampler
  samp.type=data$samp.type
  ID.marked=data$ID.marked
  n.marked=data$n.marked
  X<-as.matrix(data$X)
  J<-nrow(X)
  K<- data$K
  K1D=data$K1D
  buff<- data$buff
  M.both=M1+M2
  locs=data$locs
  G.marked=data$G.marked
  G.obs=data$G.obs
  n.cat=data$IDlist$n.cat
  IDcovs=data$IDlist$IDcovs
  n.levels=unlist(lapply(IDcovs,length))
  if(!is.matrix(G.marked)){
    G.marked=matrix(G.marked)
  }
  if(!is.matrix(G.obs)){
    G.obs=matrix(G.obs)
  }
  if(!is.list(IDcovs)){
    stop("IDcovs must be a list")
  }
  if(ncol(G.marked)!=n.cat){
    stop("G.marked needs n.cat number of columns")
  }
  if(ncol(G.obs)!=n.cat){
    stop("G.obs needs n.cat number of columns")
  }
  
  xlim<- c(min(X[,1]),max(X[,1]))+c(-buff, buff)
  ylim<- c(min(X[,2]),max(X[,2]))+c(-buff, buff)
  
  ##pull out initial values
  lam0=inits$lam0
  sigma<- inits$sigma
  
  n.samp1=sum(samp.type=="markedID")
  n.samp2=sum(samp.type=="markednoID")
  n.samp3=sum(samp.type=="unmarked")
  n.samp4=sum(samp.type=="unk")
  n.samples=length(this.j)
  
  useMarkednoID=FALSE
  if(n.samp2>0){
    useMarkednoID=TRUE
  }
  useUnk=FALSE
  if(n.samp4>0){
    useUnk=TRUE
  }
  
  #build y.marked
  y.marked=matrix(0,M1,J)
  for(l in 1:length(ID.marked)){
    y.marked[ID.marked[l],this.j[l]]=y.marked[ID.marked[l],this.j[l]]+1
  }
  
  G.type=rep(c(1,1,2,0),times=c(n.samp1,n.samp2,n.samp3,n.samp4))
  G.obs=cbind(G.type,G.obs)
  
  #initialize unknown IDs
  G.true=matrix(0,nrow=M.both,ncol=n.cat)
  G.true[1:n.marked,]=G.marked
  G.true=cbind(c(rep(1,M1),rep(2,M2)),G.true)
  ID=c(ID.marked,rep(NA,n.samples-length(ID.marked)))
  nextID=max(ID,na.rm=TRUE)+1
  
  y.true2D=apply(y.marked,c(1,2),sum)
  if(marktype=="natural"){
    y.true2D=rbind(y.true2D,matrix(0,nrow=M2,ncol=J))
  }else{
    y.true2D=rbind(y.true2D,matrix(0,nrow=M1+M2-n.marked,ncol=J))
  }
  if(M1<n.marked)stop("M1 must be larger than the number of marked individuals.")
  
  #Make sure G.obs for marked ID samples matches G.marked 
  for(l in 1:(n.samp1)){
    obsidx1=which(G.obs[l,]!=0)
    obsidx2=which(G.true[ID.marked[l],]!=0)
    obsidx.use=intersect(obsidx1,obsidx2)
    if(!all(G.true[ID.marked[l],obsidx.use]==G.obs[l,obsidx.use]))stop(paste("G.obs for sample",l,"does not match the corresponding G.true"))
    if(any(G.true[ID.marked[l],obsidx1]==0))stop(paste("G.obs for sample",l,"implies a corresponding element of G.true is not actually missing. (coded as 0 when actually known"))

  }
  
  
  
  #marked noID
  if(useMarkednoID){
    if(marktype=="natural"){
      for(l in (n.samp1+1):(n.samp1+n.samp2)){
        #can you match an ID'd guy in same trap?
        obsidx=which(G.obs[l,]!=0)
        matches=rep(FALSE,M1)
        for(i2 in 1:M1){#can only match marked
          obsidx2=which(G.true[i2,]!=0)
          sameobsidx=intersect(obsidx,obsidx2)
          matches[i2]=all(G.true[i2,sameobsidx]==G.obs[l,sameobsidx])
        }
        matches=which(matches)
        if(length(matches)==0){#must be new ID
          ID[l]=nextID
          G.true[ID[l],obsidx]=G.obs[l,obsidx]
          y.true2D[ID[l],this.j[l]]=y.true2D[ID[l],this.j[l]]+1
          nextID=nextID+1
        }else if(length(matches)==1){
          if(y.true2D[matches,this.j[l]]>0){#caught at same trap?
            ID[l]=matches
            y.true2D[ID[l],this.j[l]]=y.true2D[ID[l],this.j[l]]+1
            #new sample for this ID might fill in some missing G.true indices
            notobsidx=which(G.true[ID[l],]==0)
            G.true[ID[l],notobsidx]=G.obs[l,notobsidx]
          }else{#must be new ID
            ID[l]=nextID
            G.true[ID[l],obsidx]=G.obs[l,obsidx]
            y.true2D[ID[l],this.j[l]]=y.true2D[ID[l],this.j[l]]+1
            nextID=nextID+1
          }
        }else{
          sametrap=y.true2D[matches,this.j[l]]>0
          if(any(sametrap)){
            ID[l]=matches[which(sametrap)[1]]
            y.true2D[ID[l],this.j[l]]=y.true2D[ID[l],this.j[l]]+1
            #new sample for this ID might fill in some missing G.true indices
            notobsidx=which(G.true[ID[l],]==0)
            G.true[ID[l],notobsidx]=G.obs[l,notobsidx]
          }else{#must be new ID
            ID[l]=nextID
            G.true[ID[l],obsidx]=G.obs[l,obsidx]
            y.true2D[ID[l],this.j[l]]=y.true2D[ID[l],this.j[l]]+1
            nextID=nextID+1
          }
        }
        if(nextID>M1)stop("Need to raise M1 to initialize data (marktype=natural)")
      }
    }else{ #premarked requires different approach because number of marks known
      for(l in (n.samp1+1):(n.samp1+n.samp2)){
        obsidx=which(G.obs[l,]!=0)
        matches=rep(FALSE,M1)
        for(i2 in 1:M1){#can only match marked
          obsidx2=which(G.true[i2,]!=0)
          sameobsidx=intersect(obsidx,obsidx2)
          matches[i2]=all(G.true[i2,sameobsidx]==G.obs[l,sameobsidx])
        }
        matches=which(matches)
        #chose match with closest other capture, if you match a captured ID
        assignedIDsthatmatch=which(ID%in%matches)
        if(length(assignedIDsthatmatch)>0){
          thesetraps=X[this.j[assignedIDsthatmatch],]
          if(length(assignedIDsthatmatch)>1){
            dists=sqrt((X[this.j[l],1]-thesetraps[,1])^2+(X[this.j[l],2]-thesetraps[,2])^2)
          }else{
            dists=sqrt((X[this.j[l],1]-thesetraps[1])^2+(X[this.j[l],2]-thesetraps[2])^2)
          }
          pick=which(dists==min(dists))[1]#pick first one if multiple to choose from
          ID[l]=ID[assignedIDsthatmatch[pick]]
        }else{#if none of your matches were captured, pick first match
          ID[l]=matches[1]
        }
        #new sample for this ID might fill in some missing G.true indices
        notobsidx=which(G.true[ID[l],]==0)
        G.true[ID[l],notobsidx]=G.obs[l,notobsidx]
      }
    }
  }
  #unmarked next...
  nextID=M1+1
  for(l in (n.samp1+n.samp2+1):(n.samp1+n.samp2+n.samp3)){
    #can you match an unmarked guy in same trap already assigned an ID?
    obsidx=which(G.obs[l,]!=0)
    matches=rep(FALSE,M.both)
    for(i2 in (M1+1):M.both){#can only match unmarked
      obsidx2=which(G.true[i2,]!=0)
      sameobsidx=intersect(obsidx,obsidx2)
      matches[i2]=all(G.true[i2,sameobsidx]==G.obs[l,sameobsidx])
    }
    matches=which(matches)
    # matches=which(apply(G.true[(M1+1):(M.both),obsidx],1,function(x){all(x==G.obs[l,obsidx])}))+M1 #can only match unmarked
    if(length(matches)==0){#must be new ID
      ID[l]=nextID
      G.true[ID[l],obsidx]=G.obs[l,obsidx]
      y.true2D[ID[l],this.j[l]]=y.true2D[ID[l],this.j[l]]+1
      nextID=nextID+1
    }else if(length(matches)==1){
      if(y.true2D[matches,this.j[l]]>0){#caught at same trap?
        ID[l]=matches
        y.true2D[ID[l],this.j[l]]=y.true2D[ID[l],this.j[l]]+1
        #new sample for this ID might fill in some missing G.true indices
        notobsidx=which(G.true[ID[l],]==0)
        G.true[ID[l],notobsidx]=G.obs[l,notobsidx]
      }else{#must be new ID
        ID[l]=nextID
        G.true[ID[l],obsidx]=G.obs[l,obsidx]
        y.true2D[ID[l],this.j[l]]=y.true2D[ID[l],this.j[l]]+1
        nextID=nextID+1
      }
    }else{
      sametrap=y.true2D[matches,this.j[l]]>0
      if(any(sametrap)){
        ID[l]=matches[which(sametrap)[1]]
        y.true2D[ID[l],this.j[l]]=y.true2D[ID[l],this.j[l]]+1
        #new sample for this ID might fill in some missing G.true indices
        notobsidx=which(G.true[ID[l],]==0)
        G.true[ID[l],notobsidx]=G.obs[l,notobsidx]
      }else{#must be new ID
        ID[l]=nextID
        G.true[ID[l],obsidx]=G.obs[l,obsidx]
        y.true2D[ID[l],this.j[l]]=y.true2D[ID[l],this.j[l]]+1
        nextID=nextID+1
      }
    }
    if(nextID>M.both)stop("Need to raise M2 to initialize data.")
  }
  #then unknown...
  #keeping nextID where it is so we assign non matches to unmarked class
  if(useUnk){
    for(l in (n.samp1+n.samp2+n.samp3+1):(n.samp1+n.samp2+n.samp3+n.samp4)){
      #can you match an unmarked guy in same trap already assigned an ID?
      obsidx=which(G.obs[l,]!=0)
      matches=rep(FALSE,M.both)
      for(i2 in 1:M.both){#can match marked or unmarked
        obsidx2=which(G.true[i2,]!=0)
        sameobsidx=intersect(obsidx,obsidx2)
        matches[i2]=all(G.true[i2,sameobsidx]==G.obs[l,sameobsidx])
      }
      matches=which(matches)
      # matches=which(apply(G.true[,obsidx],1,function(x){all(x==G.obs[l,obsidx])})) #can match marked or unmarked
      if(length(matches)==0){#must be new ID
        ID[l]=nextID
        G.true[ID[l],obsidx]=G.obs[l,obsidx]
        y.true2D[ID[l],this.j[l]]=y.true2D[ID[l],this.j[l]]+1
        nextID=nextID+1
      }else if(length(matches)==1){
        if(y.true2D[matches,this.j[l]]>0){#caught at same trap?
          ID[l]=matches
          y.true2D[ID[l],this.j[l]]=y.true2D[ID[l],this.j[l]]+1
          #new sample for this ID might fill in some missing G.true indices
          notobsidx=which(G.true[ID[l],]==0)
          G.true[ID[l],notobsidx]=G.obs[l,notobsidx]
        }else{#must be new ID
          ID[l]=nextID
          G.true[ID[l],obsidx]=G.obs[l,obsidx]
          y.true2D[ID[l],this.j[l]]=y.true2D[ID[l],this.j[l]]+1
          nextID=nextID+1
        }
      }else{
        sametrap=y.true2D[matches,this.j[l]]>0
        if(any(sametrap)){
          ID[l]=matches[which(sametrap)[1]]
          y.true2D[ID[l],this.j[l]]=y.true2D[ID[l],this.j[l]]+1
          #new sample for this ID might fill in some missing G.true indices
          notobsidx=which(G.true[ID[l],]==0)
          G.true[ID[l],notobsidx]=G.obs[l,notobsidx]
        }else{#must be new ID
          ID[l]=nextID
          G.true[ID[l],obsidx]=G.obs[l,obsidx]
          y.true2D[ID[l],this.j[l]]=y.true2D[ID[l],this.j[l]]+1
          nextID=nextID+1
        }
      }
      if(nextID>M.both)stop("Need to raise M2 to initialize data.")
    }
  }
  
  #Check for correctness
  for(l in 1:n.samples){
    obsidx=which(G.obs[l,]!=0)
    if(!all(G.true[ID[l],obsidx]==G.obs[l,obsidx]))stop("Error in initialization")
  }
  
  #initialize remainder of G.true
  for(i in 1:M.both){
    for(c in 2:(n.cat+1)){
      if(G.true[i,c]==0){
        G.true[i,c]=sample(IDcovs[[c-1]],1,prob=inits$gamma[[c-1]])
      }
    }
  }
  
  #initialize match
  match=matrix(FALSE,nrow=n.samples,ncol=M1+M2)
  for(l in 1:n.samples){
    idx=which(G.obs[l,]!=0)
    if(length(idx)>1){#multiple observed
      match[l,]=apply(G.true[,idx],1,function(x){all(x==G.obs[l,idx])})
    }else if(length(idx)==1){#single observed
      match[l,]=G.true[,idx]==G.obs[l,idx]
    }else{#fully latent G.obs
      match[l,]=rep(TRUE,M.both)
    }
  }
  
  #initialize z
  z=1*(rowSums(y.true2D)>0)
  if(marktype=="premarked"){
    z[1:M1]=1
  }
  
  #intialize s
  s<- cbind(runif(M.both,xlim[1],xlim[2]), runif(M.both,ylim[1],ylim[2])) #assign random locations
  idx=which(rowSums(y.true2D)>0) #switch for those actually caught
  for(i in idx){
    trps<- matrix(X[y.true2D[i,]>0,1:2],ncol=2,byrow=FALSE)
    if(nrow(trps)>1){
      s[i,]<- c(mean(trps[,1]),mean(trps[,2]))
    }else{
      s[i,]<- trps
    }
  }
  
  #identify known ID (1), known status, unknown ID (2), unknown ID and status (3)
  G.type=rep(0,n.samples)
  G.type[samp.type=="markedID"]=1
  G.type[samp.type%in%c("markednoID","unmarked")]=2
  G.type[samp.type=="unk"]=3
  n.fixed=sum(samp.type=="markedID")

  #Stuff category level probabilities into a ragged matrix. (may have different numbers of levels)
  gammaMat=matrix(0,nrow=n.cat,ncol=max(n.levels))
  for(l in 1:n.cat){
    gammaMat[l,1:n.levels[l]]=gamma[[l]]
  }
  
  #calculate G.latent, which indicies of G currently latent
  G.latent=matrix(1,nrow=M.both,ncol=n.cat+1)
  G.latent[1:M.both,1]=0 #never change column for marked status
  for(l in 1:n.samples){
    for(m in 2:(n.cat+1)){
      if(G.obs[l,m]!=0){
        G.latent[ID[l],m]=0
      }
    }
  }
  #calculate G.latent for only marked individuals
  G.latent.marked=matrix(1,nrow=n.marked,ncol=n.cat+1)
  G.latent.marked[1:n.marked,1]=0 #never change column for marked status
  G.latent.marked[1:n.marked,2:(n.cat+1)]=(1*(G.marked==0))
  
  if(!is.null(dim(data$locs))){
    max.locs=dim(locs)[2]
    tel.inds=which(rowSums(is.na(locs[,,1]))<max.locs)
    n.locs.ind=rowSums(!is.na(locs[,,1]))
    n.locs.ind=n.locs.ind[tel.inds]
    print("using telemetry to initialize telmetered s. Remove from data if not using in the model.")
    #update s starts for telemetry guys
    for(i in tel.inds){
      if(n.locs.ind[i]>1){
        s[i,]=colMeans(locs[i,1:n.locs.ind[i],])
      }else{
        s[i,]=locs[i,1,]
      }
      #make sure s is in state space
      if(s[i,1]<xlim[1]){
        s[i,1]=xlim[1]
      }
      if(s[i,1]>xlim[2]){
        s[i,1]=xlim[2]
      }
      if(s[i,2]<ylim[1]){
        s[i,2]=ylim[1]
      }
      if(s[i,2]>ylim[2]){
        s[i,2]=ylim[2]
      }
    }
  }else{
    tel.inds=NA
    n.locs.ind=NA
  }
  
  D=e2dist(s, X)
  lamd<- lam0*exp(-D*D/(2*sigma*sigma))
  #make 3D y.true
  y.true3D=array(0,dim=c(M.both,J,3))
  for(l in 1:n.samples){
    y.true3D[ID[l],this.j[l],G.type[l]]=y.true3D[ID[l],this.j[l],G.type[l]]+1
  }
  y.true2D=apply(y.true3D,c(1,2),sum)
  
  if(obstype=="poisson"){
    ll.y=dpois(y.true2D,K1D*lamd*z,log=TRUE)
  }else if(obstype=="negbin"){
    theta.d=inits$theta.d
    ll.y=y.true2D*0
    for(i in 1:M.both){
      if(z[i]==1){
        ll.y[i,]=dnbinom(y.true2D[i,],mu=lamd[i,],size=theta.d*K1D,log=TRUE)
      }
    }
  }else{
    stop("obstype not recognized")
  }
  
  if(!is.finite(sum(ll.y)))stop("Starting observation model likelihood not finite. Possible error in K1D (if supplied by user) or problem initializing data.")
  
  
  return(list(s=s,z=z,ID=ID,y.full=y.true2D,y.event=y.true3D,K1D=K1D,
         n.samples=n.samples,n.fixed=n.fixed,samp.type=G.type,this.j=this.j,match=match,
         xlim=xlim,ylim=ylim,locs=locs,tel.inds=tel.inds,n.locs.ind=n.locs.ind,
         G.true=G.true,G.obs=G.obs,G.latent=G.latent,G.latent.marked=G.latent.marked,
         gammaMat=gammaMat))

}