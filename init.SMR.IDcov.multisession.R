e2dist<- function (x, y){
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}

init.SMR.IDcov.multisession=function(data,inits=NA,M1=NA,M2=NA,marktype="premarked",obstype="poisson"){
  N.session=nrow(data$this.j)
  n.samples=rowSums(!is.na(data$this.j))
  init.session=vector("list",N.session)
  
  #split inits by session
  inits.use=vector("list",N.session)
  parms=names(inits)
  for(g in 1:N.session){
    inits.use[[g]]=vector("list",length(parms))
    names(inits.use[[g]])=parms
    inits.use[[g]]$lam0=inits$lam0[g]
    inits.use[[g]]$sigma=inits$sigma[g]
    if(obstype=="negbin"){
      inits.use[[g]]$theta.d=inits$theta.d[g]
    }
  }
  for(g in 1:N.session){
    #append gamma inits
    gamma=vector("list",data$IDlist[[g]]$n.cat) #population frequencies of each category level. Assume equal here.
    n.levels=unlist(lapply(data$IDlist[[g]]$IDcovs,length)) #number of levels per IDcat
    for(i in 1:data$IDlist[[g]]$n.cat){
      gamma[[i]]=rep(1/n.levels[i],n.levels[i])
    }
    inits.use[[g]]$gamma=gamma
  }
  
  
  #initialize sessions one by one
  anyTelemetry=FALSE
  for(g in 1:N.session){
    if(all(is.na(data$locs))){
      locs.use=NA
    }else if(all(is.na((data$locs[g,,,])))){
      locs.use=NA
    }else{
      tlocs.sess.max=max(rowSums(!is.na(data$locs[g,1:M1[g],,1])))
      locs.use=data$locs[g,1:M1[g],1:tlocs.sess.max,1:2]
      anyTelemetry=TRUE
    }
    data.use=list(this.j=data$this.j[g,1:n.samples[g]],this.k=data$this.k[g,1:n.samples[g]],samp.type=data$samp.type[g,1:n.samples[g]],
                  ID.marked=data$ID.marked[[g]],n.marked=M1[g],locs=locs.use,X=data$X[[g]],buff=data$buff[g],
                  K1D=data$K1D[[g]],G.marked=data$G.marked[g,1:data$n.marked[g],],G.obs=data$G.obs[g,1:n.samples[g],],
                  IDlist=data$IDlist[[g]])
    init.session[[g]]=init.SMR.IDcov(data.use,inits.use[[g]],M1=M1[g],M2=M2[g],marktype=marktype,obstype=obstype)
  }
  
  J=unlist(lapply(data$X,nrow))
  M.both=M1+M2
  maxM.both=max(M.both)
  s=array(NA,dim=c(N.session,maxM.both,2))
  z=matrix(NA,N.session,maxM.both)
  ID=matrix(NA,N.session,max(n.samples))
  y.full=array(NA,dim=c(N.session,maxM.both,max(J)))
  y.event=array(NA,dim=c(N.session,maxM.both,max(J),3))
  K1D=matrix(NA,N.session,max(J))
  n.fixed=rep(NA,N.session)
  samp.type=matrix(NA,N.session,max(n.samples))
  match=array(NA,dim=c(N.session,max(n.samples),maxM.both))
  if(anyTelemetry){
    tel.inds=matrix(NA,N.session,dim(data$locs)[2],dim(data$locs)[3])
    n.locs.ind=matrix(NA,N.session,dim(data$locs)[2],dim(data$locs)[3])
    n.tel.inds=rep(NA,N.session)
  }else{
    tel.inds=n.locs.ind=n.tel.inds=NA
  }
  
  for(g in 1:N.session){
    s[g,1:M.both[g],]=init.session[[g]]$s
    z[g,1:M.both[g]]=init.session[[g]]$z
    ID[g,1:n.samples[g]]=init.session[[g]]$ID
    y.full[g,1:M.both[g],1:J[g]]=init.session[[g]]$y.full
    y.event[g,1:M.both[g],1:J[g],]=init.session[[g]]$y.event
    K1D[g,1:J[g]]=init.session[[g]]$K1D
    samp.type[g,1:n.samples[g]]=init.session[[g]]$samp.type
    match[g,1:n.samples[g],1:M.both[g]]=init.session[[g]]$match
    n.fixed[g]=init.session[[g]]$n.fixed
    if(anyTelemetry){
      n.tel.inds[g]=sum(rowSums(!is.na(data$locs[g,,,1]))>0)
      tel.inds[g,1:n.tel.inds[g]]=init.session[[g]]$tel.inds
      n.locs.ind[g,1:n.tel.inds[g]]=init.session[[g]]$n.locs.ind
    }
  }
  #put X in ragged array
  X.new=array(NA,dim=c(N.session,max(J),2))
  for(g in 1:N.session){
    X.new[g,1:J[g],]=data$X[[g]]
  }
  
  #remove unused telemetry dimensions if not all marked individuals telemetered
  if(anyTelemetry){
    rem.idx=which(colSums(is.na(n.locs.ind))==N.session)
    if(length(rem.idx)>0){
      n.locs.ind=n.locs.ind[,-rem.idx]
      tel.inds=tel.inds[,-rem.idx]
    }
  }
  
  #IDcov stuff
  n.cats=rep(NA,N.session)
  max.n.levels=rep(NA,N.session)
  for(g in 1:N.session){
    n.cats[g]=data$IDlist[[g]]$n.cat
    max.n.levels[g]=max(unlist(lapply(data$IDlist[[g]]$IDcovs,max)))
  }
  max.n.cat=max(n.cats)
  n.levels=matrix(NA,N.session,max.n.cat)
  for(g in 1:N.session){
    n.levels[g,]=unlist(lapply(data$IDlist[[g]]$IDcovs,max))
  }
  max.n.levels.all=max(n.levels)
  G.true=G.latent=array(NA,dim=c(N.session,maxM.both,max.n.cat+1))
  G.obs=array(NA,dim=c(N.session,max(n.samples),max.n.cat+1))
  G.latent.marked=array(NA,dim=c(N.session,max(data$n.marked),max.n.cat+1))
  gammaMat=array(NA,dim=c(N.session,max.n.cat,max.n.levels.all))
  for(g in 1:N.session){
    G.true[g,1:M.both[g],1:(n.cats[g]+1)]=init.session[[g]]$G.true
    G.latent[g,1:M.both[g],1:(n.cats[g]+1)]=init.session[[g]]$G.latent
    G.obs[g,1:n.samples[g],1:(n.cats[g]+1)]=init.session[[g]]$G.obs
    G.latent.marked[g,1:data$n.marked[g],1:(n.cats[g]+1)]=init.session[[g]]$G.latent.marked
    gammaMat[g,1:n.cats[g],1:max.n.levels[g]]=init.session[[g]]$gammaMat
  }
  
  return(list(s=s,z=z,ID=ID,y.full=y.full,y.event=y.event,K1D=K1D,J=J,X=X.new,
              n.samples=n.samples,n.fixed=n.fixed,samp.type=samp.type,this.j=data$this.j,match=match,
              xlim=data$xlim,ylim=data$ylim,locs=locs.use,tel.inds=tel.inds,n.locs.ind=n.locs.ind,n.tel.inds=n.tel.inds,
              G.true=G.true,G.obs=G.obs,G.latent=G.latent,G.latent.marked=G.latent.marked,
              gammaMat=gammaMat,n.levels=n.levels,n.cat=n.cats))
  
  
}
