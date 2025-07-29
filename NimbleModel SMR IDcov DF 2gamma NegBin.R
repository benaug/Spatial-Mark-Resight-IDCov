NimModel <- nimbleCode({
  #detection function priors as function of first ID cov (excluding marked status, which is 1st)
  for(i in 1:n.levels[2]){
    lam0[i] ~ dunif(0,15)
    sigma[i] ~ dunif(0,10)
  }
  theta.d ~ dunif(0,25) #careful with this prior. Too much prior mass near 0 gives very strong prior weight to high overdispersion
  #data augmentation priors for marked (1) and unmarked (2) individuals
  psi1 ~ dunif(0,1)
  psi2 ~ dunif(0,1)
  #sample type observation model priors (Dirichlet)
  alpha.marked[1] <- 1
  alpha.marked[2] <- 1
  alpha.marked[3] <- 1
  alpha.unmarked[1] <- 1
  alpha.unmarked[2] <- 1
  theta.marked[1:3] ~ ddirch(alpha.marked[1:3])
  theta.unmarked[1] <- 0
  theta.unmarked[2:3] ~ ddirch(alpha.unmarked[1:2])
  
  #categorical ID covariate priors
  for(m in 2:n.cat){ #skip first cat for marked/unmarked
    alpha[m,1:n.levels[m]] <- 1 # parameters
    gammaMat.M[m,1:n.levels[m]] ~ ddirch(alpha[m,1:n.levels[m]])
    gammaMat.UM[m,1:n.levels[m]] ~ ddirch(alpha[m,1:n.levels[m]])
  }

  #likelihoods (except for s/z priors)
  #Marked individuals first
  for(i in 1:M1) {
    z[i] ~ dbern(psi1)
    for(m in 1:n.cat){
      G.true[i,m] ~ dcat(gammaMat.M[m,1:n.levels[m]])
    }
    s[i,1] ~ dunif(xlim[1],xlim[2])
    s[i,2] ~ dunif(ylim[1],ylim[2])
    #lam depends on G.true[i,2] here, 1st ID cov
    lam[i,1:J] <- GetDetectionRate(s = s[i,1:2], X = X[1:J,1:2], J=J,
                                   sigma=sigma[G.true[i,2]], lam0=lam0[G.true[i,2]], z=z[i]) 
    p[i,1:J] <- theta.d/(theta.d+lam[i,1:J])
    y.full[i,1:J] ~ dNBVector(p=p[i,1:J],theta.d=theta.d*K1D[1:J],z=z[i]) #vectorized obs mod. trap op: sum of N NB RVs is NB with theta.d=N*theta.d
    #custom distribution that skips likelihood eval for the individuals currently with 0 captures.
    y.event[i,1:J,1:3] ~ dmulti2(y.full[i,1:J],prob=theta.marked[1:3],capcounts=capcounts[i])
  }

  #Then unmarked individuals
  for(i in (M1+1):M.both){
    z[i] ~ dbern(psi2)
    for(m in 1:n.cat){
      G.true[i,m] ~ dcat(gammaMat.UM[m,1:n.levels[m]])
    }
    s[i,1] ~ dunif(xlim[1],xlim[2])
    s[i,2] ~ dunif(ylim[1],ylim[2])
    #lam depends on G.true[i,2] here, 1st ID cov
    lam[i,1:J] <- GetDetectionRate(s = s[i,1:2], X = X[1:J,1:2], J=J,
                                   sigma=sigma[G.true[i,2]], lam0=lam0[G.true[i,2]], z=z[i]) 
    p[i,1:J] <- theta.d/(theta.d+lam[i,1:J])
    y.full[i,1:J] ~ dNBVector(p=p[i,1:J],theta.d=theta.d*K1D[1:J],z=z[i]) #vectorized obs mod. trap op: sum of N NB RVs is NB with theta.d=N*theta.d
    #custom distribution that skips likelihood eval for the individuals currently with 0 captures.
    y.event[i,1:J,2:3] ~ dmulti2(y.full[i,1:J],prob=theta.unmarked[2:3],capcounts=capcounts[i])
  }
  
  #If you have telemetry
  for(i in 1:n.tel.inds){
    for(m in 1:n.locs.ind[i]){
      locs[i,m,1] ~ dnorm(s[tel.inds[i],1],sd=sigma[G.true[tel.inds[i],2]])
      locs[i,m,2] ~ dnorm(s[tel.inds[i],2],sd=sigma[G.true[tel.inds[i],2]])
    }
  }
  
  #calculate number of marked and unmarked inds captured and abundance
  capcounts[1:M.both] <- Getcapcounts(y.full=y.full[1:M.both,1:J])
  #stuffing G.latent here so we can use it in the model. not used in function
  n.M <- Getncap(capcounts=capcounts[1:M1],ID=ID[1:n.samples],G.latent=G.latent[1:M.both,1:n.cat])
  n.UM <- Getncap(capcounts=capcounts[(M1+1):M.both],ID=ID[1:n.samples],G.latent=G.latent[1:M.both,1:n.cat])
  N.M <- sum(z[1:M1])
  N.UM <- sum(z[(M1+1):M.both])
  N.tot<- N.M + N.UM
})# end model
