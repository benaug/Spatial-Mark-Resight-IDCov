NimModel <- nimbleCode({
  #detection function priors - shared across sessions
  lam0.fixed ~ dunif(0,15)
  sigma.fixed ~ dunif(0,10)
  theta.d.fixed ~ dunif(0,25) #careful with this prior. Too much prior mass near 0 gives very strong prior weight to high overdispersion
  #Expected density/data augmentation priors for marked + unmarked individuals
  D ~ dunif(0,10) #Expected density
  
  for(g in 1:N.session){
    #plug in shared df parameter for each session. Must use lam0[g] and sigma[g] here for custom update.
    #alternatively, can be estimated separately or with random effects.
    lam0[g] <- lam0.fixed
    sigma[g] <- sigma.fixed
    theta.d[g] <- theta.d.fixed
    lambda[g] <- D*area[g] #expected N
    N[g] ~ dpois(lambda[g]) #realized N
    
    #sample type observation model priors (Dirichlet)
    alpha.marked[g,1] <- 1
    alpha.marked[g,2] <- 1
    alpha.marked[g,3] <- 1
    alpha.unmarked[g,1] <- 1
    alpha.unmarked[g,2] <- 1
    theta.marked[g,1:3] ~ ddirch(alpha.marked[g,1:3])
    theta.unmarked[g,1] <- 0
    theta.unmarked[g,2:3] ~ ddirch(alpha.unmarked[g,1:2])
    
    #categorical ID covariate priors
    for(m in 2:n.cat[g]){ #skip first cat for marked/unmarked
      for(l in 1:n.levels[g,m]){
        alpha[g,m,l] <- 1 # prior parameters
      }
      gammaMat[g,m,1:n.levels[g,m]] ~ ddirch(alpha[g,m,1:n.levels[g,m]])
    }
    
    #likelihoods (except for s/z priors)
    #Marked individuals first
    for(i in 1:M1[g]) {
      for(m in 1:n.cat[g]){
        G.true[g,i,m] ~ dcat(gammaMat[g,m,1:n.levels[g,m]])
      }
      s[g,i,1] ~ dunif(xlim[g,1],xlim[g,2])
      s[g,i,2] ~ dunif(ylim[g,1],ylim[g,2])
      lam[g,i,1:J[g]] <- GetDetectionRate(s = s[g,i,1:2], X = X[g,1:J[g],1:2], J=J[g],sigma=sigma[g], lam0=lam0[g], z=z[g,i])
      p[g,i,1:J[g]] <- theta.d[g]/(theta.d[g]+lam[g,i,1:J[g]])
      y.full[g,i,1:J[g]] ~ dNBVector(p=p[g,i,1:J[g]],theta.d=theta.d[g]*K1D[g,1:J[g]],z=z[g,i]) #vectorized obs mod. trap op: sum of N NB RVs is NB with theta.d=N*theta.d
      #custom distribution that skips likelihood eval for the individuals currently with 0 captures.
      y.event[g,i,1:J[g],1:3] ~ dmulti2(y.full[g,i,1:J[g]],prob=theta.marked[g,1:3],capcounts=capcounts[g,i])
    }
    
    #Then unmarked individuals
    for(i in (M1[g]+1):M.both[g]){
      for(m in 1:n.cat[g]){
        G.true[g,i,m] ~ dcat(gammaMat[g,m,1:n.levels[g,m]])
      }
      s[g,i,1] ~ dunif(xlim[g,1],xlim[g,2])
      s[g,i,2] ~ dunif(ylim[g,1],ylim[g,2])
      lam[g,i,1:J[g]] <- GetDetectionRate(s = s[g,i,1:2], X = X[g,1:J[g],1:2], J=J[g],sigma=sigma[g], lam0=lam0[g], z=z[g,i])
      p[g,i,1:J[g]] <- theta.d[g]/(theta.d[g]+lam[g,i,1:J[g]])
      y.full[g,i,1:J[g]] ~ dNBVector(p=p[g,i,1:J[g]],theta.d=theta.d[g]*K1D[g,1:J[g]],z=z[g,i]) #vectorized obs mod. trap op: sum of N NB RVs is NB with theta.d=N*theta.d
      #custom distribution that skips likelihood eval for the individuals currently with 0 captures.
      y.event[g,i,1:J[g],2:3] ~ dmulti2(y.full[g,i,1:J[g]],prob=theta.unmarked[g,2:3],capcounts=capcounts[g,i])
    }
    
    # #If you have telemetry
    # for(i in 1:n.tel.inds[g]){
    #   for(m in 1:n.locs.ind[g,i]){
    #     locs[g,tel.inds[g,i],m,1] ~ dnorm(s[g,tel.inds[g,i],1],sd=sigma[g])
    #     locs[g,tel.inds[g,i],m,2] ~ dnorm(s[g,tel.inds[g,i],2],sd=sigma[g])
    #   }
    # }
    
    #calculate number of marked and unmarked inds captured and abundance
    capcounts[g,1:M.both[g]] <- Getcapcounts(y.full=y.full[g,1:M.both[g],1:J[g]])
    n.M[g] <- Getncap(capcounts=capcounts[g,1:M1[g]],ID=ID[g,1:n.samples[g]],G.latent=G.latent[g,1:M.both[g],1:n.cat[g]])
    n.UM[g] <- Getncap(capcounts=capcounts[g,(M1[g]+1):M.both[g]],ID=ID[g,1:n.samples[g]],G.latent=G.latent[g,1:M.both[g],1:n.cat[g]])
    N.UM[g] <- N[g] - M1[g]
  }
})# end model
