#Multisession SMR using reversible jump MCMC instead of data augmentation
#For premarked individuals only. I will write a slightly modified sampler for natural marks.

#This testscript shows how to share lam0, sigma, and/or expected density across sessions
#can share theta.marked and theta.unmarked but not done here.

library(nimble)
source("sim.SMR.IDcov.DF.multisession.R")
source("sim.SMR.IDcov.DF.R")
source("NimbleModel SMR IDcov DF Multisession Premarked Poisson.R")
source("NimbleFunctions SMR IDcov Multisession Premarked Poisson.R") #no extra DF stuff here
source("init.SMR.IDcov.DF.multisession.R")
source("init.SMR.IDcov.DF.R")
source("sSampler Multisession.R")

#If using Nimble version 0.13.1 and you must run this line 
nimbleOptions(determinePredictiveNodesInModel = FALSE)
# #If using Nimble before version 0.13.1, run this line instead
# nimble:::setNimbleOption('MCMCjointlySamplePredictiveBranches', FALSE)

####Simulate some data####
#Here, I'll simulate 3 populations with different n.marked, K, X, and state space areas
#sharing D, lam0, sigma so they can be shared during estimation
N.session=3
D = rep(0.4,N.session) #expected density in units of sigma and X
n.marked=c(12,13,14)

#here, we have N.session x 2 group detection functions. 
#Could be sex-specific parms shared across sessions.
lam0=cbind(rep(0.5,N.session),rep(0.25,N.session))
sigma=cbind(rep(0.5,N.session),rep(0.75,N.session))

K=c(5,6,7) #number of occasions
buff=rep(2,N.session) #state space buffer
#make trapping arrays
X1=expand.grid(3:11,3:11)
X2=expand.grid(3:12,3:12)
X3=expand.grid(3:13,3:13)
X=list(X1,X2,X3) #put in a list, one for each session

#See what expected N is for these expected D and state space areas
area=getArea(X=X,buff=buff)
area #state space areas for each session resulting from X and buff
lambda=D*area
lambda #expected N in each session

#theta is probability of observing each sample type for marked and unmarked individuals
theta.marked=matrix(rep(c(0.75,0.15,0.1),N.session),nrow=N.session,byrow=TRUE) #P(ID, Marked no ID, unk status). must sum to 1
theta.unmarked=rep(0.75,N.session) #prob known marked status. #P(ID, Marked no ID, unk status)=(0,theta.unmarked,1-theta.unmar
marktype="premarked" #are individuals premarked, or naturally marked? This test script only handles premarked.
# marktype="natural"
obstype="poisson"
tlocs=c(0,0,0) #number of telemetry locs/marked individual in each session. For "premarked"

#categorical ID covariate stuff - not session-specific in data simulator
n.cat=2  #number of ID categories (not including marked status)
gamma=IDcovs=vector("list",n.cat) #population frequencies of each category level. Assume equal here.
n.levels=rep(2,n.cat) #number of levels per IDcat
if(all(n.levels==1))stop("This specification has no categorial ID covariates. Use testscript for regular SMR.")
for(i in 1:n.cat){
  gamma[[i]]=rep(1/n.levels[i],n.levels[i])
  IDcovs[[i]]=1:n.levels[i]
}
theta.cat=rep(1,n.cat)#sample-level IDcov observation probabilities. Data missing at random if <1. 
#data simulator assumes all IDcovs known for marked inds. MCMC sampler accepts missing values coded as 0.
data=sim.SMR.IDcov.DF.multisession(N.session=N.session,lambda=lambda,n.marked=n.marked,marktype=marktype,
             theta.marked=theta.marked,theta.unmarked=theta.unmarked,
             lam0=lam0,sigma=sigma,K=K,X=X,buff=buff,tlocs=tlocs,
             obstype=obstype,
             n.cat=n.cat,IDcovs=IDcovs,gamma=gamma,theta.cat=theta.cat)

#What is the observed data?
head(data$samp.type) #vector of each samp.type of each detection. Must be in this order (I think).
table(data$samp.type)
head(data$this.j) #trap of capture for each sample
head(data$this.k) #occasion of each capture for each sample (not used in this "2D" sampler)
head(data$ID.marked) #true ID's for marked and identified samples
#here are observed ID covariate data. Missing values coded with "0" instead of "NA"
dim(data$G.marked) #ID covs for marked individuals (can be missing, coded with 0, though not simulated) n.marked x n.cat
dim(data$G.obs) #ID covs for each sample. n.samples x n.cat
str(data$locs) #possibly telemetry. n.marked x tlocs x 2 array (or ragged array if number of locs/ind differ). 
#Rows are 1:n.marked individuals, columns are max telemetry points for a single
#individual, fill in NAs for inds with no telemetry and/or inds without max number of telemetry points.
#in latter case, order telemetry points first, then NAs

####Fit model in Nimble####
if(marktype=="natural"){
  # M1=40 #Augmentation level for marked.
  stop("Natural marks not handled with this testscript. There will be another one that does.")
}else{
  M1=n.marked #Set to n.marked if premarked. psi1 will be estimated, but can be ignored.
}
M2=c(155,165,175) #Augmentation level for unmarked
#Monitor N.M and N.UM, marked and unmarked ind abundance to make sure N.M does not hit M1
#and N.UM does not hit M1+M2 during sampling. If so, raise the offending M and run again.
M.both=M1+M2
#Need some inits to initialize data
#Use reasonable inits for lam0 and sigma since we check to make sure initial observation
#model likelihood is finite
#also use this function checks to make sure theta.marked and theta.unmarked inits are in
#the correct structure. 
inits=list(lam0=lam0,sigma=sigma) #initializing with 1 parameter per session, just set all to same value

#This function structures the simulated data to fit the model in Nimble (some more restructing below)
#Also checks some inits
nimbuild=init.SMR.IDcov.DF.multisession(data,inits,M1=M1,M2=M2,marktype=marktype,obstype="poisson")

#We include marked/unmarked status as the first ID category for the Nimble sampler, so add 1
n.cat.nim=nimbuild$n.cat+1

#we are going to stuff the gamma values into a ragged matrix for use in Nimble.
#Also, we are going to include gamma values for marked status in the first row, but we are
#not actually going to use them. This makes Nimble happy--I tried not including these, but 
#Nimble changes the inits I give it for marked status. By not providing a prior (in the model code) on these
#probabilities, they are not updated
n.levels.nim=cbind(2,nimbuild$n.levels)

#make G.true data and initial values. First column of G.true is the marked status, which is known
#for all individuals
G.true.init=nimbuild$G.true
G.true.init[,,1]=NA #Column 1 is data
G.true.data=nimbuild$G.true
G.true.data[,,2:max(n.cat.nim)]=NA #columns 2:n.cat.nim are latent. Fixed values for marked inds not updated during MCMC

#inits for nimble
theta.unmarked.init=matrix(c(0,0.5,0.5),N.session,3,byrow=TRUE)
N.init=rowSums(nimbuild$z,na.rm=TRUE)
N.UM.init=rep(NA,N.session)
for(g in 1:N.session){
  N.UM.init[g]=sum(nimbuild$z[g,(M1[g]+1):M.both[g]])
}
(N.init-N.UM.init)==n.marked #should be n.marked[g] individuals in initialized data
library(abind)
gammaMat.init=array(NA,dim=c(N.session,max(n.cat.nim),max(nimbuild$n.levels)))
gammaMat.init[,-1,]=nimbuild$gammaMat
gammaMat.init[,1,1:2]=0.5

Niminits <- list(N=N.init,N.UM=N.UM.init,
                 z=nimbuild$z,s=nimbuild$s,G.true=G.true.init,ID=nimbuild$ID,capcounts=apply(nimbuild$y.full,c(1,2),sum),
                 y.full=nimbuild$y.full,y.event=nimbuild$y.event,
                 gammaMat=gammaMat.init,theta.unmarked=theta.unmarked.init,G.latent=nimbuild$G.latent,
                 lam0.fixed=c(0.75,0.75),sigma.fixed=c(0.5,0.5),D=0.5)

#constants for Nimble
J=unlist(lapply(data$X,nrow))
constants<-list(N.session=N.session,M1=M1,M2=M2,M.both=M.both,J=J,K=K,K1D=nimbuild$K1D,n.samples=nimbuild$n.samples,
                n.cat=n.cat.nim,n.levels=n.levels.nim,xlim=data$xlim,ylim=data$ylim,area=area)

# Supply data to Nimble. Note, y.true and y.true.event are treated as completely latent (but known IDs enforced)
z.data=matrix(NA,N.session,max(M.both))
for(g in 1:N.session){
  z.data[g,1:data$n.marked[g]]=1
}

Nimdata<-list(y.full=array(NA,dim=c(N.session,max(M.both),max(J))),y.event=array(NA,c(N.session,max(M.both),max(J),3)),
              ID=matrix(NA,N.session,max(nimbuild$n.samples)),z=z.data,X=nimbuild$X,capcounts=matrix(NA,N.session,max(M.both)),
              G.true=G.true.data)

# #If you have telemetry use these instead. Make sure to uncomment telemetry BUGS code.
# constants<-list(N.session=N.session,M1=M1,M2=M2,M.both=M.both,J=J,K=K,K1D=nimbuild$K1D,n.samples=nimbuild$n.samples,
#                 n.cat=n.cat.nim,n.levels=n.levels.nim,xlim=data$xlim,ylim=data$ylim,area=area,
#                 tel.inds=nimbuild$tel.inds,n.tel.inds=nimbuild$n.tel.inds,n.locs.ind=nimbuild$n.locs.ind)
# Nimdata<-list(y.full=array(NA,dim=c(N.session,max(M.both),max(J))),y.event=array(NA,c(N.session,max(M.both),max(J),3)),
#               ID=matrix(NA,N.session,max(nimbuild$n.samples)),z=z.data,X=nimbuild$X,capcounts=matrix(NA,N.session,max(M.both)),
#               G.true=G.true.data,locs=data$locs)

# set parameters to monitor
parameters=c('D','lambda','lam0.fixed','sigma.fixed','theta.marked','theta.unmarked','gammaMat',
              'n.M','n.UM','N.UM','N')
#other things we can monitor with separate thinning rate
parameters2=c("ID","s")

# Build the model, configure the mcmc, and compile
start.time<-Sys.time()
Rmodel <- nimbleModel(code=NimModel, constants=constants, data=Nimdata,check=FALSE,
                      inits=Niminits)
conf <- configureMCMC(Rmodel,monitors=parameters, thin=1,
                      monitors2=parameters2, thin2=10,
                      useConjugacy = TRUE)

###Three *required* sampler replacement
##Here, we remove the default samplers for y.full and y.event, which are not correct
#and replace it with the custom "IDSampler"
conf$removeSampler("y.full")
conf$removeSampler("y.event")
n.samples=nimbuild$n.samples
for(g in 1:N.session){
  conf$addSampler(target = paste0("y.full[",g,",1:",M.both[g],",1:",J[g],"]"),
                  type = 'IDSampler',control = list(M1=M1[g],M2=M2[g],M.both=M.both[g],J=J[g],K1D=nimbuild$K1D[g,1:J[g]],
                                                    n.fixed=nimbuild$n.fixed[g],samp.type=nimbuild$samp.type[g,1:n.samples[g]],
                                                    n.samples=n.samples[g],
                                                    this.j=nimbuild$this.j[g,1:n.samples[g]],
                                                    match=nimbuild$match[g,1:n.samples[g],1:M.both[g]],
                                                    g=g,n.cat=n.cat.nim[g],G.obs=nimbuild$G.obs[g,1:n.samples[g],1:n.cat.nim[g]],
                                                    G.latent.marked=nimbuild$G.latent.marked[g,1:M1[g],1:n.cat.nim[g]]),
                  silent = TRUE)
}

#replace default G.true sampler, which is not correct, with custom sampler for G.true, "GSampler"
conf$removeSampler("G.true")
for(g in 1:N.session){
  for(i in 1:M.both[g]){
    for(m in 2:n.cat.nim[g]){ #don't need to update first cat bc it is mark status
      conf$addSampler(target = paste("G.true[",g,",",i,",",m,"]", sep=""),
                      type = 'GSampler',
                      control = list(g=g,i=i,m=m,n.levels=n.levels.nim[g,]), silent = TRUE) 
    }
  }
}

z.ups=round(M.both*0.25) # how many z proposals per iteration per session?
J=nimbuild$J
conf$removeSampler("N")
for(g in 1:N.session){
  #nodes used for update, calcNodes + z nodes
  y.nodes <- Rmodel$expandNodeNames(paste("y.full[",g,",","1:",M.both[g],",1:",J[g],"]"))
  lam.nodes <- Rmodel$expandNodeNames(paste("lam[",g,",","1:",M.both[g],",1:",J[g],"]"))
  N.node <- Rmodel$expandNodeNames(paste("N[",g,"]"))
  N.UM.node <- Rmodel$expandNodeNames(paste("N.UM[",g,"]")) #only used to update derived parameter when N updates
  z.nodes <- Rmodel$expandNodeNames(paste("z[",g,",","1:",M.both[g],"]"))
  calcNodes <- c(N.node,N.UM.node,y.nodes,lam.nodes)
  
  conf$addSampler(target = c("N"),
                  type = 'zSampler',control = list(z.ups=z.ups[g],J=J[g],M1=M1[g],M.both=M.both[g],g=g,
                                                   y.nodes=y.nodes,lam.nodes=lam.nodes,N.node=N.node,z.nodes=z.nodes,
                                                   calcNodes=calcNodes),
                  silent = TRUE)
}

###Two *optional* sampler replacements:

#replace default activity center sampler that updates x and y locations separately with a joint update
#a little more efficient. sSampler below only tunes s when z=1. Should better tune activity centers for 
#uncaptured individuals
conf$removeSampler("s")
for(g in 1:N.session){
  for(i in 1:M.both[g]){
    conf$addSampler(target = paste("s[",g,",",i,", 1:2]", sep=""),
                    type = 'sSampler',control=list(i=i,g=g,xlim=nimbuild$xlim[g,],ylim=nimbuild$ylim[g,],scale=1),silent = TRUE)
    #scale parameter here is just the starting scale. It will be tuned.
  }
}

#replace independent lam0 and sigma samplers with block sampler better accommodating for posterior covariance
#should improve mixing and increase posterior effective sample size. AF_slice works better than block RW. 
#Need to not use this update or modify it when using lam0 or sigma covariates.
#This sampler is slower, so not worth it if data is not so sparse there is strong posterior correlation
#between lam0 and sigma
#looping over number of levels of 1st ID cov (2nd in nimble due to marked status being 1)
#assuming n.levels is the same across sessions for this IDcov
conf$removeSampler(c("lam0.fixed","sigma.fixed"))
for(df in 1:n.levels.nim[1,2]){
  conf$addSampler(target = c(paste("lam0.fixed[",df,"]"),paste("sigma.fixed[",df,"]")),type = 'RW_block',
                  control = list(adaptive=TRUE),silent = TRUE)
}

# Build and compile
Rmcmc <- buildMCMC(conf)
# runMCMC(Rmcmc,niter=1) #this will run in R, used for better debugging
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

# Run the model.
start.time2<-Sys.time()
Cmcmc$run(2500,reset=FALSE) #short run for demonstration. can keep running this line to get more samples
end.time<-Sys.time()
end.time-start.time  # total time for compilation, replacing samplers, and fitting
end.time-start.time2 # post-compilation run time

library(coda)
mvSamples = as.matrix(Cmcmc$mvSamples)
plot(mcmc(mvSamples[2:nrow(mvSamples),]))

data$N #realized Ns
lambda #expected N
data$n.M #true number of captured marked individuals
data$n.UM #true number of captured unmarked individuals

#Important! If N.UM hits M2 during sampling, raise M2. 