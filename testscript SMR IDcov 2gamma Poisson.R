#This version allows for categorical ID covariates observable with or without individual ID of the sample.
#These covariates provide ID exclusions for the latent ID samples, improving precision.
#These covariates can also be used to model other parameters. See other testscripts for how to do this.
#Not done here.
#This version also estimates separate population frequencies for the ID covs for marked and unmarked.
#BUT, I haven't included a data simulator--using one that assumes they are the same. This model
#likely only makes sense for "premarked", not "natural" marks because in the latter, the different
#observed gamma frequencies are just the product of the capture process.

#This is an SMR data simulator and MCMC sampler that handles all sample types
#1) marked, known ID
#2) marked, unknown ID
#3) unmarked, unknown ID
#4) unknown marked status, unknown ID

#It handles both "premarked" scenarios where you know the number of marked individuals
#and "natural" where marks are from natural patterns on the animals so the number of marked 
#individuals is unknown. For "premarked", consider using the random thinning model with partial ID covariates.

#This sampler also handles telemetry for marked individuals in the "premarked" scenario. Don't try
#to use telemetry with "natural" scenario. Not realistic and I'm not sure what my code will do!

#y[i,j,k] ~ Poisson(lam[i,j,k])
#y.event[i,j,k,1:3] ~ Multinomial(theta.marked[1:3],y[i,j,k]) for marked i
#y.event[i,j,k,1:3] ~ Multinomial(theta.unmarked[1:3],y[i,j,k]) for unmarked i

#event 1 is you know the ID (marked known ID samples)
#event 2 is you know the mark status, but not ID (marked, unknown ID or unmarked samples)
#event 3 is you don't know mark status or ID (unknown marked status samples)

library(nimble)
source("sim.SMR.IDcov.R")
source("NimbleModel SMR IDcov 2gamma Poisson.R")
source("NimbleFunctions SMR IDcov 2gamma Poisson.R") #Need this version for correct custom GSampler
source("init.SMR.IDcov.R")
source("sSampler.R")

#If using Nimble version 0.13.1 and you must run this line 
nimbleOptions(determinePredictiveNodesInModel = FALSE)
# #If using Nimble before version 0.13.1, run this line instead
# nimble:::setNimbleOption('MCMCjointlySamplePredictiveBranches', FALSE)

####Simulate some data####
N=78
n.marked=12
lam0=0.5
sigma=0.5
K=10 #number of occasions
buff=3 #state space buffer
X<- expand.grid(3:11,3:11) #make a trapping array
#theta is probability of observing each sample type for marked and unmarked individuals
theta.marked=c(0.75,0.15,0.1) #P(ID, Marked no ID, unk status). must sum to 1
theta.unmarked=0.75 #prob known marked status. #P(ID, Marked no ID, unk status)=(0,theta.unmarked,1-theta.unmarked)
marktype="premarked" #are individuals premarked, or naturally marked?
# marktype="natural"
obstype="poisson"
tlocs=0 #number of telemetry locs/marked individual. For "premarked"
#categorical ID covariate stuff
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


data=sim.SMR.IDcov(N=N,n.marked=n.marked,marktype=marktype,
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

#note, the sample-level data is this: the type, the trap of capture, the observed ID covs
tmp=data.frame(samp.type=data$samp.type,this.j=data$this.j,G.obs=data$G.obs)
head(tmp)

####Fit model in Nimble####
if(marktype=="natural"){
  M1=40 #Augmentation level for marked.
}else{
  M1=n.marked #Set to n.marked if premarked. psi1 will be estimated, but can be ignored.
}
M2=125 #Augmentation level for unmarked
#Monitor N.M and N.UM, marked and unmarked ind abundance to make sure N.M does not hit M1
#and N.UM does not hit M1+M2 during sampling. If so, raise the offending M and run again.
M.both=M1+M2
#Need some inits to initialize data
#Use reasonable inits for lam0 and sigma since we check to make sure initial observation
#model likelihood is finite
#also use this function checks to make sure theta.marked and theta.unmarked inits are in
#the correct structure. 
#Initializing both marked and unmarked inds with same gamma
inits=list(lam0=lam0,sigma=sigma,gamma=gamma)

#This function structures the simulated data to fit the model in Nimble (some more restructing below)
#Also checks some inits
nimbuild=init.SMR.IDcov(data,inits,M1=M1,M2=M2,marktype=marktype,obstype="poisson")

#We include marked/unmarked status as the first ID category for the Nimble sampler, so add 1
n.cat.nim=n.cat+1

#we are going to stuff the gamma values into a ragged matrix for use in Nimble.
#Also, we are going to include gamma values for marked status in the first row, but we are
#not actually going to use them. This makes Nimble happy--I tried not including these, but 
#Nimble changes the inits I give it for marked status. By not providing a prior (in the model code) on these
#probabilities, they are not updated
gammaMat=rbind(c(0.5,0.5),nimbuild$gammaMat)
n.levels.nim=c(2,n.levels)

#make G.true data and initial values. First column of G.true is the marked status, which is known
#for all individuals
G.true.init=nimbuild$G.true
G.true.init[,1]=NA #Column 1 is data
G.true.data=nimbuild$G.true
G.true.data[,2:n.cat.nim]=NA #columns 2:n.cat.nim are latent. Fixed values for marked inds not updated during MCMC

#inits for nimble
#full inits. Nimble can initialize psi1 and psi2, but if sigma and lam0 initialized too far away
#from truth, it can stop adapting before convergence and mix very poorly.
Niminits <- list(z=nimbuild$z,s=nimbuild$s,G.true=G.true.init,ID=nimbuild$ID,capcounts=rowSums(nimbuild$y.full),
                 y.full=nimbuild$y.full,y.event=nimbuild$y.event,
                 gammaMat.M=gammaMat,gammaMat.UM=gammaMat, #plugging in 1 gammaMat for marked and unmarked gammaMat inits
                 theta.unmarked=c(0,0.5,0.5),G.latent=nimbuild$G.latent,
                 lam0=inits$lam0,sigma=inits$sigma)

#constants for Nimble
J=nrow(data$X)
constants<-list(M1=M1,M2=M2,M.both=M.both,J=J,K=K,K1D=data$K1D,n.samples=nimbuild$n.samples,
                n.cat=n.cat.nim,n.levels=n.levels.nim,xlim=data$xlim,ylim=data$ylim)

# Supply data to Nimble. Note, y.true and y.true.event are treated as completely latent (but known IDs enforced)
z.data=c(rep(1,data$n.marked),rep(NA,M.both-data$n.marked))

Nimdata<-list(y.full=matrix(NA,nrow=M.both,ncol=J),y.event=array(NA,c(M.both,J,3)),
              G.true=G.true.data,ID=rep(NA,nimbuild$n.samples),z=z.data,X=as.matrix(X),capcounts=rep(NA,M.both))

# #If you have telemetry use these instead. Make sure to uncomment telemetry BUGS code.
# constants<-list(M1=M1,M2=M2,M.both=M.both,J=J,K=K,K1D=data$K1D,n.samples=nimbuild$n.samples,
#                 n.cat=n.cat.nim,n.levels=n.levels.nim,xlim=data$xlim,ylim=data$ylim,tel.inds=nimbuild$tel.inds,
#                 n.tel.inds=length(nimbuild$tel.inds),n.locs.ind=nimbuild$n.locs.ind)
# Nimdata<-list(y.full=matrix(NA,nrow=M.both,ncol=J),y.event=array(NA,c(M.both,J,3)),
#               G.true=G.true.data,ID=rep(NA,nimbuild$n.samples),z=z.data,X=as.matrix(X),capcounts=rep(NA,M.both),
#               locs=data$locs)

# set parameters to monitor
parameters=c('psi1','psi2','lam0','sigma','theta.marked','theta.unmarked','gammaMat.M',
              'gammaMat.UM','n.M','n.UM','N.M','N.UM','N.tot')
#other things we can monitor with separate thinning rate
parameters2=c("ID","s")

# Build the model, configure the mcmc, and compile
start.time<-Sys.time()
Rmodel <- nimbleModel(code=NimModel, constants=constants, data=Nimdata,check=FALSE,
                      inits=Niminits)
conf <- configureMCMC(Rmodel,monitors=parameters, thin=1,
                      monitors2=parameters2, thin2=10,
                      useConjugacy = TRUE)

###Two *required* sampler replacement
##Here, we remove the default samplers for y.full and y.event, which are not correct
#and replace it with the custom "IDSampler"
conf$removeSampler("y.full")
conf$removeSampler("y.event")
conf$addSampler(target = paste0("y.full[1:",M.both,",1:",J,"]"),
                type = 'IDSampler',control = list(M1=M1,M2=M2,M.both=M.both,J=J,K1D=data$K1D,
                                                  n.fixed=nimbuild$n.fixed,n.marked=n.marked,
                                                  samp.type=nimbuild$samp.type,n.samples=nimbuild$n.samples,
                                                  this.j=nimbuild$this.j,n.cat=n.cat.nim,G.obs=nimbuild$G.obs,
                                                  G.latent.marked=nimbuild$G.latent.marked),
                silent = TRUE)
#replace default G.true sampler, which is not correct, with custom sampler for G.true, "GSampler"
conf$removeSampler("G.true")
for(i in 1:M.both){
  for(m in 2:n.cat.nim){ #don't need to update first cat bc it is mark status
    conf$addSampler(target = paste("G.true[",i,",",m,"]", sep=""),
                    type = 'GSampler',
                    control = list(i = i,m=m,M1=M1,M.both=M.both,n.cat=n.cat.nim,n.samples=nimbuild$n.samples,
                                   n.levels=n.levels.nim,G.obs=nimbuild$G.obs), silent = TRUE) 
  }
}

###Two *optional* sampler replacements:

#replace default activity center sampler that updates x and y locations separately with a joint update
#a little more efficient. sSampler below only tunes s when z=1. Should better tune activity centers for 
#uncaptured individuals
conf$removeSampler(paste("s[1:",M.both,", 1:2]", sep=""))
for(i in 1:M.both){
  # conf$addSampler(target = paste("s[",i,", 1:2]", sep=""), #do not adapt covariance bc s's not deterministically linked to unmarked individuals
  #                 type = 'RW_block',control=list(adaptive=TRUE,adaptScaleOnly=TRUE,adaptInterval=500),silent = TRUE)
  conf$addSampler(target = paste("s[",i,", 1:2]", sep=""),
                  type = 'sSampler',control=list(i=i,xlim=nimbuild$xlim,ylim=nimbuild$ylim,scale=0.25),silent = TRUE)
  #scale parameter here is just the starting scale. It will be tuned.
}

#replace independent lam0 and sigma samplers with block sampler better accommodating for posterior covariance
#should improve mixing and increase posterior effective sample size. AF_slice works better than block RW. 
#Need to not use this update or modify it when using lam0 or sigma covariates.
#This sampler is slower, so not worth it if data is not so sparse there is strong posterior correlation
#between lam0 and sigma
conf$removeSampler(c("lam0","sigma"))
conf$addSampler(target = c("lam0","sigma"),type = 'AF_slice',
                control = list(adaptive=TRUE),silent = TRUE)

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

data$n.M #true number of captured marked individuals
data$n.UM #true number of captured unmarked individuals

#Important! If N.UM hits M2 during sampling, raise M2. 
#For an unknown number of marked individuals, if N.M hits M1 during sampling, raise M1.


####Look an posterior pairwise sample match probs
#assuming you monitored ID in 2nd monitor
library(MCMCglmm)
mvSamples2 = as.matrix(Cmcmc$mvSamples2)
idx=grep("ID",colnames(mvSamples2))
burnin=10 #set appropriately...
IDpost=posterior.mode(mcmc(mvSamples2[burnin:nrow(mvSamples2),idx]))
#For simulated data sets, comparing posterior mode ID to truth.
#Numbers will not be the same (except marked individuals), but all samples with same true ID will have
#same ID in posterior mode when posterior mode is exactly correct. Numbers just don't match up.
cbind(data$ID,round(IDpost))

#calculate posterior probability of pairwise sample matches
#P(sample x belongs to same individual as sample y)
n.samples=length(data$this.j)
n.iter=nrow(mvSamples2)
pair.probs=matrix(NA,n.samples,n.samples)
for(i in 1:n.samples){
  for(j in 1:n.samples){
    count=0
    for(iter in burnin:n.iter){
      count=count+1*(mvSamples2[iter,idx[j]]==mvSamples2[iter,idx[i]])
    }
    pair.probs[i,j]=count/(n.iter-burnin+1)
  }
}

this.samp=1 #sample number to look at
round(pair.probs[this.samp,],3) #probability this sample is from same individual as all other samples
round(pair.probs[this.samp,data$ID==data$ID[this.samp]],3) #for simulated data, these are the other samples truly from same individual

#marked but no ID samples will generally be linked with correct samples with higher prob than unmarked.
#most uncertainty in unknown marked status


#Can look at activity centers. Can help make sure telmetry is formatted correctly.
#telemetry individuals should have more precise AC ests, on average.
burnin=20
mvSamples2 = as.matrix(Cmcmc$mvSamples2)
idx=grep("s",colnames(mvSamples2))
#look at chains
plot(mcmc(mvSamples2[burnin:nrow(mvSamples2),idx]))

#plot at 1 at a time
par(mfrow=c(1,1),ask=FALSE)
this.i=1
plot(mvSamples2[burnin:nrow(mvSamples2),idx[this.i]],
     mvSamples2[burnin:nrow(mvSamples2),idx[this.i+M.both]],xlim=data$xlim,ylim=data$ylim)
points(data$X,pch=4)

