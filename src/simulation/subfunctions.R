samp.slct = function(seed, fnt.pop, n, Cluster=NULL, Clt.samp=NULL, dsgn, size = NULL, size.I = NULL){
	set.seed(seed)
	N = nrow(fnt.pop)
	if(!is.null(Cluster)) size.Cluster = N/Cluster
	# one-ste sample design
	if(dsgn=="pps"){
	  fnt.pop$x=size
	  samp = sam.pps(fnt.pop,size, n)
	}
	# two-stage cluster sampling design (informative design at the second stage)
	if(dsgn == "srs-pps"){
	  
	  # calcualte the size for the second stage
	  fnt.pop$x = size              
	  
	  #-- first stage: select clusters by srs
	  # Clt.samp = 25
	  # n = 1000
	  index.psuI = sample(1:Cluster, Clt.samp, replac = F)
	  index.psuI = sort(index.psuI)  #sort selected psus
	  sample.I = fnt.pop[fnt.pop$psu %in% index.psuI,]
	  sample.I$wt.I = Cluster/Clt.samp
	  
	  #-- second stage: select subj.samp within selected psus by pps to p^a (or (1-p)^a)
	  samp=NULL
	  for (i in 1: Clt.samp){
	    popn.psu.i= sample.I[sample.I$psu==index.psuI[i, 1],]
	    size.II.i = sample.I[sample.I$psu==index.psuI[i, 1],"x"]
	    samp.i = sam.pps(popn.psu.i,size.II.i, n/Clt.samp)
	    samp.i$wt = samp.i$wt*samp.i$wt.I
	    samp = rbind(samp,samp.i)
	  }#sum(samp.cntl$wt);nrow(fnt.cntl)
	}
	if(dsgn == "pps-pps"){
	  
	  # calcualte the size for the second stage
	  fnt.pop$x = size               
	  
	  #-- first stage: select clusters by pps to size aggrated level p^a(or (1-p)^a)
	  # Clt.samp = 25
	  # n = 1000
	  index.psuI = sam.pps(matrix(1:Cluster,,1),size.I, Clt.samp)
	  index.psuI = index.psuI[order(index.psuI[,1]),]  #sort selected psus
	  sample.I = fnt.pop[fnt.pop$psu %in% index.psuI[,1],]
	  sample.I$wt.I = rep(index.psuI[,'wt'],each=size.Cluster)
	  
	  #-- second stage: select subj.samp within selected psus by pps to p^a (or (1-p)^a)
	  samp=NULL
	  for (i in 1: Clt.samp){
	    popn.psu.i= sample.I[sample.I$psu==index.psuI[i,1],]
	    size.II.i = sample.I[sample.I$psu==index.psuI[i,1],"x"]
	    samp.i = sam.pps(popn.psu.i,size.II.i, n/Clt.samp)
	    samp.i$wt = samp.i$wt*samp.i$wt.I
	    samp = rbind(samp,samp.i)
	  }#sum(samp.cntl$wt);nrow(fnt.cntl)
	}
	rownames(samp) = as.character(1:dim(samp)[1])
  return(samp)      
}

#################################################################################
#-FUNCTION sam.pps is to select a sample with pps sampling                      #
# INPUT:  popul - the population including response and covariates              #
#         MSize - the size of measure to compute the selection probabilities(>0)#
#         n -     the sample size                                               #
# OUTPUT: sample selected with pps sampling of size n,                          #
#         including response, covariates, and sample weights                    #
#################################################################################
sam.pps<-function(popul,Msize, n){
  N=nrow(popul)
  pps.samID=sample(N,n,replace=F,prob=Msize)   #F but make sure N is large!
  if (dim(popul)[2] == 1){
  	sam.pps.data=as.data.frame(popul[pps.samID,])
  	names(sam.pps.data) = names(popul)
  }else{sam.pps.data=popul[pps.samID,]}
  sam.pps.data$wt=sum(Msize)/n/Msize[pps.samID]               
  return(sam.pps = sam.pps.data)
}