### If used, please cite:
#
#   Maynard et al. (2017) Diversity begets diversity in competition for space. Nature Ecology & Evolution
#
#######

#This code contains the main sampling and model fitting functions for the analytical, patch-occupany modeling results


########## This selects random communities based on a pairwise competition matrix, and calculates the steady-state abundance 

## set up for parallel processing, but doesn't save much time here since the internal functions in sim.comms are more computationally intensive (i.e., poorly coded)

registerDoMC(8) 
getDoParWorkers()
getDoParName()

min.com<-3
max.com<-9

sim.num<-1000

com.comb<-list()

#prob.vec is a vector of sample probabilities, in case you don't want all species to be sampled equally
#comp.df is the pairwise competition matrix, in same order as prob.vec


final<-foreach(k = min.com:max.com)%dopar%{

   com.comb[[k]]<-sim.comms(k,big.cmat=comp.df,snum=sim.num,pvec=prob.vec)

}



## get intransitivity, for example, for the community with 9 species
apply(com.comb[[9]]$namemat,1,calc.intran,big.comp=comp.df,allp=permutations(9,9))




#### MISC FUNCTIONS CALLED IN THE ABOVE CODE



## calculates the p_ij value for each cell in the transition matrix. calls 'getp' which does the heavy lifting
get.pmat<-function(cmat){

	omat<-matrix(NA,nrow=nrow(cmat),ncol=ncol(cmat))

	rownames(omat)<-colnames(omat)<-rownames(cmat)

	for(ti in 1:nrow(cmat)){
		for(tj in 1:ncol(cmat)){
			omat[ti,tj]<-getp(i=ti,j=tj,tind=seq(1:nrow(cmat)),cmat=cmat,rep(nrow(cmat),nrow(cmat)))
		}
	}

	omat
}

## iteratively calculates P_ij from Ulrich et al. 2014, Oikos, equations 8-9
getp<-function(i,j,tind,cmat,new.m){

	tsum<-0

	tm<-max(2,new.m[j])

	if(i==j){
		## probability a species replaces itself
		fbit<-1
		for(k in tind){
			if(k!=i & k!=j){
				fbit<-fbit*(1-cmat[k,i]) 
			}
		}
	}
	else{
		## probability i defeats j
		fbit<-1/(tm-1)*cmat[i,j]
		if(length(tind)>3){
			for(k in tind){
				if(k!=i & k!=j){
					tsum<-tsum+cmat[j,k]*getp(i=i,j=j,tind=tind[tind!=k],cmat=cmat,new.m-1)
				}
			}
		}
		else{
			for(k in tind){
				if(k!=i & k!=j){
					tsum<-tsum+cmat[j,k]*cmat[i,j]
				}
			}
		}
	}
	fbit+1/(tm-1)*tsum
}



### sample communities and calculate steady state species abundances	
sim.comms<-function(N,big.cmat,snum,pvec,max.try=5000){

	namemat<-data.frame(matrix("hello",nrow=snum,ncol=N))
	abmat<-data.frame(matrix(0,nrow=snum,ncol=N))
	n.cycle<-data.frame(matrix(0,nrow=snum,ncol=1))
	for(j in 1:ncol(namemat)){
		namemat[,j]<-as.character(namemat[,j])
	} 
	
	for(i in 1:snum){


		## force unique species in the sample
		go<-TRUE
		count<-0
		while(go){
			tid<-sample(seq(1:nrow(big.cmat)),N,replace=F,prob=pvec)
			go<-I(sum(unlist(as.list(by(substr(names(big.cmat)[tid],1,5),substr(names(big.cmat)[tid],1,5),length)))>1)>0)
	
			count<-count+1
			if(count>max.try){
				go<-FALSE
			}
		}

		##calculate equilibrium abundances
		tout1<-steadyStates(new("markovchain", states=rownames(big.cmat[tid,tid]), transitionMatrix=t(get.pmat(big.cmat[tid,tid]))))

		namemat[i,]<-as.character(colnames(tout1))
		abmat[i,]<-colSums(tout1)/sum(colSums(tout1))

		print(c(N,i))
	}

	#get rid of duplicate entries
	row.names(namemat)<-row.names(abmat)<-seq(1:nrow(namemat))
	namemat<-unique(namemat)
	abmat<-abmat[row.names(namemat),]

	list(namemat=namemat,abmat=abmat)
}


## the main function for getting intran
calc.val<-function(tmat,allp){

	tmat.keep<-tmat
	n.keep<-0

	# get the number of nonzero in the upper triangle for each permutation
	all.vec<-apply(allp,1,get.intran.val,ptmat=tmat) 
	
	##get the max entry
	max.vec<-allp[all.vec==max(all.vec),]
	if(!is.null(nrow(max.vec))){
		max.vec<-max.vec[1,]
	}
	tmat.keep<-tmat[max.vec,max.vec]

	# calculate intransitivity for matrix with most nonzero in upper tri
	retv<-get.int.val(tmat.keep)
	
	retv
}


## calculate number of entries where upper triangle is greater than lower triangle
get.int.val<-function(tvc,ptmat){
	tmat2<-ptmat[tvc,tvc]	
	tmat2<-tmat2-t(tmat2)
	sum(tmat2[upper.tri(tmat2)]>0)
}


## calculates the intran value for a specific matrix ordering, as in Ulrich et al 2014, eqn 18.	
get.intran.val<-function(tmat){

	tn<-nrow(tmat)
	tmat<-ptmat[tvc,tvc]	
	tmat<-tmat-t(tmat)
	tsum<-sum(tmat[upper.tri(tmat)]>0)
	1-(2*tsum/(tn*(tn-1))) 
	### zero means complete hierarchical, 1 means totally intransitive 
}

## the wrapper function that first subsets the data based on the apply call
calc.intran<-function(namevec,big.comp,allp=NA){

	tmat<-as.matrix(big.comp[as.character(namevec),as.character(namevec)])

	cval<-calc.val(tmat,allp=allp)
	
}


