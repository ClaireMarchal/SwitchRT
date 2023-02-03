require("foreign",quietly=TRUE)
require("fields",quietly=TRUE)
require("qvalue",quietly=TRUE)


##################################################
############### Input common fuctions  ###########
##################################################

getDistances = function(RT, set1, set2) {  														
  R1 = RT[,c(set1, set2)]
  R2 = RT[,c(set1, set2)]
  R1 = as.data.frame(t(na.omit(R1)))
  R2 = as.data.frame(t(na.omit(R2)))
  R1test = as.matrix(R1)
  R2test = as.matrix(R2)
  datNum  = nrow(R1)
  compMat = matrix(1:(datNum*datNum), nrow=datNum, ncol=datNum)		
  set1Ind = 1:length(set1)
  set2Ind = (max(set1Ind)+1) : (length(set2)+length(set1))
  LCTa = compMat[set1Ind,set1Ind]						# LCT, fingerprint set
  LCTb = compMat[set2Ind,set2Ind]						# LCT, non-fingerprint set
  DCT = c(compMat[set1Ind,set2Ind], compMat[set2Ind,set1Ind])		# DCT (combined)
  LCT = sort(c(LCTa, LCTb))						# LCT, fingerprint vs. fingerprint, non vs. non
  DCT = sort(DCT)								# DCT, fingerprint vs. non-fingerprint comparisons
  require(foreign)							# For rdist matrix euclidean distance
  require(fields)	    
  
  DRs = NULL							
  reg = dim(R1)[2]
  
  # First pass for distances								
  for(p in 1:reg) {	
    Distances = rdist(R1test[,p],R2test[,p])			# Distances among sets
    within  = Distances[c(LCT)]					# Distances within set1, set2
    within  = within[within > 0.0001]				# Remove diagonal 0s
    between = Distances[c(DCT)]					# Distances between set1, set2
    DRs$local.p.value[p] = t.test(between, within, alternative="greater")$p.value	# P-value for current set agains individual probe/region/gene
    DRs$within[p]	= mean(within)					# Average distances within set1, set2 
    DRs$between[p]	= mean(between)					# Average distances between set1, set2
    DRs$DR[p]		= mean(between)/mean(within)		# Distance ratio
    DRs$Dom[p]		= p		
  }
  return(DRs)
}

getPvalues = function(DRs) {										#Second pass for global tests of between vs. within group distances 
  withinDist = ecdf(DRs$within)
  for(p in 1:length(DRs$Dom)) {	
    DRs$global.p.value[p]= 1-withinDist(DRs$between[p])			# P-value for current set against global distribution of within set distances
  }
  return(DRs)
}

getQvalues = function(DRs) {										
  library(qvalue)
  qobj <- qvalue(DRs$global.p.value, gui=F)
  #hist(qobj$qvalues)
  #qplot(qobj)
  #qsummary(qobj)
  DRs$q.value = qobj$qvalues
  return(DRs)
}

