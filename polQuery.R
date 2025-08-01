Arguments <- commandArgs(trailingOnly = T)
#Arguments <- c('-D2025/01/10/04:28:20', '-F343.500000', 'J0922-3959')
timeWindow <- 60	# Days
#-------- Function to return residuals in RM fit
residEVPA <- function(x, y, w){	# Optimization function for RM fit
	return(function(para){
		RM <- para[1]
		EVPA <- para[2]
		return( sum(w* sin(y - RM*x - EVPA)^2 ))
	})
}
#-------- Parse arguments
parseArg <- function( args ){
	srcNum <- argNum <- length(args)
	for( index in 1:argNum ){
		if(substr(args[index], 1,2) == "-D"){ refDate <- as.Date(substring(args[index], 3));    srcNum <- srcNum - 1 }
		if(substr(args[index], 1,2) == "-F"){ refFreq <- as.numeric(substring(args[index], 3)); srcNum <- srcNum - 1}
	}
	srcList = args[(argNum - srcNum + 1):argNum]
	return(list(refDate = refDate, refFreq = refFreq, srcList = srcList[grep('^J[0-9]',srcList )]))
}

argList <- parseArg(Arguments)
refDate <- argList$refDate
refFreq <- argList$refFreq
srcList <- argList$srcList

#-------- Load Flux.Rdata from web
#if( file.exists('Flux.Rdata') == FALSE ){ download.file("https://www.alma.cl/~skameno/AMAPOLA/Flux.Rdata", "Flux.Rdata") }
download.file("https://www.alma.cl/~skameno/AMAPOLA/Flux.Rdata", "Flux.Rdata")
load('Flux.Rdata')
#load(url("https://www.alma.cl/~skameno/Grid/Stokes/Flux.Rdata")) 
FLDF$timeDiff <- as.numeric(difftime(FLDF$Date, refDate, units='days'))
#-------- For each source
IatRef <- QatRef <- UatRef <- numeric(0)
for(sourceName in srcList){
	srcDF <- FLDF[((FLDF$Src == sourceName) & (abs(FLDF$timeDiff) < timeWindow)),]
	#if(nrow(srcDF) < 6){ srcList <- srcList[-which(srcList %in% sourceName)]; next }
	if(nrow(srcDF) < 4){ srcList <- srcList[-which(srcList %in% sourceName)]; next }
	if(min(abs(srcDF$timeDiff)) > timeWindow){ srcList <- srcList[-which(srcList %in% sourceName)]; next }
    #for( freq in freqList ){
    #    if( nrow(srcDF[srcDF$Freq == freq,]) < 2 ){ srcDF <- srcDF[srcDF$Freq != freq,]}
    #}
	freqList <- as.numeric(unique(srcDF$Freq))
    freqNum <- length(freqList)
	if(freqNum < 2){ srcList <- srcList[-which(srcList %in% sourceName)]; next }
	if(diff(range(freqList)) < 50.0){ srcList <- srcList[-which(srcList %in% sourceName)]; next }  # frequenc range should be > 50 GHz
	freqList <- freqList[order(freqList)]
	estI <- errI <- estQ <- errQ <- estU <- errU <- numeric(freqNum)
	#-------- For each frequency
    if(freqNum > 1){
	    for(freq_index in 1:freqNum){
	    	srcFreqDF <- srcDF[srcDF$Freq == freqList[freq_index],]
	    	if(((nrow(srcFreqDF) >= 3) & (diff(range(srcFreqDF$timeDiff)) > min(srcFreqDF$timeDiff)))){
			    fit <- lm(data=srcFreqDF, formula=I ~ timeDiff, weights=(I/(eI + 1.0e-3)) / abs(timeDiff + 5))
			    estI[freq_index] <- summary(fit)$coefficients[1,'Estimate'];  errI[freq_index] <- summary(fit)$coefficients[1,'Std. Error']
			    fit <- lm(data=srcFreqDF, formula=Q ~ timeDiff, weights=(I/(eQ + 1.0e-3)) / abs(timeDiff + 5))
			    estQ[freq_index] <- summary(fit)$coefficients[1,'Estimate'];  errQ[freq_index] <- summary(fit)$coefficients[1,'Std. Error']
			    fit <- lm(data=srcFreqDF, formula=U ~ timeDiff, weights=(I/(eU + 1.0e-3)) / abs(timeDiff + 5))
			    estU[freq_index] <- summary(fit)$coefficients[1,'Estimate'];  errU[freq_index] <- summary(fit)$coefficients[1,'Std. Error']
		    } else {
			    estI[freq_index] <- median(srcFreqDF$I); errI[freq_index] <- median(srcFreqDF$eI) * 10.0
			    estQ[freq_index] <- median(srcFreqDF$Q); errQ[freq_index] <- median(srcFreqDF$eQ) * 10.0
			    estU[freq_index] <- median(srcFreqDF$U); errU[freq_index] <- median(srcFreqDF$eU) * 10.0
            }
            if(estI[freq_index] < 0.0){ estI[freq_index] <- 1.0e-6; estQ[freq_index] <- estU[freq_index] <- 1.0e-8 }
	    }
	    lambdaSQ <- (0.299792458 / freqList)^2; lambdasqSpan <- diff(range(lambdaSQ))
	    estP <- sqrt(estQ^2 + estU^2); errP <- 0.01*estI + sqrt(errQ^2 + errU^2); estEVPA <- 0.5*atan2(estU, estQ)
	    fit <- lm(log(estI) ~ log(freqList/100.0), weights=1.0/errI^2); I100 <- exp(as.numeric(coef(fit)[1])); spixI <- as.numeric(coef(fit)[2])
	    fit <- lm(log(estP + 1.0e-6) ~ log(freqList/100.0), weights=1.0/errP^2); P100 <- exp(as.numeric(coef(fit)[1])); spixP <- as.numeric(coef(fit)[2])
	    estEVPAend <- estEVPA[freqNum]
	    if(estEVPAend - estEVPA[1] >  pi/2){ estEVPAend <- estEVPAend - pi }
	    if(estEVPAend - estEVPA[1] < -pi/2){ estEVPAend <- estEVPAend + pi }
	    RMinit <- (estEVPAend - estEVPA[1]) / (lambdaSQ[freqNum] - lambdaSQ[1])
	    fit <- optim(par=c(RMinit, estEVPA[freqNum]), fn=residEVPA(lambdaSQ, estEVPA, 1.0/errP^2), method='Nelder-Mead')
	    RM <- fit$par[1]; EVPAintercept <- fit$par[2]
	    PatFreqRef <- P100*(refFreq/100)^spixP
	    EVPAatFreqRef <- RM* (0.299792458 / refFreq)^2 + EVPAintercept
	    IatRef <- append(IatRef, I100*(refFreq/100)^spixI)
	    QatRef <- append(QatRef, PatFreqRef * cos(2.0* EVPAatFreqRef))
	    UatRef <- append(UatRef, PatFreqRef * sin(2.0* EVPAatFreqRef))
    } else {
        spix <- -0.7
        IatRef <- append(IatRef, sum(srcDF$I/srcDF$eI^2)/sum(1.0/srcDF$eI^2)* (refFreq/median(srcDF$Freq))^spix)
        QatRef <- append(QatRef, sum(srcDF$Q/srcDF$eQ^2)/sum(1.0/srcDF$eQ^2)* (refFreq/median(srcDF$Freq))^spix)
        UatRef <- append(UatRef, sum(srcDF$U/srcDF$eU^2)/sum(1.0/srcDF$eU^2)* (refFreq/median(srcDF$Freq))^spix)
    }
}
DF <- data.frame(Src=srcList, I=IatRef, Q=QatRef, U=UatRef)
write.table(na.omit(DF), file="CalQU.data", append=F, quote=F, col.names=F, row.name=F)
