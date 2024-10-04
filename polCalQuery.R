#-------- Parse arguments
parseArg <- function( args ){
	argList <- list(NA, NA, 3600, 3, 0.05, NA, FALSE)
	names(argList) <- c('startTime', 'RA', 'execDuration', 'Band', 'threshFlux', 'refFreq', 'load')
	if( !file.exists("Flux.Rdata") ){ argList$load <- TRUE }
	argNum <- length(args)
	for( index in 1:argNum ){
		if(substr(args[index], 1,2) == "-s"){ argList$startTime <- as.character(substring(args[index], 3)) }
		if(substr(args[index], 1,2) == "-r"){ argList$RA <- as.numeric(substring(args[index], 3)) }
		if(substr(args[index], 1,2) == "-d"){ argList$execDuration <- as.numeric(substring(args[index], 3)) }
		if(substr(args[index], 1,2) == "-b"){ argList$Band <- as.integer(substring(args[index], 3)) }
		if(substr(args[index], 1,2) == "-t"){ argList$threshFlux <- as.numeric(substring(args[index], 3))}
		if(substr(args[index], 1,2) == "-#"){ argList$refFreq <- as.numeric(substring(args[index], 3))}
		if(substr(args[index], 1,2) == "-L"){ argList$load <- TRUE}
	}
    #-------- Automated start time setting
    if( is.na(argList$startTime) ){ argList$startTime <- as.character(as.POSIXct(Sys.time())) }  # Set NOW
    if( is.na(argList$RA) == FALSE ){
        startLST <- argList$RA - argList$execDuration / 3600    # Straddling transit
        NowLST <- 12* (mjd2gmst(ISO8601mjdSec(argList$startTime) / SEC_PER_DAY) + ALMA_LONG)/pi # LST in [hour]
        startSEC <- (ISO8601mjdSec(argList$startTime) + 3600* (startLST - NowLST)) %% SEC_PER_DAY
        startTime <- sprintf('%sT%02d:%02d:%02d', strsplit(argList$startTime, '[T| ]')[[1]][1], as.integer(startSEC %/% 3600), as.integer(startSEC %/% 60 %% 60), as.integer(startSEC %% 60))
        argList$startTime <- startTime
        cat(sprintf('Because RA is specified, start time is set to %s straddling transit.\n', argList$startTime))
    }
	return(argList)
}
#-------- Time Constants
SEC_PER_DAY <- 86400
MJD_1901 <- 15384
DAY_PER_YEAR <- 365
DAY_PER_4YEAR <- 1461
DOY_MON <- c(0, 31, 59, 90, 120, 151, 181, 212, 242, 273, 303, 334)
#-------- FE-specific PA
#         Band1      2     3      4     5      6     7      8      9   10
BandPA <- c(45.0, -45.0, 80.0, -80.0, 45.0, -45.0, 36.45, 90.0, -90.0, 0.0)*pi/180
BandFreq <- c(43.0, 75.0, 97.5, 132.0, 183.0, 233.0, 343.5, 400.0, 650.0, 800.0)
ALMA_LAT <- -23.029* pi/180.0   # [radian]
ALMA_LONG <- -67.755* pi/180.0  # [radian]
cos_phi <- cos(ALMA_LAT)
sin_phi <- sin(ALMA_LAT)
maxSinEL <- sin(86/180*pi)
minSinEL <- sin(20/180*pi)
#-------- Calculate Day of year from Month and date
md2doy <- function(year, month, date){
	is_leap <- ((year%%4 == 0) && (month > 2))	# Leap Year Frag
	DOY_MON[month] + date + is_leap
}
#-------- Calculate (fractional) Modified Julian Date in unit of second. This unit is used in CASA
doy2mjdSec <- function(year, doy, hour, minute, second){
	if((year < 1901) || (year > 2099)){ return(-1)}
	year <- year - 1901
	mjd <- year %/%4 * DAY_PER_4YEAR + year%%4 * DAY_PER_YEAR + doy + MJD_1901
	sod <- (hour*60 + minute)*60 + second
	return(mjd* SEC_PER_DAY + sod)
}
#-------- Convert ISO8601 format string into mjdSec
ISO8601mjdSec <- function( timeString ){				# Input string YYYY-MM-DDTHH:MM:SS.S
	year <- as.integer(substring(timeString, 1, 4))
	month <- as.integer(substring(timeString, 6, 7))
	day <- as.integer(substring(timeString, 9, 10))
	UT_hour <- as.integer(substring(timeString, 12, 13))
	minute <- as.integer(substring(timeString, 15, 16))
	second <- as.integer(substring(timeString, 18, 19))
	return(doy2mjdSec(year, md2doy(year, month, day), UT_hour, minute, second))
}
#-------- MJD to Greenwich mean sidereal time [radian]
mjd2gmst <- function(mjd, ut1utc = 0){
	# mjd2gmst : Calculate Greenwidge Mean Sidereal Time
	# mjd : Modified Julian Date
	# ut1utc : UT1 - UTC [sec]
	FACT <- c(24110.54841, 8640184.812866, 0.093104, 0.0000062)
	MJD_EPOCH <- 51544.5  # MJD at 2000 1/1 12:00:00 UT
	TU_UNIT <- 36525.0
	SEC_PER_DAY   <- 86400
	
	tu <- (mjd - MJD_EPOCH) / TU_UNIT
	ut1 <- (mjd%%1)* SEC_PER_DAY + ut1utc
	gmst <- (ut1 + FACT[1] + ((FACT[4]* tu + FACT[3])* tu + FACT[2])* tu) / SEC_PER_DAY
	return(2* pi* (gmst %% 1))
}
#-------- # cos(hour angle) when it passes the given EL
EL_HA <- function(sinEL, dec){
	cosHA <- (sinEL - sin_phi* sin(dec)) / (cos_phi* cos(dec))
	return(cosHA)
}
#-------- # filtering by source polarization
sourceDataFrame <- function(DF, refFreq=100.0, refDate=Sys.Date()){
    DateRange <- 60    # 60-day window
    #-------- Filter by observing date
	DF <- DF[abs(as.Date(DF$Date) - as.Date(refDate)) < DateRange,]
    sourceList <- unique(DF$Src)
    sourceList <- sourceList[grep('^J[0-9]', sourceList)]  # Filter SSOs out
    SDF <- data.frame( matrix(rep(NA, 9), nrow=1))[numeric(0),]
    colnames(SDF) <- c('Src', 'RA', 'DEC', 'I', 'Q', 'U', 'V', 'P', 'EVPA')
	for(src in sourceList){
		srcDF <- DF[((DF$Src == src) & (abs(as.Date(DF$Date) - as.Date(refDate)) < DateRange)),] 
        srcDF$P  <- sqrt(srcDF$Q^2 + srcDF$U^2)
        srcDF$eP <- sqrt(srcDF$eQ^2 + srcDF$eU^2)
        srcDF$EVPA  <- 0.5* atan2(srcDF$U, srcDF$Q)
        srcDF$eEVPA <- 0.5* sqrt(srcDF$Q^2 * srcDF$eU^2 + srcDF$U^2 * srcDF$eQ^2) / (srcDF$P)^2
        srcDF$relTime <- as.numeric(srcDF$Date) - as.numeric(as.POSIXct(refDate))
        if( (diff(range(srcDF$Freq)) > 100.0) & ( min(abs(srcDF$relTime)) / diff(range(srcDF$relTime)) < 2 ) ){
            tempDF <- estimateIQUV(srcDF, refFreq)
            tempDF$RA  <- pi* (60.0* as.numeric(substring(src, 2, 3)) + as.numeric(substring(src, 4, 5))) / 720.0
            tempDF$DEC <- pi* sign(as.numeric(substring(src, 6, 10)))* (as.numeric(substring(src, 7, 8)) + as.numeric(substring(src, 9, 10))/60.0) / 180.0
            SDF <- rbind(SDF, tempDF)
        }
    }
    return( SDF )
}
#-------- Estimate Stokes parameters by freqneyc and date 
estimateIQUV <- function(DF, refFreq){
    DF$relFreq <- DF$Freq / refFreq
    if(length(unique(DF$relFreq)) < 4 ){ return( data.frame( Src=DF$Src[1], I=0.0, Q=0.0, U=0.0, V=0.0, P=0.0, EVPA=0.0))}
    if(length(unique(DF$relTime)) < 4 ){
        fitI <- lm(formula=log(I) ~ log(relFreq), data=DF, weight=(I / eI)^2 * (864000 / abs(relTime + 864000)) )
	    fitP <- lm(formula=log(P) ~ log(relFreq), data=DF, weight=(DF$P/DF$eP)^2 * (864000 / abs(relTime + 864000)))
    } else {
        fitI <- lm(formula=log(I) ~ relTime + log(relFreq), data=DF, weight=(I / eI)^2 * (864000 / abs(relTime + 864000)) )
	    fitP <- lm(formula=log(P) ~ relTime + log(relFreq), data=DF, weight=(P/eP)^2 * (864000 / abs(relTime + 864000)))
    }
    weight <- 1.0/(abs(DF$eEVPA)^2 * abs(log(DF$relFreq) + 1.0)^2 * (864000 / abs(DF$relTime + 864000)))
    Twiddle <- sum( weight* exp((0.0 + 2.0i)*DF$EVPA) ) / sum(weight)
    IQUV <- data.frame(Src=DF$Src[1], I=exp(coef(fitI)[[1]]), Q=0.0, U=0.0, V=0.0, P=exp(coef(fitP)[[1]]), EVPA=0.5*Arg(Twiddle))
    IQUV$Q <- IQUV$P* Re(Twiddle)
    IQUV$U <- IQUV$P* Im(Twiddle)
	return( IQUV )
}
#-------- Input parameters
#Arguments <- strsplit("-s2024-10-01T03:24:00 -r12.5 -d7200 -b3", ' ')[[1]]     # for debugging
#Arguments <- strsplit("-d7200 -b3", ' ')[[1]]     # for debugging
Arguments <- commandArgs(trailingOnly = TRUE)
argList <- parseArg(Arguments)
setwd('./')
if(is.na(argList$refFreq)){ argList$refFreq <- BandFreq[argList$Band] }
startmjdSec <- ISO8601mjdSec(argList$startTime)
endmjdSec <- startmjdSec + argList$execDuration
startLST <- mjd2gmst(startmjdSec/SEC_PER_DAY) + ALMA_LONG
endLST   <- mjd2gmst(endmjdSec/SEC_PER_DAY) + ALMA_LONG
if( endLST < startLST){ endLST <- endLST + 2*pi}
#-------- Load Flux.Rdata from web
if( argList$load ){
	cat('--- Loading Flux.Rdata from the web\n')
	FluxDataURL <- "https://www.alma.cl/~skameno/Grid/Stokes/"
	load(url(paste(FluxDataURL, "Flux.Rdata", sep='')))     # Data frame of FLDF
    save(FLDF, file='Flux.Rdata')
} else {
	load('Flux.Rdata')
}
#-------- Filter out single-sideband Band-3 frequency 
pos <- regexpr("RB",FLDF$File)
FLDF$Band <- as.integer(substr(FLDF$File, pos+3, pos+4))
FLDF$BandPA <- BandPA[FLDF$Band]
FLDF <- FLDF[-which((FLDF$Band == 3) & (abs(FLDF$Freq - 97.45) > 1.0)),]
SDF <- sourceDataFrame( FLDF, argList$refFreq, argList$startTime)
#---- Filter by P > 0.05 Jy
SDF <- SDF[SDF$P > argList$threshFlux,]
#---- Filter by m > 3%
SDF <- SDF[SDF$P / SDF$I > 0.03,]
#---- Filter by EL
if( is.na(argList$RA) ){
    SDF$startHA <- startLST - SDF$RA
    SDF$endHA   <- endLST - SDF$RA
} else {
    SDF$startHA <- pi* (argList$RA/12.0 - argList$execDuration/SEC_PER_DAY) - SDF$RA
    SDF$endHA   <- pi* (argList$RA/12.0 + argList$execDuration/SEC_PER_DAY) - SDF$RA
}
#cat(sprintf('%s: HA = %.1f - %.1f\n', SDF$Src, SDF$startHA, SDF$endHA))
SDF <- SDF[ which(cos(SDF$startHA) > EL_HA(minSinEL, SDF$DEC)), ]	# EL >20º at the beggining
SDF <- SDF[ which(cos(SDF$endHA)   > EL_HA(minSinEL, SDF$DEC)), ] # EL >20º at the end
SDF <- SDF[ which(cos(SDF$startHA) < EL_HA(maxSinEL, SDF$DEC)), ] # EL < 86º at the beggining
SDF <- SDF[ which(cos(SDF$endHA)   < EL_HA(maxSinEL, SDF$DEC)), ] # EL < 86º at the end
SDF <- SDF[ which( (EL_HA(maxSinEL, SDF$DEC) > 1.0) | (sin(SDF$startHA)* sin(SDF$endHA) > 0)),] # transit for EL>86º
#-------- For each source, calculate polarizaiton responses
sourceList <- unique(SDF$Src)
SDF$Score <- SDF$QC_H <- SDF$QC_L <- SDF$QC_0 <- SDF$UC_H <- SDF$UC_L <- SDF$UC_0 <- numeric(length(sourceList))
for(src in sourceList){
	index <- which(SDF$Src == src)
	HA <- seq(SDF$startHA[index], SDF$endHA[index], by=0.01)
	cos_HA  <- cos(HA)
	sin_HA  <- sin(HA)
	cos_dec <- cos(SDF$DEC[index])
	sin_dec <- sin(SDF$DEC[index])
	PA <- atan2(sin_HA, sin_phi*cos_dec/cos_phi - sin_dec* cos_HA) + BandPA[argList$Band]
	QCpUS <- SDF$Q[index] * cos(2.0* PA) + SDF$U[index] * sin(2.0* PA)
	UCmQS <- SDF$U[index] * cos(2.0* PA) - SDF$Q[index] * sin(2.0* PA)
	SDF$UC_0[index] <- (HA[which.min( abs(UCmQS) )] + SDF$RA[index])
	SDF$QC_0[index] <- (HA[which.min( abs(QCpUS) )] + SDF$RA[index])
	SDF$QC_H[index] <- max(QCpUS)
	SDF$QC_L[index] <- min(QCpUS)
	SDF$UC_H[index] <- max(UCmQS)
	SDF$UC_L[index] <- min(UCmQS)
}
SDF$Score <- (0.55 -0.45* sign(SDF$UC_H * SDF$UC_L))* (SDF$UC_H - SDF$UC_L)         # Scoring by range of liner polaration responses
SDF <- SDF[order(SDF$Score, decreasing=TRUE),]                                      # Sort by score
#-------- Output results
text_sd <- 'Source      Score   I      P      EVPA  UC_L   UC_H   UC_0   QC_0\n'
cat(text_sd)
cat(text_sd, file='CalQU.data')
for(index in 1:nrow(SDF)){
    text_sd <- sprintf('%s %6.3f %6.3f %6.3f %6.1f %6.3f %6.3f %6.2f %6.2f\n', SDF[index,]$Src, SDF[index,]$Score, SDF[index,]$I, SDF[index,]$P, 180*SDF[index,]$EVPA/pi, SDF[index,]$UC_L, SDF[index,]$UC_H, SDF[index,]$UC_0*12/pi, SDF[index,]$QC_0*12/pi)
    cat(text_sd)
    cat(text_sd, file='CalQU.data', append=T)
}
