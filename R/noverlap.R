#  noverlap.r
#  ----------
#  Peter J. Rousseeuw, Andreas Christmann
#  Version         MAY/2004
#  Microsoft Windows XP
#  !!! NO WARRANTY !!!
#  This program computes an approximation of n_overlap
#  of a fit T in an (NP+1)-dimensional data set of size N.
#  (actually: the number of points which can be removed such that
#  the reduced data set has no longer overlap, i.e. that the MLE does not exist)
#  Missing values are not allowed.
#  The data set Z has to be a matrix with
#  ncol(Z) -1 < nrow(Z) <= 10000 
#  The first ncol(Z)-1 columns of Z are the design matrix X.
#  The last column of Z is the binary response vector y (0/1).
#  Missing values are not allowed.
#  The paths in the "dos" line and in the "dyn.load" line 
#  have to be specified appropriately.
###########################################################################

noverlap1.for <- function(Z,NDIR=10000) {
  N <- nrow(Z); NP <- ncol(Z)-1; NP1 <- NP+1;
  MA <- matrix(rep(0,N),ncol=1);           RESID <- matrix(rep(0,N),ncol=1);
  JRES <- matrix(rep(0,N),ncol=1);         JSAMP <- matrix(rep(0,NP),ncol=1);
  TVEC <- matrix(c(rep(0,NP),0.5),ncol=1); R <- matrix(rep(0,NP),ncol=1); 
  XN <- matrix(rep(0,N),ncol=1);           EVECS <- matrix(rep(0,NP*NP),ncol=NP);
  EVALS <- matrix(rep(0,NP),ncol=1);       AVE <- matrix(rep(0,NP),ncol=1);  
  COV <- matrix(rep(0,NP*NP),ncol=NP);     RLOC <- matrix(rep(0,NP),ncol=1);
  RSCA <- matrix(rep(0,NP),ncol=1);        XA <- matrix(rep(0,N*NP1),ncol=NP1);
  COEFFS <- matrix(rep(0,NP),ncol=1);      NSIN <- 0;
  NQUASI <- N
  JLV <- matrix(rep(0,N),ncol=1)
  JRV <- matrix(rep(0,N),ncol=1)
  ETAS <- matrix(rep(0,N),ncol=1) 
  Y <- matrix(rep(0,N),ncol=1)
  M <- matrix(rep(0,N),ncol=1)
  YS <- matrix(rep(0,N),ncol=1)
  MS <- matrix(rep(0,N),ncol=1)
  Index <- matrix(rep(0,N),ncol=1)
  FF <- matrix(rep(0,N),ncol=1)
  S <- matrix(rep(0,N),ncol=1)

  .Fortran("novrlp", as.double(Z), 
    as.integer(N), as.integer(NP), as.integer(NP1), as.integer(NDIR), 
    as.integer(MA), as.integer(RESID), as.integer(JRES), as.integer(JSAMP),
    as.double(TVEC), as.double(R), as.double(XN), as.double(EVECS),
    as.double(EVALS), as.double(AVE), as.double(COV), as.double(RLOC), 
    as.double(RSCA), as.double(XA), as.double(COEFFS), as.integer(NSIN), 
    as.integer(NQUASI),
    as.integer(JLV), as.integer(JRV), as.double(ETAS),
    as.integer(Y), as.integer(M), as.integer(YS), as.integer(MS), 
    as.integer(Index), as.integer(FF), as.integer(S),
    PACKAGE="noverlap"
   )} 

noverlap.for <- function(Z,NDIR=10000,PLOT=FALSE) {
    if (nrow(Z) > 10000) {
    stop(message="ERROR: Number of cases > 10,000")}
    if ( nrow(Z) < ncol(Z)-1) {
    stop(message="ERROR: Number of cases < Number of explanatory variables")}
    tmptmp1 <- noverlap1.for(Z,NDIR)
    NP <- ncol(Z)-1
    eta <- Z[,1:NP] %*% as.matrix(tmptmp1[[20]])
    tmptmp2 <- cbind(c(1:length(eta)),Z[,NP+1],eta)
    tmptmp3 <- tmptmp2[order(tmptmp2[,3]),]
    dimnames(tmptmp3) <- list(NULL, c("id", "y", "xu'"))
    if (PLOT == TRUE) {
      plot(tmptmp3[,3],tmptmp3[,2], xlab="xu'", ylab="y", pch="|", cex=1.3,
           lab=c(5,1,7),las=1,ylim=c(0,1),
           main=paste("Plot of response y versus xu'.  n_overlap=",tmptmp1[[22]]))} 

    dimnames(tmptmp3) <- list(NULL, c("id", "y", "xu'"))
    print(tmptmp1[[22]])
    return(tmptmp2 <- list(NOVERLAP=tmptmp1[[22]],COEFFICIENTS=tmptmp1[[20]],
                     NSIN=tmptmp1[[21]],DETAILS=tmptmp3))}


### EXAMPLES
# "Z2" <- matrix(c(-1.5, -1, 0, 0, 1, 1, 2, 3, 3, 3.5, 
#                   0, 3, 1, 2, 2, 4, 2, 1, 3, 4,
#                   0, 1, 0, 0, 0, 0, 1, 1, 1, 1), ncol=3)
#data(Z2)
#noverlap.for(Z2)
#noverlap.for(Z2,NDIR=100000)
#
#par(mfrow=c(2,1))
#noverlap.for(Z2,NDIR=10000,PLOT=TRUE)
#tmp <- noverlap.for(Z2)
#tmp$NOVERLAP
#tmp$COEFFICIENTS
#tmp$NSIN
#tmp$DETAILS
#Z3 <- as.data.frame(Z2)
#names(Z3) <- c("x1","x2","y")
#plot(x2 ~x1, data=Z3,pch=as.character(y),main="Scatterplot")
#abline(c(0,1.5),col="blue")
#
