\name{noverlap.for}
\alias{noverlap.for}
\alias{noverlap}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Minimal number of overlap based on an affine hyperplane in
       binary regression (complete or quasicomplete separation)}
\description{
Applies the regression depth method (RDM) to binary regression.
This method computes approximately the number of data points that can be removed
from a data set such that the remaining data set has NO overlap, i.e.
that the remaining data set has complete separation or quasi-complete separation.
If Noverlap=0 then these maximum likelihood estimates for the parameter vector
in many binary regression models such as logistic regression or probit regression do not exist.
}
\usage{
noverlap.for(Z,NDIR=10000,PLOT=FALSE)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{Z}{The data set Z has to be a matrix with ncol(Z)-1 <= nrow(Z) <= 10000.
           The first ncol(Z)-1 columns of Z are the design matrix X.
           The last column of Z is the binary response vector y (0/1).}
  \item{NDIR}{Maximal number of directions (integer)}
  \item{PLOT}{logical, if TRUE then draw a plot}
}
\value{
A list with components
  \item{Z}{data matrix}
  \item{NDIR}{10000 directions}
  \item{PLOT}{FALSE, i.e. make no plot}
}
\author{Andreas Christmann, Peter J. Rousseeuw
\email{Christmann@statistik.uni-dortmund.de}}
\references{
Christmann, A., Rousseeuw, P.J. (2001).
  Measuring overlap in logistic regression.
  \emph{Computational Statistics and Data Analysis}, \bold{37}, 65-75.

Christmann, A. (2002).
  \emph{Classification based on the support vector machine and on regression depth.}
  In: Y. Dodge (Ed.):
  Statistical Data Analysis Based on the L1-Norm and Related Methods.
  Series: Statistics for industry and technology. Birkhaeuser, Basel, pp. 341-352.

Christmann, A., Fischer, P., Joachims, T. (2002).
    Comparison between various regression depth methods and the
    support vector machine to approximate the minimum
    number of misclassifications.
    \emph{Computational Statistics}, \bold{17}, 273-287.
}
\seealso{\code{\link{noverlap1.for}}}
\keyword{regression}
\keyword{robust}
\examples{
data(Z2)
noverlap.for(Z2)
noverlap.for(Z2,NDIR=100000)
# x11()
postscript(file="tmp1.ps")
par(mfrow=c(2,1))
noverlap.for(Z2,NDIR=10000,PLOT=TRUE)
tmp <- noverlap.for(Z2)
tmp$NOVERLAP
tmp$COEFFICIENTS
tmp$NSIN
tmp$DETAILS
Z3 <- as.data.frame(Z2)
names(Z3) <- c("x1","x2","y")
plot(x2 ~ x1, data=Z3,pch=as.character(y),main="Scatterplot")
abline(c(0,1.5),col="blue")
points(Z3[2,1],Z3[2,2],pch=as.character(Z3[2,3]),col="red")
dev.off()

# NO OVERLAP: maximum likelihood estimates do NOT exist
data(Z1)
Z1
# X11()
postscript(file="tmp2.ps")
noverlap.for(Z1)
tmp <- noverlap.for(Z1)
tmp$NOVERLAP
tmp$COEFFICIENTS
tmp$NSIN
tmp$DETAILS
Z3 <- as.data.frame(Z1)
names(Z3) <- c("x1","y")
plot(y ~ x1, data=Z3,pch=as.character(y),main="Scatterplot")
summary(glm(y ~ x1, data=Z3, family=binomial(link=logit), trace=TRUE, maxit=30))
dev.off()

# NO OVERLAP: maximum likelihood estimates in the logistic regression model
# do NOT exist for the banknotes data set
data(Banknotes)
Banknotes
# X11()
postscript(file="tmp3.ps")
tmp <- noverlap.for(Banknotes,PLOT=TRUE)
dev.off()
tmp$NOVERLAP
tmp$COEFFICIENTS
Z3 <- as.data.frame(Banknotes)
names(Z3) <- c("x1","x2", "x3", "x4", "x5", "x6","y")
summary(glm(y ~ x1+x2+x3+x4+x5+x6, data=Z3, family=binomial(link=logit), trace=TRUE, maxit=30))
}
