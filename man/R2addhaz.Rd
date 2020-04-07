\name{R2addhaz}
\alias{R2addhaz}
\title{
Estimate R^2 for additive hazards model
}
\description{
The function computes R^2 measure of explained variation under the semiparametric additive hazards model.
}
\usage{
R2addhaz(data)
}
\arguments{
  \item{data}{a data.frame with survival data. The first column needs to be the censored failure time. The second column needs to be the event indicator, 1 if the event is observed, 0 if it is censored. The other columns are covariates.
}
}
\details{
The semiparametric hazards model\deqn{\lambda(t | Z)=\lambda_0(t) + \beta Z}{\lambda(t | Z)=\lambda_0(t) + \beta Z} is fitted to the data. The R^2 measure of explained variation is then computed. }
\value{
\item{R}{R^2 measure of explained variation.}
}
\references{
Rava, D., Xu, R. "Explained Variation under the Additive Hazards Model", March 2020, arXiv:2003.09460}
\author{
Denise Rava
}

\examples{
Z=runif(100,0,sqrt(3)) #generate covariates
u=runif(100,0,1)
t=-log(u)/as.vector((1+Z)) #generate failure time
status=rep(1,100) #censoring indicator
sd<-as.data.frame(cbind(t,status,Z)) #data frame of survival data
R2addhaz(sd)}
\keyword{ ~R2 }
\keyword{ ~AdditiveHazardsModel }
