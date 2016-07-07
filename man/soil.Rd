\name{SOIL}
\alias{SOIL}
\title{Sparsity Oriented Importance Learning (SOIL)}
\description{
The package calculates the Sparsity Oriented Importance Learning (SOIL) measure, which provides 
an objective and informative profile of variable importances for high dimensional 
regression and classification models.}
\usage{
SOIL(x, y, n_train = ceiling(n/2), no_rep = 100,
                 n_train_bound = n_train - 2, n_bound = n - 2,
                 psi = 1, family = c("gaussian",
                 "binomial"), method = c("union", "customize"),
                 candidate_models, weight_type = c("BIC", "AIC",
                 "ARM"), prior = TRUE, reduce_bias = FALSE)
}
\arguments{
  \item{x}{Matrix of predictors.}

  \item{y}{Response variable.}

  \item{n_train}{Size of training set when the weight function is \code{ARM} or \code{ARM} with \code{prior=TRUE}. The default value is \code{n_train=ceiling(n/2).}}

  \item{no_rep}{Number of replications when the weight function is \code{ARM} and \code{ARM} with \code{prior=TRUE}. The default value is \code{no_rep=100}.}

  \item{n_train_bound}{When computing the weights using \code{ARM}, the candidate models with the size larger than \code{n_train_bound} will be dropped. The default value is \code{n_train-2}.}
  \item{n_bound}{When computing the weights using \code{AIC} or \code{BIC}, the candidate models with the size larger than \code{n_train_bound} will be dropped. The default value is \code{n-2}.}
  \item{psi}{A positive number to control the improvement of the prior weight. The default value is 1.}
  \item{family}{Choose the family for GLM models. So far \code{gaussian} and \code{binomial} are implemented. The default is \code{gaussian}.}
  \item{method}{Users can choose either \code{union} or \code{customize}. If \code{method=="union"}, then the program automatically provides the candidate models as a union of solution paths of Lasso, SCAD, and MCP; If \code{method="customize"}, users must provide their own set of candidate models in the input argument \code{candidate_models} as a matrix, each row of which is a 0/1 index vector representing whether each variable is included/excluded in the model. For details see Example section.}
  \item{candidate_models}{Only available when \code{method="customize"}. It is a matrix of candidate models, each row of which is a 0/1 index vector representing whether each variable is included/excluded in the model. For details see Example section.}
  \item{weight_type}{Options for computing weights for SOIL measure. Users can choose among \code{ARM}, \code{AIC} and \code{BIC}. The default is \code{BIC}.}
  \item{prior}{Whether to use prior in the weighting function. The default is \code{TRUE}.}
  \item{reduce_bias}{If the binomial model is used, occasionally the algorithm might has convergence issue when the problem of so-called complete separation or quasi-complete separation happens. Users can set \code{reduce_bias=TRUE} to solve the issue. The algorithm will use an adjusted-score approach when ftting the binomial model for computing the weights. This method is developed in Firth, D. (1993). Bias reduction of maximum likelihood estimates. Biometrika 80, 27-38.}
}

\details{
See the paper provided in Reference section. 
}

\value{
A "SOIL" object is retured. The components are:
\item{importance}{SOIL importance values for each variable.}
\item{weight}{The weight for each candidate model.}
\item{candidate_models_cleaned}{Cleaned candidate models: the duplicated candidate models are cleaned; When computing SOIL weights using AIC and BIC, the models with more than n-2 variables are removed (n is the number of observaitons); When computing SOIL weights using ARM, the models with more than n_train-2 variables  are removed (n_train is the number of training observations).}
}

\references{
Ye, Chenglong, Yang, Y., Yang, Y. (2016), "Sparsity Oriented Importance Learning," \emph{preprint}.\cr
BugReport: \url{https://github.com/emeryyi/SOIL}\cr
}

\examples{
# REGRESSION CASE

# generate simulation data
n <- 50
p <- 8
beta <- c(3,1.5,0,0,2,0,0,0)
sigma <- matrix(0,p,p)
for(i in 1:p){
   for(j in 1:p) sigma[i,j] <- 0.5^abs(i-j)
}
x <- mvrnorm(n, rep(0,p), sigma)
e <- rnorm(n)
y <- x \%*\% beta + e


# compute SOIL using ARM with prior
v_ARM <- SOIL(x, y, family = "gaussian",
weight_type = "ARM", prior = TRUE)

# compute SOIL using AIC
v_AIC <- SOIL(x, y, family = "gaussian", 
weight_type = "AIC", prior = TRUE)

# compute SOIL using BIC
v_BIC <- SOIL(x, y, family = "gaussian", weight_type = "BIC")

# user supplied candidate models
candidate_models = rbind(c(0,0,0,0,0,0,0,1), 
c(0,1,0,0,0,0,0,1), c(0,1,1,1,0,0,0,1), 
c(0,1,1,0,0,0,0,1), c(1,1,0,1,1,0,0,0), 
c(1,1,0,0,1,0,0,0))

v1_BIC <- SOIL(x, y, 
psi=1, 
family = "gaussian",
method = "customize", 
candidate_models = candidate_models, 
weight_type = "BIC", prior = TRUE)

# CLASSIFICATION CASE

# generate simulation data
n = 300 
p = 8
b <- c(1,1,1,-3*sqrt(2)/2) 
x=matrix(rnorm(n*p, mean=0, sd=1), n, p) 
feta=x[, 1:4]\%*\%b 
fprob=exp(feta)/(1+exp(feta))
y=rbinom(n, 1, fprob)


# compute SOIL for model_check using BIC with prior
b_BIC <- SOIL(x, y, family = "binomial", weight_type = "BIC")

candidate_models = 
rbind(c(0,0,0,0,0,0,0,1), 
c(0,1,0,0,0,0,0,1), 
c(1,1,1,1,0,0,0,0), 
c(0,1,1,0,0,0,0,1), 
c(1,1,0,1,1,0,0,0), 
c(1,1,0,0,1,0,0,0), 
c(0,0,0,0,0,0,0,0),
c(1,1,1,1,1,0,0,0))

# compute SOIL for model_check using AIC
# user supplied candidate models
b_AIC <- SOIL(x, y, family = "binomial",
method = "customize", candidate_models = candidate_models, 
weight_type = "AIC")
}

