\name{mcmc_nhpp_erlmix}
\alias{do.mcmc.nhpp.erlmix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{MCMC on Erlang Mixture Model
%%  ~~function to do ... ~~
}
\description{Run MCMC on an Erlang mixture model for nonhomogeneous Poisson process intensities.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
do.mcmc.nhpp.erlmix(nsim,t,maxT,pri.param,init=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nsim}{a positive scalar, the number of MCMC iterations.}
  \item{t}{a numeric vector, the data used for intensity estimation.}
  \item{maxT}{a positive scalar, the upper limite of an observation window.}
  \item{pri.param}{list of numeric values, prior choices for an Erlang mixture model  (refer to Details for the full list of (hyper)paramters).}
  \item{init}{list of numeric values (theta,c0.H,b.H0), initial values for model parameters.}
%%     ~~Describe \code{x} here~~
}
\details{The function implements the MCMC scheme described in Kim and Kottas (2022), which draws posterior samples for parameters of Erlang mixture models for NHPP intensities. The Erlang mixture model is defined as
\loadmathjax
\mjdeqn{\lambda(t) = \sum_{j=1}^J \omega_j \text{ga}(t|j,\theta^{-1}), \quad t \in \mathbb{R}^+,}{}
\mjdeqn{\omega_j = H(j\theta)-H((j-1)\theta), \quad H \sim \mathcal{G}(H_0,c_0).}{}
where \mjeqn{\mathcal{G}(H_0,c_0)}{} is a gamma process specified through \mjeqn{H_0}{}, a (parametric) cumulative intensity function, and \mjeqn{c_0}{}, a positive scalar parameter. For any \mjeqn{t \in \mathbb{R}^+}{}, \mjeqn{\text{E}(H(t))= H_0(t)}{} and \mjeqn{\text{Var}(H(t)) = H_0(t)/c_0}{}, and thus \mjeqn{H_0}{} plays the role of the centering cumulative intensity, whereas \mjeqn{c_0}{} is a precision parameter. For simplicity, the centering cumulative intensity is chosen as \mjeqn{H_0(t) = t/b}{}. Still, the random precision parameter which can take small values close to 0, ensures the model flexibility.

The model takes the following priors for the model parameters \mjeqn{(J,\theta,c_0,b)}{}: Lomax prior with shape \mjeqn{2}{}, scale \mjeqn{d_\theta}{}, and thus median \mjeqn{d_\theta(\sqrt{2}-1)}{} along with fixed \mjeqn{J}{}. Exponential priors with rates \mjeqn{a_{c_0}}{} and \mjeqn{a_{b}}{} for \mjeqn{c_0}{} and \mjeqn{b}{}, respectively. Therefore, the \mjeqn{\texttt{pri.param}}{} argument requires a lits of compoents: with fixed J,

a.theta, b.theta, sig2.theta for \mjeqn{\theta}{},

a.c0.H, sig2.c0.H for \mjeqn{c_0}{},

a.b.H0, sig2.b.H0 for \mjeqn{b}{}.

sig2.theta, sig2.c0.H, sig2.b.H0 are tuning parameterers for M-H alogirhms to sample each parameter.


%%  ~~ If necessary, more details than the description above ~~
}
\value{A list is returned with the following components:
  \item{time}{system date and time at every \mjeqn{(\texttt{nsim}/5)}{}th MCMC iteration.}
  \item{ar.theta}{ratio of the number of accepted \mjeqn{\theta}{}s to \mjeqn{\texttt{nsim}}{}.}
  \item{ar.c0.H}{ratio of the number of accepted \mjeqn{c_0}{}s to \mjeqn{\texttt{nsim}}{}.}
  \item{ar.b.H0}{ratio of the number of accepted \mjeqn{b}{}s to \mjeqn{\texttt{nsim}}{}.}
  \item{H}{matrix of sampled \mjeqn{\omega}{}.}
  \item{theta}{vector of sampled \mjeqn{\theta}{}.}
  \item{c0.H}{vector of sampled \mjeqn{c_0}{}.}
  \item{b.H0}{vector of sampled \mjeqn{b}{}.}

%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Hyotae Kim and Athanasios Kottas (2022). Erlang mixture modeling for Poisson process intensities. Statistics and Computing, 32 (1), 1-15.
}
\author{Hyotae Kim
%%  ~~who you are~~
}
\examples{
# synthetic data generation
maxT <- 2e1
param <- NULL
param$shape <- 0.5
param$scale <- 0.01
t <- sampler.nhpp(maxT,param,"D")

nsim <- 1000

# a list of (hyper)parameters
pri.param <- NULL
pri.param$J <- 50
pri.param$a.theta <- 2; pri.param$b.theta <- 1; pri.param$sig2.theta <- 0.001
pri.param$a.c0.H <- 0.1; pri.param$sig2.c0.H <- 0.15
pri.param$a.b.H0 <- 2; pri.param$sig2.b.H0 <- 0.1

# a list of initival values for each parameter
init <- NULL
init$theta <- 1; init$c0.H <- 1; init$b.H0 <- 1

### Defaults for the list, init, are set to prior means.
# post.samp <- do.mcmc.nhpp.erlmix(nsim,t,maxT,pri.param)

post.samp <- do.mcmc.nhpp.erlmix(nsim,t,maxT,pri.param,init)

# acceptance ratios of each parameter
print(c(post.samp$ar.theta,post.samp$ar.c0.H,post.samp$ar.b.H0))

# trace plots for the model parameters
par(mfrow=c(2,2))
ts.plot(post.samp$theta)
ts.plot(post.samp$c0.H)
ts.plot(post.samp$b.H0)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
