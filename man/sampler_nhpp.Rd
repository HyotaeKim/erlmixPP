\name{sampler_nhpp}
\alias{sampler.nhpp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Nonhomogeneous Poisson Process Sampler
%%  ~~function to do ... ~~
}
\description{Take a sample from a nonhomogeneous Poisson process (NHPP) with a certain-formed intensity function.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{sampler.nhpp(maxT, param, type)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{maxT}{a positive number, the upper bound of an interval, (0,maxT), from which temporal sample points are taken.}
  \item{param}{parameters of a Weubll hazard or a mixture of Weibull hazard functions for the intensity function.}
  \item{type}{a character, either "D" for decreasing, "I" for increasing, or "B" for bimodal intensity.}
}
\details{Decreasing and increasing intensity functions are defined by the Weibull hazard function,
\loadmathjax
\mjdeqn{h(t) = a/b(t/b)^{a-1}}{} for \mjeqn{b > 0}{} and either \mjeqn{a < 1}{} for a decreasing or \mjeqn{a > 1}{} for an increasing intensity. The bimodal intensity function is obtained by a mixture of two Weibull hazard functions,
\mjdeqn{h(t) = w_1*a_1/b_1(t/b_1)^{a_1-1}+w_2*a_2/b_2(t/b_2)^{a_2-1}}{} for positive vectors \mjeqn{(a_1,a_2)}{}, \mjeqn{(b_1,b_2)}{}, and \mjeqn{(w_1,w_2)}{}.
%%  ~~ If necessary, more details than the description above ~~
}
\value{A realization of NHPP with the specified intensity function.
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
maxT <- 2e1
param <- NULL
param$shape <- 0.5
param$scale <- 0.01
samp <- sampler.nhpp(maxT,param,"D")
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
