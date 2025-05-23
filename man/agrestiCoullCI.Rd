% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/agrestiCoullCI.R
\name{agrestiCoullCI}
\alias{agrestiCoullCI}
\title{Agresti-Coull confidence interval for a binomial proportion}
\usage{
agrestiCoullCI(n1, n0, q)
}
\arguments{
\item{n1}{number of successes/ones}

\item{n0}{number of failures/zeroes}

\item{q}{quantile for eventual CI (e.g. 0.95 for a 95 percent binomial CI)}
}
\value{
the approximate (q x 100) percent confidence interval for (p|n1,n0,q)
}
\description{
Agresti-Coull confidence interval for a binomial proportion
}
\details{
\eqn{z_\alpha = \Phi^{-1}(1 - \frac{\alpha}{2})}

\eqn{\tilde{n} = n_{\text{successes}} + n_{\text{failures}} + z^2_\alpha}

\eqn{\tilde{p} = \frac{1}{\tilde{n}}(n_{\text{success}} + \frac{z^2_\alpha}{2})}

\eqn{p \approx \tilde{p} \pm z_\alpha \times \sqrt{\frac{\tilde{p}}{\tilde{n}} \times (1 - \tilde{p})}}
}
\examples{
agrestiCoullCI(10, 3, 0.95)
}
