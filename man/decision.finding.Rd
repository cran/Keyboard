% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/decision.finding.R
\name{decision.finding}
\alias{decision.finding}
\title{Return Dose Escalation or De-escalation Decision Given Number of Toxicity and Responses}
\usage{
decision.finding(out.matrix, n, t, r)
}
\arguments{
\item{out.matrix}{A decision matrix.}

\item{n}{The number of enrolled subjects.}

\item{t}{The number of subjects experienced dose limiting toxicities (DLT).}

\item{r}{The number of responses.}
}
\value{
The function \code{decision.finding()} returns the decision.
}
\description{
This function returns the decision given the decision table, the total number of patients, the number of patients who experienced toxicity and the number of responses.
}
\details{
\code{decision.finding()} is a helper function to \code{get.decision.obd.kb()}  and \code{get.decision.obd2.kb()}. \code{decision.finding()} returns the decision given the decision table, the total number of patients, the number of patients who experienced toxicity and the number of responses. If the number of subjects enrolled is n, then the number of possible combination for subjects who experienced toxicity and efficacy is (n+1) x (n+1).  For n = 30, the number of entries is 961. Finding the decision given a scenario is tedious.  \code{decision.finding()} resolves this challenge by helping users to obtain the decision easily with the input of the decision table (output from \code{get.decision.obd.kb()}  and \code{get.decision.obd2.kb()}), the total number of patients, the number of patients who experienced toxicity and the number of responses.
}
\examples{
 output.matrix <- get.decision.obd.kb( toxicity.low = 0.15,
                 toxicity.moderate= 0.25,
                 toxicity.high = 0.35,
                 efficacy.low = 0.25,
                 efficacy.moderate = 0.45,
                 efficacy.high = 0.65,
                 target.toxicity=0.20,
                 target.efficacy=0.40,
                 cohortsize=3, ncohort=10)$decision.matrix
decision <- decision.finding (out.matrix=output.matrix, n=6, t=1, r=3)
decision
}
\author{
Xiaomeng Yuan, Chen Li, Hongying Sun, Li Tang and Haitao Pan
}
