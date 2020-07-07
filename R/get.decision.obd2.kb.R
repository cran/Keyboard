#' Automatic Boundary and Decision Tables for Phase I/II Trials.
#' 
#' This function automatically generates a boundary table and a decision table for single-agent Phase I/II trials using Keyboard design.
#'
#' @details 
#' The Keyboard design utilizes the posterior distributions of the toxicity and efficacy to guide dosage transition. To determine whether to escalate or de-escalate given the observed data at the current dose, we must generate a boundary table and a decision matrix. The boundary table requires a lower boundary and an upper boundary for both toxicity and efficacy and a decision matrix, which provides dose escalation or de-escalation suggestions. 
#' 
#' A Bayesian toxicity probability interval design from Yuan's group is used to determine the lower and upper boundaries for toxicity  based on the target toxicity rate. 
#' 
#' The boundaries for efficacy are calculated based on the efficacy failure rate, which is \eqn{1-target.efficacy}.  The efficacy failure rate is defined as 1 - target.efficacy and  treated similarly to toxicity. The lower and upper boundaries for efficacy failure are determined using a Bayesian toxicity probability interval design.  To compute the lower and upper boundaries for efficacy, we use the following:
#' \deqn{efficacy.lower.boundary = 1 - efficacy.failure.upper.boundary}
#'  \deqn{efficacy.upper.boundary = 1 - efficacy.failure.lower.boundary}
#'  
#' There are 3 subintervals for toxicity:
#'  \enumerate{
#' \item The low interval for toxicity is (0, toxicity.lower.boundary). 
#' \item The moderate interval for toxicity is (toxicity.lower.boundary, toxicity.upper.boundary).
#' \item The high interval for toxicity is (toxicity.upper.boundary, 1).
#' }
#' There are 3 subintervals for efficacy:
#'  \enumerate{
#' \item The low interval for efficacy is (0, efficacy.lower.boundary). 
#' \item The moderate interval for efficacy is (efficacy.lower.boundary, efficacy.upper.boundary).
#' \item The high interval for efficacy is (efficacy.upper.boundary, 1).
#' }
#' Assuming that the toxicity probability is \eqn{p_i} and the efficacy probability is \eqn{q_i} at dose level i, the probability unit intervals for toxicity and efficacy can be partitioned in to subintervals (a, b) and (c,d). (a, b)  is a subinterval for toxicity probability and (c, d) is a subinterval for efficacy probability. (a, b) x (c,d) is a combination interval. 3 toxicity intervals and 3 efficacy intervals result in 9 combination intervals. Besides the target toxicity rate and the target efficacy rate, a matrix is required that contains 9 decisions for these 9 combination intervals ordered as ((low toxicity, low efficacy),(low toxicity, moderate efficacy),(low toxicity, high efficacy), (moderate toxicity, low efficacy),(moderate toxicity, moderate efficacy),(moderate toxicity, high efficacy),(high toxicity, low efficacy),(high toxicity, moderate efficacy),(high toxicity, high efficacy)).  Possible decisions are "E", "D","S". Decision "D" denotes de-escalation, so the next cohort of patients will be treated at the next lower dose level.  Decision "E" denotes escalation, so the next cohort of patients will be treated at the next higher dose level.  Decision "S" denotes stay, so the next cohort of patients will be treated at the current dose level. 
#' 
#' Shown here is an example of boundary tables designed using this method with a target toxicity rate of 0.2 and a target efficacy rate of 0.4:
#' 
#' \tabular{rrrrrr}{
#'        \tab   \tab Efficacy.low \tab Efficacy.moderate \tab Efficacy.high  \cr
#'        \tab   \tab  (0,0.27) \tab  (0.27,0.52) \tab (0.52,1) \cr
#'   Toxicity.low \tab  (0,0.16)\tab  E\tab  E\tab  S \cr
#'   Toxicity.moderate \tab  (0.16,0.24)\tab  S\tab  S \tab S  \cr
#'   Toxicity. high \tab  (0.24,1)\tab  D\tab  D\tab  D \cr

#' }
#' For example, the interval combination (0, 0.16) x (0.27,0.52) corresponds to a decision "S". This means that the next cohort of patients will be treated at the current dose level if the observed toxicity rate of current dose falls in (0, 0.16)  and the observed efficacy rate falls in (0.27,0.52) . 
#' 
#' Suppose that there are d doses in the trial and that the current dose is i. Define the number of patients as \eqn{n_i}, the number of patients who experienced toxicity as \eqn{x_i} and the number of responses as \eqn{y_i}.  The trial data can be shown as follows:
#' \deqn{D= {(n_i, x_i, y_i), i = 1, ..., d}} 
#' 
#' Bayesian rule is used to calculate the joint unit probability mass (JUPM) for the toxicity and efficacy combination intervals. For a given combination interval, JUPM is calculated as follows:
#' \deqn{JUPM_{(a,b)}^{(c,d)} =Pr{p_j \in (a,b), q_j \in (c,d) | D} / (b-a)*(d-c) }
#' \eqn{Pr{p_j \in (a,b), q_j \in (c,d) | D} } are the posterior probabilities of \eqn{p_i} and \eqn{q_i} falling in the subinterval (a,b) and (c,d).  Assume the prior for both \eqn{p_i} and \eqn{q_i} follows independent beta distributions \eqn{beta(\alpha_p,\beta_p)} and \eqn{beta(\alpha_q,\beta_q)} respectively. The posterior distributions for \eqn{p_i} and \eqn{q_i} are \eqn{beta(\alpha_p + x_i,\beta_p + n_i - x_i)} and \eqn{beta(\alpha_q + y_i,\beta_q + n_i - y_i)}. Using these posterior distributions, calculate the JUPM for all  16 combination intervals and find the winning combination interval (a*,b*) and (c*,d*), which achieves the largest JUPM. The decision that corresponds to this winning combination interval is used to treat the next cohort of patients.
#' 
#' 
#' Two dose exclusion rules are applied in the trial design: safety rule and futility rule.
#' Safety rule:
#' if at least 3 patients have been treated at a given dose and given 
#' the observed data indicate that the probability of the toxicity rate of
#' the current dose exceeding the target toxicity rate is more
#' than 95\%, we eliminate the current and any higher dose from the trial to prevent
#' exposing future patients to with unacceptable high toxicity. The probability
#' threshold can be specified with \code{cutoff.eli.toxicity}. When a dose is
#' eliminated, the design recommends the next lower dose for treating the next cohort of 
#' patients. If the lowest dose is overly toxic, then the trial terminates early and
#' no dose is selected as the OBD. This corresponds to a dose assignment of "DUT", which means de-escalation due to unacceptable high toxicity and excluding the current dose and any dose higher than this dose from the trial.
#' 
#' Futility rule:
#' if at least 3 patients have been treated at a given dose and 
#' the observed data indicate that the probability of the  current dose's efficacy rate  exceeding the target efficacy rate is less than 30\%, then we eliminate this dose from the trial  to avoid
#' exposing future patients to these futile doses. The probability
#' threshold can be specified with \code{cutoff.eli.efficacy}.  This corresponds to two possible dose assignments: "EUE" and "DUE", both excluding this dose from the trial. "EUE" denotes escalation due to unacceptable low efficacy and "DUE" denotes de-escalation due to unacceptable low efficacy.
#' 
#' An attractive feature of the Keyboard design is that its dose escalation and
#' de-escalation rule can be tabulated before implementing the trial. Thus,
#' when conducting the trial, no real-time calculation or model fitting is needed, and one 
#'  needs to count only the number of patients, the number of DLTs, and the number of responses observed at the current dose and the desicion to escalate or de-escalate is  based on the pre-tabulated
#' decision rules.
#'
#' @param target.toxicity  The target toxicity rate.
#' @param target.efficacy  The target efficacy rate.
#' @param cohortsize  The number of patients in the cohort.
#' @param ncohort The total number of cohorts.
#' @param cutoff.eli.toxicity The cutoff value to eliminate a dose with unacceptable high toxicity for safety. The default value is 0.95.
#' @param cutoff.eli.efficacy The cutoff value to eliminate a dose with unacceptable low efficacy. The default value is 0.3.
#' @param decision The pre-specified decision matrix for the 9 combination intervals: (low toxicity, low efficacy),(low toxicity, moderate efficacy),(low toxicity, high efficacy), (moderate toxicity, low efficacy),(moderate toxicity, moderate efficacy),(moderate toxicity, high efficacy),(high toxicity, low efficacy),(high toxicity, moderate efficacy),(high toxicity, high efficacy)
#' @return \code{get.decision.obd2.kb()} returns a prespecified boundary table and a dose assignment decision table:
#' \enumerate{
#' \item Boundary table (\code{$boundary.table})
#' \item Decision matrix (\code{$decision.matrix})
#' 
#' }
#' @author Hongying Sun, Li Tang, and Haitao Pan
#' @examples
#' \donttest{
#' decision.obd2.kb <- get.decision.obd2.kb(target.toxicity=0.2,
#'            target.efficacy=0.4, cohortsize=3, ncohort=10)
#'  print(decision.obd2.kb)
#' }
#' @family single-agent phase I/II functions
#' 
#'
#' @references 
#' Liu S. and Yuan, Y.  Bayesian Optimal Interval Designs for Phase I Clinical Trials, \emph{Journal of the Royal Statistical Society: Series C}. 2015; 64, 507-523.
#'
#' Yuan Y., Hess K.R., Hilsenbeck S.G. and Gilbert M.R.  Bayesian Optimal Interval Design: A Simple and Well-performing Design for Phase I Oncology Trials. \emph{Clinical Cancer Research}. 2016; 22, 4291-4301.
#' 
#' 
#'  Li DH, Whitmore JB, Guo W, Ji Y.  Toxicity and efficacy probability interval design for phase I adoptive cell therapy dose-finding clinical trials.
#' \emph{Clinical Cancer Research}. 2017; 23:13-20.
#'https://clincancerres.aacrjournals.org/content/23/1/13.long 
#' 
#' @export
get.decision.obd2.kb <- function(target.toxicity, target.efficacy, cohortsize, ncohort,cutoff.eli.toxicity= 0.95, cutoff.eli.efficacy=0.3, decision=c("E","E","S","S","S","S","D","D","D")){
  ## simple error checking
  if(target.toxicity<0.05) {warning("Error: the target toxcicity rate is too low! \n"); return();}
  if(target.toxicity>0.6)  {warning("Error: the target toxicity rate is too high! \n"); return();}
  if(target.efficacy<= 0.19) {warning("Error: the target efficacy rate is too low! \n");  return();}

  #write function to generate boundaries, can't use BOIN package because BOIN package restricts the toxicity not greater than 0.6.
  #finds the boundary
  find.boundary <- function(target,  p.saf=0.6*target, p.tox=1.4*target){

    lambda1  = log((1-p.saf)/(1-target))/log(target*(1-p.saf)/(p.saf*(1-target)))
    lambda2  = log((1-target)/(1-p.tox))/log(p.tox*(1-target)/(target*(1-p.tox)))
    out <- list(lambda_e = lambda1, lambda_d=lambda2)
  }

  find.efficacy.boundary <- function(target.efficacy,  p.saf=0.6*(1-target.efficacy), p.tox=1.4*(1-target.efficacy)){

    gamma2  = 1-log((1-p.saf)/(1-(1-target.efficacy)))/log((1-target.efficacy)*(1-p.saf)/(p.saf*(1-(1-target.efficacy))))
    gamma1  = 1-log((1-(1-target.efficacy))/(1-p.tox))/log(p.tox*(1-(1-target.efficacy))/((1-target.efficacy)*(1-p.tox)))
    out <- list(gamma_e = gamma1, gamma_d=gamma2)
  }
  #generate decision table
  #get boundary for toxicity
  bound.toxicity <- find.boundary(target=target.toxicity)
  toxicity.e <- round(bound.toxicity$lambda_e,2)
  toxicity.d <- round(bound.toxicity$lambda_d,2)


  # get boundary for efficacy
  bound.efficacy <- find.efficacy.boundary(target=target.efficacy)
  efficacy.e <- round(bound.efficacy$gamma_e,2)
  efficacy.d <- round(bound.efficacy$gamma_d,2)





  #create toxicity matrix
  toxicity.matrix<- matrix(c(0,toxicity.e, 0,toxicity.e,0,toxicity.e,toxicity.e, toxicity.d,toxicity.e, toxicity.d,toxicity.e, toxicity.d,  toxicity.d,1, toxicity.d,1,toxicity.d,1),ncol=2, byrow=TRUE)
  #create efficacy failure matrix
  #create decision matrix
  decision.matrix <- matrix(decision, ncol=1)


  #create efficacy failure matrix
  efficacy.matrix <- matrix(rep(c(0,efficacy.e,efficacy.e,efficacy.d,efficacy.d,1 ), 3), ncol=2, byrow=TRUE)
  #create decision table

  decision.table <- cbind(toxicity.matrix,efficacy.matrix, decision.matrix)
  colnames(decision.table) <- c("T1", "T2", "EF1","EF2","Decision")
  #create output tables to store results
  #to store the total number of patients.
  ntvector <- c()
  ndvector <-  c()
  nrvector <-  c()
  for (i in 1:ncohort) {
    nt<-rep(i*cohortsize,(i*cohortsize+1)^2)
    ntvector <- append(ntvector, nt)
    for (j in 0:(i*cohortsize)){
      nd <- rep(j, (i*cohortsize+1))
      #print(nd)
      ndvector <- append(ndvector, nd)
      nr <- c(0:(i*cohortsize))
      nrvector <- append(nrvector,nr)
    }

  }
  ntmatrix <- matrix(ntvector,ncol=1)
  ndmatrix <- matrix(ndvector,ncol=1)
  nrmatrix <- matrix(nrvector,ncol=1)


  #output table matrix
  outmatrix <- cbind(ntmatrix,ndmatrix,nrmatrix)
  colnames(outmatrix)<- c("N", "T","R")
  outmatrix <- as.data.frame(outmatrix)

  #this function is to return the total number of patients, the number of patients with DLT and the number of patients with reponse given a matrix and an index.
  pinfor <- function(mout, i){
    out <- list(  n = mout[i,1],
                  d = mout[i,2],
                  r = mout[i,3])
    # return (out)
  }
  # given decision table and an index, output the intervals for the toxicity and efficacy and also the dicisions.
  pdecision <- function(decisont, i){
    out <- list (
      t1 = as.numeric(decisont[i,1]),
      t2 = as.numeric(decisont[i,2]),
      e1 = as.numeric(decisont[i,3]),
      e2 = as.numeric(decisont[i,4]),
      d = decisont[i,5]

    )
  }
  # finds the decision that provides the largest jupm
  maxjupm <- function(decision.table, outmatrix, i){
    outputoutm <- pinfor(outmatrix,i)
    nt = outputoutm$n
    dt = outputoutm$d
    rt = outputoutm$r
    jupmv <- c()
    for (rowd in 1:nrow(decision.table)){
      #outdecision <- pdecision(decision.table,rowd)
      a = as.numeric(decision.table[rowd,1])
      b = as.numeric(decision.table[rowd,2])
      c = as.numeric(decision.table[rowd,3])
      d = as.numeric(decision.table[rowd,4])
      #D = outdecision$d
      jupm = round((pbeta(b, (1+dt), (1+nt-dt))-pbeta(a, (1+dt), (1+nt-dt)))*(pbeta(d, (1+rt), (1+nt-rt))-pbeta(c, (1+rt), (1+nt-rt)))/((b-a)*(d-c)),2)
      jupmv <- append(jupmv, jupm)
    }
    jupmv_matrix <- matrix(jupmv, ncol=1)
    colnames(jupmv_matrix)<- "jupm"
    decision.table_jupm <- cbind(decision.table,jupmv_matrix)
    decision.table_jupm <- as.data.frame(decision.table_jupm)
    decision.table_jupm$jupm <- as.numeric(as.character(decision.table_jupm$jupm))
    jupm_max <-max(decision.table_jupm$jupm)
    max_index <- which.max(decision.table_jupm$jupm)
    decision <- decision.table_jupm$D[max_index]

    out <- list(decision.table_jupm=decision.table_jupm, chosen_decision=decision, jupm_max_value =jupm_max)

  }

  #select decision
  decision_picking <-c()
  jupm_picking <-c()
  for (rowoutm in 1:nrow(outmatrix)){
    outjupm <- maxjupm(decision.table, outmatrix, rowoutm)
    decision_p <- outjupm$chosen_decision
    decision_picking <- append(decision_picking,as.character(decision_p))
    jupm_mv <- outjupm$jupm_max_value
    jupm_picking <- append(jupm_picking,jupm_mv)

  }
  decision_matrix <- matrix(decision_picking, ncol=1)
  colnames(decision_matrix)<- "Decision"
  jupm_max_matrix <- matrix(jupm_picking, ncol=1)
  colnames(jupm_max_matrix) <- "Jupm"

  output.matrix <- cbind(outmatrix, decision_matrix,jupm_max_matrix)
  # remove jupm information.
  output.matrix <- output.matrix[,-5]
  output.matrix <- as.data.frame(output.matrix)
  i <- sapply(output.matrix, is.factor)
  output.matrix[i] <- lapply(output.matrix[i], as.character)
  for (rowindex in 1:nrow(output.matrix)){
    if (output.matrix$Decision[rowindex] == "D") {
      if (1- pbeta(target.toxicity,output.matrix$T[rowindex]+1, output.matrix$N[rowindex]-output.matrix$T[rowindex]+1) > cutoff.eli.toxicity) {
        output.matrix$Decision[rowindex] = "DUT"

      } else if (1- pbeta(target.efficacy,output.matrix$R[rowindex]+1, output.matrix$N[rowindex]-output.matrix$R[rowindex]+1)< cutoff.eli.efficacy){
        output.matrix$Decision[rowindex] = "DUE"

      }
    } else if (output.matrix$Decision[rowindex] == "E"){
      if (1- pbeta(target.efficacy,output.matrix$R[rowindex]+1, output.matrix$N[rowindex]-output.matrix$R[rowindex]+1)< cutoff.eli.efficacy){
        output.matrix$Decision[rowindex] = "EUE"
      }

    }
  }



  outlist =list(
        boundary.table=decision.table,
        decision.matrix=output.matrix)

  return (outlist)

}






