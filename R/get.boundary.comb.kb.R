#' Dose Escalation or De-escalation Boundaries for Drug-combination Trials
#'
#' This function generates the optimal dose escalation or de-escalation boundaries when
#' conducting a drug-combination trial with the Keyboard design.
#'
#' @details
#' The Keyboard design relies on the posterior distribution of the toxicity
#' probability to guide dosage. To determine whether to escalate or de-escalate the dose, given the observed data at the current dose, we first identify an
#' interval that has the highest posterior probability, referred to as
#' the "strongest key". This key represents where the true dose-limiting
#' toxicity (DLT) rate of the current dose is most likely located. If the
#' strongest key is to the left of the "target key", then we escalate
#' the dose because the data suggest that the current
#' dose is likely to underdose patients; if the strongest key is
#' to the right of the target key, then we de-escalate the dose
#' because the observed data suggest that the current dose is likely to overdose the patients; and
#' if the strongest key is the target key, then we retain the current dose because
#' the observed data support that the current dose is most likely to be in the
#' proper dosing interval.
#' Graphically, the strongest key is the one with the largest area under the
#' posterior distribution curve of the DLT rate of the current dose.
#'
#' \figure{Keyboard.jpg}
#' 
#' An attractive feature of the Keyboard design is that its dose escalation and
#' de-escalation rules can be tabulated before the onset of the trial. Thus,
#' when conducting the trial, no calculation or model fitting is needed, and we
#' need to count only the number of DLTs observed at the current dose; 
#' the decision to escalate or de-escalate the dose is based on the pre-tabulated
#' decision rules.
#'
#' Given all observed data, we use matrix isotonic regression to obtain an
#' estimate of the toxicity rate of the combination of dose level j of drug A
#' and dose level k of drug B and  to  select as the MTD  the combination with the
#' toxicity estimate that is closest to the target. When there are ties, we
#' randomly choose one as the MTD.
#'
#' For patient safety, we apply the following Bayesian overdose control rule
#' after each cohort:
#' if at least 3 patients have been treated at the given dose and
#' the observed data indicate that the probability of the current combination dose's toxicity rate being above the target toxicity rate is more
#' than 95\%,  then we exclude this dose and beyond to avoid 
#' exposing future patients to these overly toxic doses. The probability
#' threshold can be specified with \code{cutoff.eli}. If the lowest dose
#' combination (1, 1) is overly toxic, then the trial terminates early, and no dose
#' is selected as the MTD.
#'
#' @param target The target dose-limiting toxicity (DLT) rate.
#' @param ncohort A scalar specifying the total number of cohorts in the trial.
#' @param cohortsize The number of patients in the cohort.
#' @param marginL The difference between the target and the lower bound of the
#'                "target key" (proper dosing interval) to be defined.\cr
#'                The default is 0.05.
#' @param marginR The difference between the target and the upper bound of the
#'                "target key" (proper dosing interval) to be defined.\cr
#'                The default is 0.05.
#' @param cutoff.eli The cutoff to eliminate an overly toxic dose and all
#'                   higher doses for safety.\cr
#'                   The recommended value is 0.95.
#' @param n.earlystop The early stopping parameter. If the number of patients treated at
#'                    the current dose reaches \code{n.earlystop}, then stop the trial
#'                    and select the MTD based on the observed data. The default
#'                    value is 100.
#' @param extrasafe Set \code{extrasafe=TRUE} to impose a stricter stopping rule for extra safety, expressed as the stopping boundary value in the result.
#' @param offset  A small positive number (between 0 and 0.5) to control how strict
#'               the stopping rule is when \code{extrasafe=TRUE}. A larger value leads
#'               to a stricter stopping rule. The default value is 0.05.
#'
#' @return The function returns a matrix, including the dose escalation, de-escalation, and elimination boundaries.
#'
#' @note In most clinical applications, the target DLT rate is often a rough
#'   guess, but finding a dose level with a DLT rate reasonably close to the
#'   target rate (which ideally would be the MTD) is what interests the
#'   investigator.
#'
#' @examples
#' ### Drug-combination trial ###
#'
#' bound <- get.boundary.comb.kb(target=0.3, ncohort=10, cohortsize=3)
#' print(bound)
#' @family drug-combination functions
#'
#' @references
#'
#' Yan F, Mandrekar SJ, Yuan Y. Keyboard: A Novel Bayesian Toxicity Probability
#' Interval Design for Phase I Clinical Trials.
#' \emph{Clinical Cancer Research}. 2017; 23:3994-4003.
#' http://clincancerres.aacrjournals.org/content/23/15/3994.full-text.pdf
#' 
#' 
#' Pan H, Lin R, Yuan Y. Keyboard design for phase I drug-combination trials.
#' \emph{Contemporary Clinical Trials}. 2020.
#' https://doi.org/10.1016/j.cct.2020.105972
#'@export
get.boundary.comb.kb <- function(target, ncohort, cohortsize,n.earlystop=100,
                                 marginL=0.05, marginR=0.05, cutoff.eli=0.95,offset=0.05,extrasafe=TRUE) {
  ## Get cutoffs for keys
  getkeys <- function(a1, a2)
  {
    delta=a2-a1
    lkey=NULL; rkey=NULL
    i=0
    cutoff=0.3
    while (cutoff>0)
    {
      i=i+1
      cutoff = a1-i*delta
      lkey = c(cutoff, lkey)
    }
    lkey[lkey<0]=0

    i=0; cutoff=0.3
    while (cutoff<1)
    {
      i=i+1
      cutoff = a2+i*delta
      rkey = c(rkey, cutoff)
    }
    rkey[rkey>1]=1
    key=c(lkey, a1, a2, rkey)

    return(key)
  }

  ## Identify the key with the highest posterior tocutoff.elicity probability
  keys = getkeys(target - marginL, target + marginR)
  nkeys = length(keys) - 1

  npts = ncohort * cohortsize
  targetp = rep(NULL, nkeys)

  a = b = 1  # Hyperparameters used in beta prior
  decision_table = matrix(NA, nrow = npts + 1, ncol = npts)
  for (ntr in 1:npts) {
    elim = 0
    for (ntox in 0:ntr) {

        if (ntox >= 3 & 1 - pbeta(target, ntox + a, ntr - ntox + b) > cutoff.eli) {
          elim = 1
          break
        }

      # targetp = rep(NULL, nkeys)
      for(i in 1:nkeys) {
        comp = 1
        if (i == 1 || i == nkeys) {
          # Compensation factor for incompleted keys:
          comp = (marginL + marginR) / (keys[i+1] - keys[i])
        }
        targetp[i] = (pbeta(keys[i+1], ntox + a, ntr - ntox + b) -
                      pbeta(keys[i], ntox + a, ntr - ntox + b)) * comp
      }
      highkey = max(which(targetp == max(targetp)))
      targetkey = which(keys == (target - marginL))

      if (highkey > targetkey) {
        decision_table[ntox+1, ntr] <- "D"
      }
      if (highkey == targetkey) {
        decision_table[ntox+1, ntr] <- "S"
      }
      if (highkey < targetkey) {
        decision_table[ntox+1, ntr] <- "E"
      }
    }
    if (elim == 1) {
      decision_table[(ntox+1):(ntr+1), ntr] <- rep("DU", ntr - ntox + 1)
    }
  }
  colnames(decision_table) <- 1:npts
  rownames(decision_table) <- 0:npts

  boundary = matrix(NA, nrow = 4, ncol = npts)
  boundary[1, ] = 1:npts
  for (i in 1:npts) {
    if (length(which(decision_table[, i] == "E"))) {
      boundary[2, i] = max(which(decision_table[, i] == "E")) - 1
    }
    else {
      boundary[2, i] = -1
    }
    if (length(which(decision_table[, i] == "D"))) {
      boundary[3, i] = min(which(decision_table[, i] == "D")) - 1
    }
    else if (length(which(decision_table[, i] == "DU"))) {
      boundary[3, i] = min(which(decision_table[, i] == "DU")) - 1
    }
    if (length(which(decision_table[, i] == "DU"))) {
      boundary[4, i] = min(which(decision_table[, i] == "DU")) - 1
    }
  }
  colnames(boundary) <- c(rep("", npts))
  rownames(boundary) <- c("Number of patients treated",
                          "Escalate if # of DLT <=",
                          "de-escalate if # of DLT >=",
                          "Eliminate if # of DLT >=")


  if (extrasafe) {
    stopbd = NULL
    ntrt = NULL
    for (n in 1:npts) {
      ntrt = c(ntrt, n)
      if (n < 3) {
        stopbd = c(stopbd, "NA")
      }
      else {
        for (ntox in 1:n) {
          if (1 - pbeta(target, ntox + 1, n - ntox +
                        1) > cutoff.eli - offset) {
            stopneed = 1
            break
          }
        }
        if (stopneed == 1) {
          stopbd = c(stopbd, ntox)
        }
        else {
          stopbd = c(stopbd, "NA")
        }
      }
    }
    stopboundary = data.frame(rbind(ntrt, stopbd)[, 1:min(npts, n.earlystop)])
    rownames(stopboundary) = c("The number of patients treated at the lowest dose  ",
                               "Stop the trial if # of DLT >=        ")
    stopboundary = data.frame(stopboundary)
    #colnames(stopboundary) = rep("", min(npts, n.earlystop))
    colnames(stopboundary) = 1:dim(stopboundary)[2]

  }



 if(extrasafe){
   return(list(boundary=boundary,safe=stopboundary))
 }else{
   return(list(boundary=boundary))
 }


}

###  obtain dose escalation/de-escalation boundaries up to n=20 for conducting the trial
get.boundary.comb.kb(target=0.2, ncohort=8, cohortsize=3, marginL=0.05, marginR=0.05,
                    cutoff.eli=0.95, offset=0.05,extrasafe=TRUE)
get.boundary.comb.kb(target=0.2, ncohort=8, cohortsize=3, marginL=0.05, marginR=0.05,
                     cutoff.eli=0.95, offset=0.05,extrasafe=FALSE)



#   target=0.2; ncohort=8; cohortsize=3; marginL=0.05; marginR=0.05; cutoff.eli=0.95; offset=0.05;extrasafe=TRUE
