#' Dose Escalation or De-escalation Boundaries for Single-agent Trials
#'
#' This function generates the optimal dose escalation or de-escalation boundaries when
#' conducting a single-agent trial with the Keyboard design.
#'
#' @details
#' The Keyboard design relies on the posterior distribution of the toxicity
#' probability to guide dosage. To make a decision of dose escalation or
#' de-escalation, given the observed data at the current dose, we first identify an
#' interval that has the highest posterior probability, referred to as
#' the "strongest key". This key represents where the true dose-limiting
#' toxicity (DLT) rate of the current dose is most likely located. If the
#' strongest key is to the left  of the "target key", then we escalate
#' the dose because the data suggest that the current
#' dose is most likely too low; if the strongest key is
#' to the right of the target key, then we de-escalate the dose
#' because the observed data suggest that the current dose is likely too high; and
#' if the strongest key is the target key, then we retain the current dose because
#' the observed data support the notion that the current dose is most likely to be in the
#' proper dosing interval.
#' Graphically, the strongest key is the one with the largest area under the
#' posterior distribution curve of the DLT rate of the current dose.
#' 
#' \figure{Keyboard.jpg}
#'
#' An attractive feature of the Keyboard design is that its dose escalation and
#' de-escalation rules can be tabulated before the onset of the trial. Thus,
#' when conducting the trial, no calculation or model fitting is needed, and we
#'  need to count only the number of DLTs observed at the current dose; 
#' the decision to escalate or de-escalate the dose is based on the pre-tabulated
#' decision rules.
#'
#' Given all observed data, the Keyboard design uses isotonic regression to obtain an efficient statistical estimate of
#' the maximum tolerated dose (MTD) by utilizing the fact that toxicity
#' presumably increases with the dose.
#'
#' For patient safety, we apply the following Bayesian overdose control rule
#' after each cohort:
#' if at least 3 patients have been treated at the given dose and
#' the observed data indicate that the probability of the current dose's toxicity rate being above the target toxicity rate is more
#' than 95\%, then we stop the trial to avoid 
#' exposing future patients to these overly toxic doses. The probability
#' threshold can be specified with \code{cutoff.eli}. When a dose is
#' eliminated, the design recommends the next lower dose for treating the next cohort of patients. If the lowest dose is overly toxic, then the trial terminates early and no dose is selected as the MTD.
#'
#' @param target The target dose-limiting toxicity (DLT) rate.
#' @param marginL The difference between the target and the lower bound of the
#'                "target key" (proper dosing interval) to be defined.\cr
#'                The default is 0.05.
#' @param marginR The difference between the target and the upper bound of the
#'                "target key" (proper dosing interval) to be defined.\cr
#'                The default is 0.05.
#' @param cutoff.eli The cutoff value to eliminate an overly toxic dose and all
#'                   higher doses for safety.\cr
#'                   The recommended value is 0.95.
#' @param ncohort The total number of cohorts.
#' @param cohortsize The number of patients in the cohort.
#' @param n.earlystop The early stopping parameter. If the number of patients treated at
#'                    the current dose reaches \code{n.earlystop}, then stop the trial
#'                    and select the MTD based on the observed data. The default
#'                    value is 100.
#' @param extrasafe Set \code{extrasafe=TRUE} to impose a stricter stopping rule for extra safety, expressed as the stopping boundary value in the result.
#' @param offset  A small positive number (between 0 and 0.5) to control how strict
#'               the stopping rule is when \code{extrasafe=TRUE}. A larger value leads
#'               to a stricter stopping rule. The default value is 0.05.
#' @return The function returns a matrix, which includes the dose escalation,
#'   de-escalation, and elimination boundaries.
#'
#' @note In most clinical applications, the target DLT rate is often a rough
#'   guess, but finding a dose level with a DLT rate reasonably close to the
#'   target rate (which ideally would be the MTD) is of interest.
#'
#' @examples
#' ### Single-agent trial ###
#'
#' bound <- get.boundary.kb(target=0.3, ncohort=10, cohortsize=3)
#' print(bound)
#' @family single-agent functions
#'
#' @references
#'
#' Yan F, Mandrekar SJ, Yuan Y. Keyboard: A Novel Bayesian Toxicity Probability
#' Interval Design for Phase I Clinical Trials.
#' \emph{Clinical Cancer Research}. 2017; 23:3994-4003.
#' http://clincancerres.aacrjournals.org/content/23/15/3994.full-text.pdf
#' @export
get.boundary.kb <- function(target, ncohort, cohortsize, marginL=0.05, marginR=0.05,cutoff.eli=0.95, n.earlystop=100, extrasafe=FALSE, offset=0.05) {
  ### simple error checking
  if(target<0.05) {warning("Error: the target is too low! \n"); return();}
  if(target>0.6)  {warning("Error: the target is too high! \n"); return();}
  if(offset>=0.5) {warning("Error: the offset is too large! \n"); return();}
  if(n.earlystop<=6) {warning("Warning: the value of n.earlystop is too low to ensure good operating characteristics. Recommend n.earlystop = 9 to 18 \n"); return();}

    ## get cutoffs for keys
    getkey <- function(a1, a2)
    {
        delta=a2-a1
        lkey=NULL; rkey=NULL
        i=0; cutoff=0.3
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

    ## Identify the key with the highest posterior to cutoff.elicity probability
    keys = getkey(target-marginL, target+marginR)
    nkeys = length(keys)-1

    targetp = NULL
    npts = ncohort*cohortsize

    a=1; b=1;  # hyperparameters used in beta prior
    decision_table = matrix(NA,nrow=npts+1,ncol=npts)
    for (ntr in 1:npts)
    {
        elim=0
        for (ntox in 0:ntr)
        {
            if (1-pbeta(target, ntox+a, ntr-ntox+b)>cutoff.eli) {elim=1; break;}

            for (i in 1:nkeys)
            {
                comp = 1;
                if (i==1 || i==nkeys) {
                    comp = (marginL+marginR)/(keys[i+1]-keys[i])  #compensation factor for incompleted keys
                }
                targetp[i] = (pbeta(keys[i+1], ntox+a, ntr-ntox+b) - pbeta(keys[i], ntox+a, ntr-ntox+b))*comp
            }
            highkey = max(which(targetp==max(targetp)))
            targetkey = which(keys==(target-marginL))

            if (highkey>targetkey)  {decision_table[ntox+1,ntr] <- "D"}
            if (highkey==targetkey) {decision_table[ntox+1,ntr] <- "S"}
            if (highkey<targetkey)  {decision_table[ntox+1,ntr] <- "E"}
        }
        if (elim==1) {
            decision_table[(ntox+1):(ntr+1),ntr] <- rep("DU",ntr-ntox+1)
        }
    }
    colnames(decision_table) <- 1:npts
    rownames(decision_table) <- 0:npts

    boundary = matrix(NA, nrow=4, ncol=npts)
    boundary[1,] = 1:npts
    for (i in 1:npts)
    {
        if (length(which(decision_table[,i]=="E"))) {
            boundary[2,i] = max(which(decision_table[,i]=="E"))-1
        }
        else {
            boundary[2,i] = -1
        }
        if (length(which(decision_table[,i]=="D"))) {
            boundary[3,i] = min(which(decision_table[,i]=="D"))-1
        }
        else if (length(which(decision_table[,i]=="DU"))) {
            boundary[3,i] = min(which(decision_table[,i]=="DU"))-1
        }
        if (length(which(decision_table[,i]=="DU"))) {
            boundary[4,i] = min(which(decision_table[,i]=="DU"))-1
        }
    }
    colnames(boundary) <- c(rep("", npts))
    rownames(boundary) <- c("Number of patients treated",
                            "Escalate if # of DLT <=",
                            "de-escalate if # of DLT >=",
                            "Eliminate if # of DLT >=")
    out = list( boundary_tab=boundary[, (1:floor(min(npts, n.earlystop)/cohortsize))*cohortsize],
      full_boundary_tab=boundary)


    ## if extrasafe, add more info, like in BOIN function
    if (extrasafe){
      stopbd=NULL;
      ntrt=NULL;
      for ( n in 1:npts){
        ntrt = c(ntrt, n);
        if (n <3){ stopbd = c(stopbd, NA);}
        else{
          for (ntox in 1:n){
            if (1-pbeta(target, ntox+1, n-ntox+1)>cutoff.eli-offset) {stopneed=1; break;}
          }
          if (stopneed==1) {stopbd=c(stopbd, ntox);} else (stopbd=c(stopbd,NA))

        }
      }
      stopboundary = rbind(ntrt, stopbd)[, 1:min(npts, n.earlystop)]
      rownames(stopboundary) = c("The number of patients treated at the lowest dose  ", "Stop the trial if # of DLT >=        ");
      colnames(stopboundary) = rep("", min(npts, n.earlystop));
      out = c(out,list(target=target, cutoff=cutoff.eli-offset, stop_boundary=stopboundary))
  }

    return(out)
}
