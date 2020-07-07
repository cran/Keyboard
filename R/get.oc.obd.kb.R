#' Generate Operating Characteristics to Find the Optimal Biological Dose
#'
#' This function generates operating characteristics to find the  optimal biological dose (OBD). 
#' 
#' @param toxicity.low The upper boundary for the low toxicity interval.
#' @param toxicity.moderate The upper boundary for the moderate toxicity interval.
#' @param toxicity.high The upper boundary for the high toxicity interval.
#' @param efficacy.low The upper boundary for the low efficacy interval.
#' @param efficacy.moderate The upper boundary for the moderate efficacy interval.
#' @param efficacy.high  The upper boundary for the high efficacy interval.
#' @param target.toxicity The target DLT rate.
#' @param target.efficacy The target efficacy rate.
#' @param ncohort The total number of cohorts.
#' @param cohortsize The number of patients in the cohort.
#' @param n.early The early stopping parameter. If the number of patients treated at
#'                the current dose reaches \code{n.early}, then we stop the trial
#'                and select the optimal biological dose (OBD) based on the observed data. The default value is 100.
#' @param startdose The starting dose level.
#' @param p.true A vector containing the true toxicity probabilities of the
#'              investigational dose levels.
#' @param q.true A vector containing the true efficacy probabilities of the
#'              investigational dose levels.
#' @param ntrial The total number of trials to be simulated.
#' @param seed The random seed for simulation.
#' @param p1 The cutoff lower limit for safety utility function 1, described in the Details section.
#' @param p2 The cutoff upper limit for safety utility function 1, described in the Details section.
#' @param q1 The cutoff lower limit for efficacy utility function 1, described in the Details section.
#' @param q2 The cutoff upper limit for efficacy utility function 1, described in the Details section.
#' @param cutoff.eli.toxicity The cutoff value to eliminate a dose with unacceptable high toxicity for safety. The default value is 0.95.
#' @param cutoff.eli.efficacy The cutoff value for the futility rule, the acceptable lowest efficacy. The default value is 0.30.
#' @param w1.toxicity The weight for toxicity utility function 2 and 3,described in the Details section.
#' @param w2.toxicity The weight for toxicity utility function 3, described in the Details section.
#' @param indicator The indicator cutoff value for utility function 3, described in the Details section.
#'
#' @details Large trials are simulated to characterize the operating characteristics of the Keyboard design under the prespecified true toxicity probabilities and true efficacy probabilities of the investigational doses. The dose assignment follows the rules described in the function \code{get.decision.obd.kb()}.
#'
#' The following stopping rules are built in the Keyboard design:
#' \enumerate{
#' \item Stop the trial if the lowest dose is eliminated from the trial due to high unacceptable toxicity.
#' \item Stop the trial if the number of patients treated at the current dose exceeds \code{n.earlystop}.
#' }
#'
#'
#'
#'
#' @return \code{get.oc.obd.kb()} returns the operating characteristics of the Keyboard design as a list, including:\cr
#' \enumerate{
#'    \item the selection percentage at each dose level using utility function 1 (\code{$selpercent1}), \cr
#'    \item the selection percentage at each dose level using utility function 2 (\code{$selpercent2}), \cr
#'    \item the selection percentage at each dose level using utility function 3 (\code{$selpercent3}), \cr
#'    \item the number of patients treated at each dose level (\code{$npatients}), \cr
#'    \item the number of dose-limiting toxicities (DLTs) observed at each dose level (\code{$ntox}), \cr
#'    \item the number of responses observed at each dose level (\code{$neff}), \cr
#'    \item the average number of DLTs (\code{$totaltox}), \cr
#'    \item the average number of responses (\code{$totaleff}), \cr
#'    \item the average number of patients (\code{$totaln}), \cr
#'    \item the percentage of early stopping without selecting the OBD using utility function 1 (\code{$percentstop1}), \cr
#'    \item the percentage of early stopping without selecting the OBD using utility function 2 (\code{$percentstop2}), \cr
#'    \item the percentage of early stopping without selecting the OBD using utility function 3 (\code{$percentstop3}), \cr
#'    \item data.frame (\code{$simu.setup}) containing simulation parameters, such as target, p.true, etc.
#'
#'
#' }
#' 
#' @author Hongying Sun, Li Tang, and Haitao Pan
#' @examples
#' \donttest{
#' toxicity.low <- 0.15
#' toxicity.moderate <- 0.25
#' toxicity.high <- 0.35
#' efficacy.low <- 0.25
#' efficacy.moderate <- 0.45
#' efficacy.high <- 0.65
#' target.toxicity <- 0.30
#' target.efficacy <- 0.40
#' p.true <-c(0.08,0.30,0.60,0.80)
#' q.true <- c(0.25,0.40,0.25,0.50)
#' oc.obd.kb <- get.oc.obd.kb(toxicity.low = toxicity.low,
#'              toxicity.moderate= toxicity.moderate,
#'              toxicity.high = toxicity.high,
#'              efficacy.low = efficacy.low,
#'              efficacy.moderate = efficacy.moderate,
#'              efficacy.high = efficacy.high,
#'              target.toxicity=target.toxicity,
#'              target.efficacy= target.efficacy,
#'              p.true= p.true, q.true= q.true)
#' oc.obd.kb
#' summary.kb(oc.obd.kb)
#' plot.kb(oc.obd.kb)
#' plot.kb(oc.obd.kb$selpercent1)
#' plot.kb(oc.obd.kb$selpercent2)
#' plot.kb(oc.obd.kb$selpercent3)
#' plot.kb(oc.obd.kb$npatients)
#' plot.kb(oc.obd.kb$ntox)
#' plot.kb(oc.obd.kb$neff)
#' }
#' @family single-agent phase I/II functions
#'
#' @references
#' Li DH, Whitmore JB, Guo W, Ji Y.  Toxicity and efficacy probability interval design for phase I adoptive cell therapy dose-finding clinical trials.
#' \emph{Clinical Cancer Research}. 2017; 23:13-20.
#'https://clincancerres.aacrjournals.org/content/23/1/13.long
#'
#' 
#' Liu S, Johnson VE.  A robust Bayesian dose-finding design for phase I/II clinical trials. \emph{Biostatistics}. 2016; 17(2):249-63.
#' https://academic.oup.com/biostatistics/article/17/2/249/1744018
#'
#' Zhou Y, Lee JJ, Yuan Y.  A utility-based Bayesian optimal interval (U-BOIN) phase I/II design to identify the optimal biological dose for targeted and immune therapies. \emph{Statistics in Medicine}. 2019; 38:S5299-5316.
#' https://onlinelibrary.wiley.com/doi/epdf/10.1002/sim.8361
#' @export
get.oc.obd.kb <- function( toxicity.low, toxicity.moderate,toxicity.high, efficacy.low, efficacy.moderate,
                           efficacy.high,target.toxicity, target.efficacy,ncohort=10, cohortsize=3, n.early=100,
                           startdose=1, p.true, q.true, ntrial = 1000, seed = 6, p1=0.15, p2=0.40, q1=0.3, q2=0.6,cutoff.eli.toxicity= 0.95,
                           cutoff.eli.efficacy=0.3, w1.toxicity =0.33, w2.toxicity=1.09, indicator =  target.toxicity){



    set.seed(seed)
    ndose = length(p.true);
    npts = ncohort * cohortsize;
    #toxicity outcome matrix
    Y=matrix(rep(0, ndose*ntrial), ncol=ndose);
    #efficacy outcome matrix
    E=matrix(rep(0, ndose*ntrial), ncol=ndose);
    # matrix to store the total number of patients
    N=matrix(rep(0, ndose*ntrial), ncol=ndose);
    # matrix to store selected dose level
    # dselect = rep(0, ntrial)
    dselect1 = rep(0, ntrial)
    dselect2 = rep(0, ntrial)
    dselect3 = rep(0, ntrial)
    dearlystop = rep(0, ntrial)

    # get decision table
    decision.matrix.output <- get.decision.obd.kb(toxicity.low = toxicity.low, toxicity.moderate= toxicity.moderate, toxicity.high = toxicity.high,
                                           efficacy.low = efficacy.low, efficacy.moderate = efficacy.moderate, efficacy.high = efficacy.high, target.toxicity=target.toxicity, target.efficacy=target.efficacy, cohortsize=cohortsize, ncohort=ncohort)$decision.matrix

    # this function outputs the decision given the dose-finding table, the number of total patients, the number of patients who experienced toxicity and the number of responses
    decision.finding <- function(out.matrix, n, t, r){
        rowindex <- which(out.matrix$N==n & out.matrix$T==t & out.matrix$R == r)
        decision <- out.matrix$Decision[rowindex]
        decision <- as.character(decision)
        return (decision)
    }


    # simulation trials
    for (trial in 1:ntrial){
        # the number of patients who experienced toxity at each level
        y <- rep(0, ndose);
        # the number of patients who experienced reponse at each level
        e <- rep(0, ndose);
        # the number of total patients treated at each level
        n <- rep(0, ndose);
        # an indicator to check whehter the trial terminates early

        d=startdose;
        earlystop = 0;
        # an indicator to check whehter the doses are eliminated
        elimi = rep(0, ndose);

        for (i in 1:ncohort){
            y[d] = y[d] + sum(runif(cohortsize) < p.true[d]);
            e[d] = e[d] + sum(runif(cohortsize) < q.true[d]);
            n[d] = n[d] + cohortsize;
            if(n[d] >= n.early) break;
            #get dose decision matrix
            decision.result = decision.finding(decision.matrix.output, n[d], y[d], e[d])
            # early stop rules
            if (n[d] >=3 ){
                if (1- pbeta(target.toxicity, y[d]+1, n[d]-y[d]+1) > cutoff.eli.toxicity){
                    elimi[d:ndose]=1;
                    if( d==1) {
                        earlystop=1;
                        break;}
                }
            }
            if (n[d] >=3 ){
                if (1- pbeta(target.efficacy, e[d]+1, n[d]-e[d]+1) < cutoff.eli.efficacy){
                    elimi[d]=1;
                }
            }




            if(!is.null(decision.result)){

                # dose transition
                if(decision.result == "E" && d != ndose) {
                    if(elimi[d+1]==0) {d=d+1;}
                } else if(decision.result == "D" && d != 1) {
                    if(elimi[d-1]==0) {d=d-1;}
                } else if (decision.result == "S") {
                    d = d ;
                } else if (decision.result == "DUT" && d != 1){
                    if(elimi[d-1]==0) {d=d-1;}
                    elimi[d:ndose]=1;
                } else if (decision.result=="DUE" && d != 1){
                    if(elimi[d-1]==0) {d=d-1;}
                    elimi[d]=1;
                } else if (decision.result == "EUE" && d != ndose){
                    if(elimi[d+1]==0) {d=d+1;}
                    elimi[d]=1;
                } else { d=d; }
            }
        }
        Y[trial,] = y;
        E[trial,] = e;
        N[trial,] = n;
        if(earlystop ==1){
            dselect1[trial]=99;
            dselect2[trial]=99;
            dselect3[trial]=99;
            dearlystop[trial]=199;

        } else{
            dselect1[trial]= select.obd.kb( target.toxicity =target.toxicity, target.efficacy=target.efficacy, npts = n, ntox = y, neff = e)$obd1
            dselect2[trial]= select.obd.kb( target.toxicity =target.toxicity, target.efficacy=target.efficacy, npts = n, ntox = y, neff = e)$obd2
            dselect3[trial]= select.obd.kb( target.toxicity =target.toxicity, target.efficacy=target.efficacy, npts = n, ntox = y, neff = e)$obd3
        } # use select.obd function.
    }

    # output results
    selpercent1 = rep(0, ndose);
    selpercent2 = rep(0, ndose);
    selpercent3 = rep(0, ndose);
    nptsdose = apply(N,2, mean);
    ntoxdose = apply(Y, 2, mean);
    neffdose = apply(E,2, mean);
    for (i in 1:ndose){
        selpercent1[i] = sum(dselect1==i)/ntrial*100;
        selpercent2[i] = sum(dselect2==i)/ntrial*100;
        selpercent3[i] = sum(dselect3==i)/ntrial*100;
    }
    out=list(name = "get.oc.obd.kb",  ## to identify object for summary.kb() function.
        selpercent1=selpercent1, selpercent2=selpercent2,selpercent3=selpercent3,npatients=nptsdose, ntox=ntoxdose, neff=neffdose, totaltox=sum(Y)/ntrial,totaleff=sum(E)/ntrial,
             totaln=sum(N)/ntrial, earlystop = sum(dearlystop==199)/ntrial*100,percentstop1=sum(dselect1== 99)/ntrial*100,percentstop2=sum(dselect2== 99)/ntrial*100,percentstop3=sum(dselect3== 99)/ntrial*100,
             simu.setup=data.frame(target.toxicity=target.toxicity, target.efficacy=target.efficacy,p.true=p.true, q.true =q.true,
                                   ncohort=ncohort, cohortsize = cohortsize,
                                   startdose = startdose,
                                   ntrial = ntrial, dose=1:ndose)
    );
    return (out)
}












