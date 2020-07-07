#' Select the Optimal Biological Dose (OBD) for Single-agent Phase I/II Trials
#'
#' This function selects the optimal biological dose (OBD) at the end of a single-agent phase I/II trial.
#' @param target.toxicity The target dose-limiting toxicity (DLT) rate.
#' @param target.efficacy  The target efficacy rate.
#' @param npts The vector containing the total number of patients treated at each dose level.
#' @param ntox The vector containing the number of subjects who experienced toxicities at each dose level.
#' @param neff The vector containing the number of subjects who experienced efficacies at each dose level.
#' @param p1 The cutoff lower limit for safety utility function 1, described in the Details section.
#' @param p2 The cutoff upper limit for safety utility function 1, described in the Details section.
#' @param q1 The cutoff lower limit for efficacy utility function 1, described in the Details section.
#' @param q2 The cutoff upper limit for efficacy utility function 1, described in the Details section.
#' @param cutoff.eli.toxicity The cutoff value to eliminate a dose with unacceptable high toxicity for safety.
#'                            The default value is 0.95.
#' @param cutoff.eli.efficacy The cutoff value for futility rule, the acceptable lowest efficacy.
#'                            The default value is 0.3.
#' @param w1.toxicity The weight for toxicity utility function 2 and 3, described in the Details section.
#' @param w2.toxicity The weight for toxicity utility function 3, described in the Details section.
#' @param indicator The indicator cutoff value for utility function 3, described in the Details section.
#' @details \code{select.obd.kb()} selects the OBD that is the most desirable based on benefit-risk tradeoff considering both toxicity and efficacy outcomes.  A utility score is used to quantify the desirability of all admissible doses. Calculation of utility scores requires the posterior probabilities for toxicity \eqn{p_i} and efficacy \eqn{q_i}, which can be computed by using \eqn{beta(\alpha_p + x_i,\beta_p + n_i - x_i)} and \eqn{beta(\alpha_q + y_i,\beta_q + n_i - y_i) } assuming that the prior for both \eqn{p_i} and \eqn{q_i} follows independent beta distributions \eqn{beta(\alpha_p,\beta_p)} and \eqn{beta(\alpha_q,\beta_q)}. Three criteria are used to calculate the desirability in this function.
#'
#' The first criterion relies on a utility function for toxicity \eqn{f_1(p)}, where p denotes the toxicity rate, and on a utility function for efficacy \eqn{f_2(q)}, where q denotes the efficacy rate. \eqn{f_1(p)} is 1 if \eqn{p \in [0, p1)}; \eqn{f_1(p)} is 0 if \eqn{p \in [p2, 1]}; \eqn{f_1(p)} is \eqn{1- (p-p1)/(p2-p1)} if \eqn{p \in [p1, p2)}. \eqn{f_2(p)} is 1 if \eqn{p \in (0, p1)}. Here, p1 is the cutoff lower limit and p2 is the cutoff upper limit for safety utility function 1\eqn{f_1(p)}.
#'
#' Similarly, \eqn{f_2(q)} is 1 if \eqn{p \in [0, q1)}; \eqn{f_2(q)} is 0 if \eqn{p \in [q2, 1]}; \eqn{f_2(q)} is \eqn{1- (p-q1)/(q2-q1)} if \eqn{p \in [q1, q2)}. \eqn{f_2(p)} is 1 if \eqn{p \in (0, q1)}. Here, q1 is the cutoff lower limit and q2 is the cutoff upper limit for safety utility function \eqn{f_2(q)}.
#'
#' The utility score that quantifies benefit-risk tradeoff at the current dose i is calculated as follows:
#' \deqn{ U(p_i, q_i) =f_1(p) * f_2(q) }
#'
#'
#' The second criterion depends on a marginal toxicity probability \eqn{\pi_{T,i} = beta(\alpha_p + x_i,\beta_p + n_i - x_i) } and a marginal efficacy probability \eqn{\pi_{E,i} = beta(\alpha_q + y_i,\beta_q + n_i - y_i) }.  Then the utility score is calculated as follows:
#' \deqn{U_i= \pi_{E,i} - w_1*\pi_{T,i} }
#'
#' The third criterion also depends on a marginal toxicity probability \eqn{\pi_{T,i} = beta(\alpha_p + x_i,\beta_p + n_i - x_i) } and a marginal efficacy probability \eqn{\pi_{E,i} = beta(\alpha_q + y_i,\beta_q + n_i - y_i) }, but it has an additional penalty when the posterior toxicity probability is high by using an indicator function.  Then utility score using this function is calculated as follows:
#'
#' \deqn{U_j= \pi_{E,i} - w_1*\pi_{T,i}-w_2*\pi_{T,i}*I(\pi_{T,i}>\rho)}
#' Here, the recommended \eqn{\rho} is the target toxicity rate.
#'
#' Once the utility score is computed for all the doses, the optimal biological dose is calculated as follows:
#' \deqn{ d = argmax_i[ U(p_i, q_i) | D]}
#'
#' @return \code{select.obd.kb()} returns the selected dose: \cr
#' \enumerate{
#' \item Selected OBD level using utility function 1 (\code{$obd1}), as described in the Details section. \cr
#' \item Selected OBD level using utility function 2 (\code{$obd2}), as described in the Details section.  \cr
#' \item Selected OBD level using utility function 3 (\code{$obd3}), as described in the Details section.  
#' }
#' @author Hongying Sun, Li Tang, and Haitao Pan
#' @examples
#' \donttest{
#' target.toxicity<-0.3
#' target.efficacy<-0.4
#' npts <- c(3,6,12,3,3)
#' ntox <-  c(1,2,4,2,3)
#' neff <-  c(0,0,5,1,1)
#' obd <- select.obd.kb (target.toxicity=target.toxicity,
#'        target.efficacy= target.efficacy, npts = npts,
#'        ntox = ntox, neff =  neff)
#' print(obd)
#' }
#'
#' @family single-agent phase I/II functions
#'
#' @references
#'  Li DH, Whitmore JB, Guo W, Ji Y.  Toxicity and efficacy probability interval design for phase I adoptive cell therapy dose-finding clinical trials.
#' \emph{Clinical Cancer Research}. 2017; 23:13-20.
#'https://clincancerres.aacrjournals.org/content/23/1/13.long 
#' 
#' Liu S, Johnson VE.  A robust Bayesian dose-finding design for phase I/II clinical trials. \emph{Biostatistics}. 2016; 17(2):249-63.
#' https://academic.oup.com/biostatistics/article/17/2/249/1744018
#'
#' Zhou Y, Lee JJ, Yuan Y. A utility-based Bayesian optimal interval (U-BOIN) phase I/II design to identify the optimal biological dose for targeted and immune therapies.
#' \emph{Statistics in Medicine}. 2019; 38:S5299-5316.
#' https://onlinelibrary.wiley.com/doi/epdf/10.1002/sim.8361
#' @export
select.obd.kb <- function (target.toxicity, target.efficacy, npts , ntox , neff , p1=0.15, p2=0.40, q1=0.3, q2=0.6,
                           cutoff.eli.toxicity= 0.95, cutoff.eli.efficacy=0.3,
                           w1.toxicity =0.33, w2.toxicity=1.09, indicator =target.toxicity){
    ## isotonic transformation using the pool adjacent violator algorithm (PAVA)
    pava <- function(x, wt = rep(1, length(x))) {
        n <- length(x)
        if (n <= 1)
            return(x)
        if (any(is.na(x)) || any(is.na(wt))) {
            stop("Missing values in 'x' or 'wt' not allowed")
        }
        lvlsets <- (1:n)
        repeat {
            viol <- (as.vector(diff(x)) < 0)
            if (!(any(viol)))
                break
            i <- min((1:(n - 1))[viol])
            lvl1 <- lvlsets[i]
            lvl2 <- lvlsets[i + 1]
            ilvl <- (lvlsets == lvl1 | lvlsets == lvl2)
            x[ilvl] <- sum(x[ilvl] * wt[ilvl])/sum(wt[ilvl])
            lvlsets[ilvl] <- lvl1
        }
        x
    }

    # check whether the dose has been elimited during the trial
    y=ntox;
    n=npts;
    e=neff;
    ndose = length(npts);
    elimi = rep(0,ndose);
    for (i in 1:ndose){
        if (n[i] >=3 ){
            if (1- pbeta(target.toxicity, y[i]+1, n[i]-y[i]+1) > cutoff.eli.toxicity){
                elimi[i:ndose]=1;
                break;
            }
        }
    }
    for (i in 1:ndose){
        if (n[i] >=3 ){
            if ( 1-pbeta(target.efficacy, e[i]+1, n[i]-e[i]+1) <  cutoff.eli.efficacy){
                elimi[i]=1;

            }
        }
    }

    # no dose should be selected  all doses are never used to treat patients.
    if(sum(n[elimi==0])==0){
        selectdose1=99;
        selectdose2=99;
        selectdose3=99;
    }
    adm.set =(n!=0) & (elimi==0);
    if(sum(adm.set)==0) {
      selectdose1=99;
      selectdose2=99;
      selectdose3=99;
      } else{
        adm.index = which(adm.set == T);
        n.adm = n[adm.set];
        y.adm = y[adm.set];
        e.adm = e[adm.set];
        # posterior mean and variance of toxicity probabilities using beta(0.05,0.05) as the prior
        prior.alpha = 0.05;
        prior.beta = 0.05;
        posterior.alpha = prior.alpha + y.adm;
        posterior.beta = prior.beta + n.adm - y.adm;
        phat = posterior.alpha/(posterior.alpha + posterior.beta);
        phat.var = posterior.alpha * posterior.beta/((posterior.alpha + posterior.beta)^2*(posterior.alpha + posterior.beta+1));

        # perform the isotonic transformation using PAVA
        phat = pava(phat, wt = 1/phat.var);
        phat = phat + (1:length(phat)) * 1e-10  ## break ties by adding an increasingly small number

        # posterior mean of efficacy probabilities using beta(0.05,0.05)
        q.posterior.alpha = prior.alpha + e.adm;
        q.posterior.beta = prior.beta + n.adm - e.adm;
        qhat = q.posterior.alpha/(q.posterior.alpha + q.posterior.beta);


        # define f(p) and f(q) function
        f.p = rep(0, length(phat));
        f.q = rep(0,length(qhat));
        for (i in 1:length(phat)){
            # print (i)
            # print(phat)
            # print (phat[i])
            if( !is.na( phat[i]) && (phat[i] <= p1) && (phat[i] >= 0 )){
                f.p[i] = 1;
            } else if ( !is.na( phat[i]) &&(phat[i]< p2) && (phat[i] > p1) ){
                f.p[i] = 1- (phat[i]-p1)/(p2-p1);
            } else if (!is.na( phat[i]) && (phat[i] <= 1) && (phat[i] >= p2) ){
                f.p[i] =0;
            }


        }
        for (i in 1:length(qhat)) {
            if ( !is.na( qhat[i]) &&(qhat[i] < q1) &&  (qhat[i] >=0 )){
                f.q[i] = 0;
            } else if ( !is.na( qhat[i]) &&(qhat[i] >=q1) && (qhat[i] < q2) ){
                f.q[i] = (qhat[i] - q1)/(q2-q1);
            } else if (!is.na( qhat[i]) && (qhat[i] <= 1) && ( qhat[i]>= q2)) {
                f.q[i] = 1;
            }
        }

        utility1 = f.p*f.q;
        if(sum(utility1) == 0){
          selectdose1 =99;
        } else {
          selectd1 = sort(-utility1, index.return = T)$ix[1];
          selectdose1 = adm.index[selectd1];
        }
        utility2 = qhat - w1.toxicity*phat;
        if(sum(utility2) == 0){
          selectdose2 =99;
        } else {
          selectd2 = sort(-utility2, index.return = T)$ix[1];
          selectdose2 = adm.index[selectd2];
        }

        utility3 = qhat - w1.toxicity*phat - w2.toxicity*phat*I(phat>indicator);
        utility3 = as.numeric(utility3);
        if(sum(utility3) == 0){
          selectdose3 =99;
        } else {
          selectd3 = sort(-utility3, index.return = T)$ix[1];
          selectdose3 = adm.index[selectd3];
        }











    }



    # out = list(obd = selectdose);
    out = list(name = "select.obd.kb",  ## to identify object for summary.kb() function.
               obd1=selectdose1, obd2=selectdose2, obd3=selectdose3);

    return(out)
}
