#' Maximum Tolerated Dose (MTD) Selection for Single-agent Trials
#'
#' This function selects the maximum tolerated dose (MTD) after the single-agent trial is
#' completed.
#'
#' @details
#' The Keyboard design starts by specifying a target toxicity interval
#' (referred to as the "target key") such that any dose with a toxicity
#' probability within that interval can be practically viewed as the MTD.
#' Based on this interval's width, the Keyboard design forms a series of
#' equally wide keys that span the rest of range from 0 to 1.
#'
#' This function selects the MTD based on isotonic estimates of toxicity
#' probabilities, selecting dose level \eqn{j*} for which the isotonic estimate
#' of the DLT rate is closest to the target. If there are ties, then we select from
#' the ties the highest dose level when the estimate of the DLT rate is smaller
#' than the target, or the lowest dose level when the estimate of the DLT rate
#' is greater than the target. The isotonic estimates are obtained by applying the
#' pooled-adjacent-violators algorithm (PAVA) [Barlow, 1972].
#'
#'
#' For some applications, investigators may prefer a stricter stopping rule
#' to ensure the lowest dose is not overly toxic. This can be achieved
#' by setting \code{extrasafe=TRUE}, which imposes the following stricter
#' safety stopping rule:\cr
#' Stop the trial if \cr
#' (i) the number of patients treated at the lowest dose \eqn{\ge 3}, and \cr
#' (ii) \deqn{Pr((toxicity rate of the lowest dose > target) | data)
#'            > cutoff.eli - offset}
#' As a tradeoff, the strong stopping rule will decrease the MTD selection
#' percentage when the lowest dose actually is the MTD.
#'
#' @param target The target dose-limiting toxicity (DLT) rate.
#' @param npts A vector containing the number of patients treated at each dose level.
#' @param ntox A vector containing the number of patients at each dose level who experienced a DLT at each dose level.
#' @param cutoff.eli The cutoff to eliminate an overly toxic dose and all
#'                   higher doses for safety.\cr
#'                   The default value is 0.95.
#' @param extrasafe Set \code{extrasafe=TRUE} to impose a stricter
#'                  stopping rule.\cr
#'                  The default is FALSE.
#' @param offset A small positive number (between 0 and 0.5) to control how
#'               strict the stopping rule is when \code{extrasafe=TRUE}. A
#'               larger value leads to a stricter stopping rule.\cr
#'              The default value is 0.05.
#'
#' @return The function returns a list with: \cr
#' \enumerate{
#'   \item the target toxicity probability (\code{$target}),\cr
#'   \item the selected MTD (\code{$MTD}),\cr
#'   \item the isotonic estimates of the DLT probability at each dose and corresponding \code{95\%} credible interval (\code{$p_est}),\cr
#'   \item the probability of overdosing defined as\cr
#'       \eqn{Pr(toxicity > target | data)} (\code{$p_overdose}).
#' }
#'
#' @note The MTD selection and dose escalation/de-escalation rules are two
#'   independent components of the trial design. When appropriate, another dose
#'   selection procedure (e.g., one based on a fitted logistic model) can be used
#'   to select the MTD after completing the trial using the Keyboard design.
#' @author Hongying Sun, Li Tang, and Haitao Pan
#' @examples
#' ### Single-agent trial ###
#'
#' n <- c(3, 3, 15, 9, 0)
#' y <- c(0, 0, 4, 4, 0)
#'
#' selmtd <- select.mtd.kb(target=0.3, npts=n, ntox=y)
#' 
#' selmtd
#'
#' @family single-agent functions
#'
#' @references
#'
#' Yan F, Mandrekar SJ, Yuan Y. Keyboard: A Novel Bayesian Toxicity Probability
#' Interval Design for Phase I Clinical Trials.
#' \emph{Clinical Cancer Research}. 2017; 23:3994-4003.
#' http://clincancerres.aacrjournals.org/content/23/15/3994.full-text.pdf
#' @export
select.mtd.kb <- function(target, npts, ntox, cutoff.eli = 0.95,
                          extrasafe = FALSE, offset = 0.05) {
    ## isotonic transformation using the pool adjacent violator algorithm (PAVA)
    pava <- function(x, wt = rep(1, length(x))) {
        n <- length(x)
        if (n <= 1) {
            return(x)
        }
        if (any(is.na(x)) || any(is.na(wt))) {
            stop("Missing values in 'x' or 'wt' not allowed")
        }
        lvlsets <- (1:n)
        repeat {
            viol <- (as.vector(diff(x)) < 0)
            if (!(any(viol))) {
                break
            }
            i <- min((1:(n - 1))[viol])
            lvl1 <- lvlsets[i]
            lvl2 <- lvlsets[i + 1]
            ilvl <- (lvlsets == lvl1 | lvlsets == lvl2)
            x[ilvl] <- sum(x[ilvl] * wt[ilvl])/sum(wt[ilvl])
            lvlsets[ilvl] <- lvl1
        }
        x # ?? printing x ??
    }
    ## determine whether the dose has been eliminated during the trial
    y = ntox
    n = npts
    ndose = length(n)
    elimi = rep(0, ndose)
    for (i in 1:ndose) {
        if (n[i] >= 3) {
            if (1 - pbeta(target, y[i] + 1, n[i] - y[i] + 1) > cutoff.eli) {
                elimi[i:ndose] = 1
                break
            }
        }
    }

    if (extrasafe) {
        if (n[1] >= 3) {
            if (1 - pbeta(target, y[1] + 1, n[1] - y[1] + 1) > cutoff.eli - offset) {
                elimi[1:ndose] = 1
            }
        }
    }

    ## no dose should be selected (i.e., selectdose=99) if the first dose is already very toxic
    ## or all uneliminated doses are never used to treat patients
    if (elimi[1] == 1 || sum(n[elimi == 0]) == 0) {
        selectdose = 99
    }
    else {
        adm.set = (n != 0) & (elimi == 0)
        adm.index = which(adm.set == T)
        y.adm = y[adm.set]
        n.adm = n[adm.set]

        ## poster mean and variance of toxicity probabilities using beta(0.05, 0.05) as the prior
        phat = (y.adm + 0.05)/(n.adm + 0.1)
        phat.var = (y.adm + 0.05) * (n.adm - y.adm + 0.05)/((n.adm + 0.1)^2 * (n.adm + 0.1 + 1))

        ## perform the isotonic transformation using PAVA
        phat = pava(phat, wt = 1/phat.var)
        phat = phat + (1:length(phat)) * 1e-10  ## break ties by adding an increasingly small number
        selectd = sort(abs(phat - target), index.return = T)$ix[1]  ## select dose closest to the target as the MTD
        selectdose = adm.index[selectd]
    }

    # old: if (verbose == TRUE) {remaining of function}:
    trtd = (n != 0)
    poverdose = pava(1 - pbeta(target, y[trtd] + 0.05, n[trtd] - y[trtd] + 0.05))
    phat.all = pava((y[trtd] + 0.05)/(n[trtd] + 0.1),
                    wt = 1/((y[trtd] + 0.05) * (n[trtd] - y[trtd] + 0.05)/((n[trtd] + 0.1)^2 * (n[trtd] + 0.1 + 1))))

    A1 = A2 = A3 = A4 = NULL
    k=1
    ## output summary statistics
    for (i in 1:ndose) {
        if (n[i] > 0) {
            A1 = append(A1, formatC(phat.all[k], digits = 2, format = "f"))
            A2 = append(A2, formatC(qbeta(0.025, y[i] + 0.05, n[i] - y[i] + 0.05), digits = 2, format = "f"))
            A3 = append(A3, formatC(qbeta(0.975, y[i] + 0.05, n[i] - y[i] + 0.05), digits = 2, format = "f"))
            A4 = append(A4, formatC(poverdose[k], digits = 2, format = "f"))
            k = k+1
        }
        else {
            # no estimate output for doses never used to treat patients
            A1 = append(A1, "----")
            A2 = append(A2, "----")
            A3 = append(A3, "----")
            A4 = append(A4, "----")
        }
    }
    p_est = data.frame(cbind('dose'=1:length(npts), 'phat'=A1, 'CI'=paste("(", A2,",", A3,")",sep="")))

    out = list(target = target, MTD = selectdose, p_est = p_est, p_overdose = A4)

    return(out)
}
