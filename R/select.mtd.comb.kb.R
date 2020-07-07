#' Maximum Tolerated Dose (MTD) Selection for Drug-combination Trials
#'
#' This function selects the maximum tolerated dose (MTD) after the drug-combination trial is
#' completed.
#'
#' @details
#' The Keyboard design starts by specifying a target toxicity interval
#' (referred to as the "target key") such that any dose with a toxicity
#' probability within that interval can be practically viewed as the MTD.
#' Based on this interval's width, the Keyboard design forms a series of
#' equally wide keys that span the rest of the range from 0 to 1.
#'
#' Given all the observed data once the trial is completed, this function
#' selects the MTD based on matrix isotonic estimates of the toxicity
#' probability of each dose combination, selecting the dose whose estimate is
#' closest to the target. When there are ties, one is randomly chosen.
#' These (matrix) isotonic estimates are obtained with
#' \code{\link[Iso]{biviso}}.
#'
#' For patient safety, the following dose elimination rule is evaluated after
#' each cohort:
#' if at least 3 patients have been treated at the given dose and the
#' observed data indicate that there is more than a 95\% chance that the current
#' dose is above the MTD, then we eliminate this dose and beyond from the trial to
#' prevent exposing future patients to these overly toxic doses. The
#' probability threshold for elimination can be specified with
#' \code{cutoff.eli}. When a dose is eliminated, the design recommends the next
#' lower dose for treating the next patient. If the lowest dose combination
#' (1, 1) is overly toxic, then the trial terminates early and no dose is selected
#' as the MTD.
#'
#' For some applications, investigators may prefer a stricter stopping rule
#' for extra safety when the lowest dose is overly toxic. This can be achieved
#' by setting \code{extrasafe=TRUE}, which imposes the following stricter
#' safety stopping rule:\cr
#' Stop the trial if \cr
#' (i) the number of patients treated at the lowest dose \eqn{\ge 3}, and \cr
#' (ii) \deqn{Pr[(toxicity rate of the lowest dose > target) | data)]
#'            > cutoff.eli - offset}
#' As a tradeoff, the strong stopping rule will decrease the MTD selection
#' percentage when the lowest dose is the true MTD.
#'
#' @param target The target dose-limiting toxicity (DLT) rate.
#' @param npts A \code{J*K} matrix \code{(J<=K)} containing the number of
#'             patients treated at each dose combination.
#' @param ntox A \code{J*K} matrix \code{(J<=K)} containing the number of
#'             patients who experienced a dose-limiting toxicity at each dose
#'             combination.
#' @param cutoff.eli The cutoff value to eliminate an overly toxic dose and all
#'                   higher doses for safety.\cr
#'                   The default value is 0.95.
#' @param extrasafe Set \code{extrasafe=TRUE} to impose a stricter
#'                  stopping rule.\cr
#'                  The default is FALSE.
#' @param offset A small positive number (between 0 and 0.5) to control how
#'               strict the stopping rule is when \code{extrasafe=TRUE}. A
#'               larger value leads to a stricter stopping rule.\cr
#'               The default value is 0.05.
#'
#' @return The function returns a list with: \cr
#' \enumerate{
#'   \item the target toxicity probability (\code{$target}),\cr
#'   \item the selected MTD (\code{$MTD}),\cr
#'   \item a matrix with the isotonic estimates of the DLT probability at each
#'       dose (\code{$p_est}).
#' }
#'
#' @note The MTD selection and dose escalation/de-escalation rules are two
#'   independent components of the trial design. When appropriate, another dose
#'   selection procedure (e.g., one based on a fitted logistic model) can be used
#'   to select the MTD after completing the trial using the Keyboard design.
#' @author Hongying Sun, Li Tang, and Haitao Pan
#' @examples
#' ### Drug-combination trial ###
#'
#' ## Select the MTD based on the data from a 3 x 5 combination trial
#' n <- matrix(c(3, 5, 0, 0, 0,
#'               7, 6, 15, 0, 0,
#'               0, 0, 4, 0, 0), ncol=5, byrow=TRUE)
#' y <- matrix(c(0, 1, 0, 0, 0,
#'               1, 1, 4, 0, 0,
#'               0, 0, 2, 0, 0), ncol=5, byrow=TRUE)
#'
#' sel.comb <- select.mtd.comb.kb(target=0.3, npts=n, ntox=y)
#'
#' summary.kb(sel.comb)
#' plot.kb(sel.comb)
#'
#' @section Uses:
#' This function uses \code{\link[Iso]{biviso}}.
#'
#' @family drug-combination functions
#'
#' @references
#'
#' @references
#'
#' Yan F, Mandrekar SJ, Yuan Y. Keyboard: A Novel Bayesian Toxicity Probability
#' Interval Design for Phase I Clinical Trials.
#' \emph{Clinical Cancer Research}. 2017; 23:3994-4003.
#' http://clincancerres.aacrjournals.org/content/23/15/3994.full-text.pdf
#' 
#' Pan H, Lin R, Yuan Y. Keyboard design for phase I drug-combination trials.
#' \emph{Contemporary Clinical Trials}. 2020.
#' https://doi.org/10.1016/j.cct.2020.105972
#'
#' @export
select.mtd.comb.kb <- function(target, npts, ntox, cutoff.eli = 0.95,
                               extrasafe = FALSE, offset = 0.05) {

    if (!requireNamespace("Iso", quietly = TRUE)) {
        stop("Package \"Iso\" needed for this function to work.",
             "Please install it.", call. = FALSE)
    }

    y = ntox; n = npts;
    if (nrow(n) > ncol(n) | nrow(y) > ncol(y)) {
        warning("Error: npts and ntox should be arranged in a way (i.e., rotated)",
            "such that for each of them, the number of rows is less than or",
            "equal to the number of columns.")
        return()
    }

    elimi = matrix(0, dim(n)[1], dim(n)[2])

    if (extrasafe) {
        if (n[1,1] >= 3) {
            if (1 - pbeta(target, y[1,1] + 1, n[1,1] - y[1,1] + 1) > cutoff.eli - offset) {
                elimi[,] = 1
            }
        }
    }

    for (i in 1:dim(n)[1]) {
        for (j in 1:dim(n)[2]) {
            if (n[i,j] >= 3) {
                if (1 - pbeta(target, y[i,j] + 1, n[i,j] - y[i,j] + 1) > cutoff.eli) {
                    elimi[i:dim(n)[1], j] = 1
                    elimi[i, j:dim(n)[2]] = 1
                    break
                }
            }
        }
    }

    if (elimi[1] == 1) { ## no dose should be selected if the first dose is already very toxic
        selectdose = c(99, 99)
        selectdoses = matrix(selectdose, nrow = 1)
    }
    else {
        phat = (y + 0.05) / (n + 0.1)
        ## perform the isotonic transformation using PAVA
        phat = Iso::biviso(phat, n + 0.1, warn = TRUE)[,]
        phat.out = phat
        phat.out[n == 0] = NA
        phat[elimi == 1] = 1.1 # to aviod selecting eliminated dose
        ## break the ties
        phat = phat * (n != 0) + (1E-5) * (matrix(rep(1:dim(n)[1],
                                                      each = dim(n)[2],
                                                      len = length(n)),
                                                  dim(n)[1],byrow=T) +
                                           matrix(rep(1:dim(n)[2],
                                                      each = dim(n)[1],
                                                      len = length(n)),
                                                  dim(n)[1]))
        ## select dose closest to the target as the MTD
        phat[n == 0] = 10 ## so that the dose without treating patients will not be selected

        selectdose = which(abs(phat - target) == min(abs(phat - target)), arr.ind = TRUE)
        if (length(selectdose) > 2) { ##if there are still ties, randomly pick the first one.
            selectdose = selectdose[1,]
        }

        ## mtd.contour == TRUE will activate the option of multiple MTDs selection [COMMENTED OUT BELOW]
        aa = function(x) as.numeric(as.character(x))
        ## previously "if (mtd.contour == FALSE)":
        selectdoses = matrix(99, nrow=1, ncol=2)
        selectdoses[1,] = matrix(selectdose, nrow = 1)

        #else {
        #    selectdoses = cbind('row' = 1:dim(n)[1], 'col' = rep(99,dim(n)[1]))
        #    for (k in dim(n)[1]:1) {
        #        kn = n[k,]
        #        ky = y[k,]
        #        kelimi = elimi[k,]
        #        kphat = phat[k,]
        #
        #        if (kelimi[1] == 1 || sum(n[kelimi == 0]) == 0) {
        #            kseldose=99
        #        }
        #        else {
        #            adm.set = (kn != 0) & (kelimi == 0)
        #            adm.index = which(adm.set == T)
        #            y.adm = ky[adm.set]
        #            n.adm = kn[adm.set]
        #            ## select dose closest to the target as the MTD
        #            selectd = sort(abs(kphat[adm.set] - target), index.return = T)$ix[1]
        #            kseldose = adm.index[selectd]
        #        }
        #        selectdoses[k, 2] = ifelse(is.na(kseldose), 99, kseldose)
        #        if (k < dim(n)[1]) {
        #            if (selectdoses[k + 1, 2] == dim(n)[2]) {
        #                selectdoses[k,2] = dim(n)[2]
        #            }
        #        }
        #        #if(k<dim(n)[1]) if(aa(selectdoses[k+1,2])==aa(selectdoses[k,2])) selectdoses[k,2] = 99
        #        if (k < dim(n)[1]) {
        #            if (aa(selectdoses[k + 1, 2]) == dim(n)[2] &
        #            aa(selectdoses[k + 1, 2]) == aa(selectdoses[k,2])) {
        #                selectdoses[k,2] = 99
        #            }
        #        }
        #    }
        #}
        selectdoses = matrix(selectdoses[selectdoses[,2] != 99,], ncol = 2)
        colnames(selectdoses) = c('DoseA', 'DoseB')
    }


    ## previously "if (mtd.contour == FALSE)":
    if (selectdoses[1,1] == 99 && selectdoses[1,2] == 99) {
        warning("All tested doses are overly toxic. No MTD is selected! \n")
        return(list(target = target, MTD = 99, p_est = matrix(NA,
                                                              nrow = dim(npts)[1],
                                                              ncol = dim(npts)[2])))
    }
    else {
        return(list(target = target, MTD = selectdoses, p_est = round(phat.out,2)))
    }

    ## For mtd.contour == TRUE
    #else {
    #    if (length(selectdoses) == 0) {
    #        warning("All tested doses are overly toxic. No MTD is selected! \n")
    #        return(list(target = target, MTD = 99, p_est = matrix(NA,
    #                                                              nrow = dim(npts)[1],
    #                                                              ncol = dim(npts)[2])))
    #    }
    #    else {
    #        return(list(target = target, MTD = selectdoses, p_est = round(phat.out,2)))
    #    }
    #}
}
