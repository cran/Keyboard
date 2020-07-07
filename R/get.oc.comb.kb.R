#' Operating Characteristics for Drug-combination Trials
#'
#' This function generates the operating characteristics of the Keyboard design for
#' drug-combination trials.
#'

#' @param target The target dose-limiting toxicity (DLT) rate.
#' @param p.true A \code{J*K} matrix \code{(J<=K)} containing the true toxicity
#'               probabilities of combinations with \code{J} dose levels of
#'               agent A and \code{K} dose levels of agent B.
#' @param ncohort A scalar specifying the total number of cohorts in the trial.
#' @param cohortsize The number of patients in the cohort.
#' @param n.earlystop The early stopping parameter. If the number of patients
#'                    treated at the current dose reaches \code{n.earlystop},
#'                    then we stop the trial and select the MTD based on the observed
#'                    data.\cr
#'                    The default value is 100.
#' @param marginL The difference between the target and the lower bound of the
#'                "target key" (proper dosing interval) to be defined.\cr
#'                The default is 0.05.
#' @param marginR The difference between the target and the upper bound of the
#'                "target key" (proper dosing interval) to be defined.\cr
#'                The default is 0.05.
#' @param startdose The starting dose combination level for the
#'                  drug-combination trial.\cr
#'                  The default is c(1, 1).
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
#' @param ntrial The total number of trials to be simulated. \cr
#'               The default value is 1000.
#'
#' @return The function returns the operating characteristics of the Keyboard
#'   combination design as a list: \cr
#'  \enumerate{
#'   \item the true toxicity probability at each dose level (\code{$p.true}),\cr
#'   \item the selection percentage at each dose level (\code{$selpercent}),\cr
#'   \item the percentage of correct selection (\code{$pcs}),\cr
#'   \item the number of patients treated at each dose level (\code{$nptsdose}),\cr
#'   \item the number of toxicities observed at each dose level (\code{$ntoxdose}),\cr
#'   \item the total number of toxicities observed in the trial (\code{$totaltox}),\cr
#'   \item the total number of patients in the trial (\code{$totaln}),\cr
#'   \item the total percentage of patients treated at the MTD (\code{$npercent}).
#' }
#' @author Hongying Sun, Li Tang, and Haitao Pan
#' @examples
#' ### Drug-combination trial ###
#'
#' p.true <- matrix(c(0.01, 0.03, 0.10, 0.20, 0.30,
#'                    0.03, 0.05, 0.15, 0.30, 0.60,
#'                    0.08, 0.10, 0.30, 0.60, 0.75), byrow=TRUE, ncol=5)
#'
#' oc.comb <- get.oc.comb.kb(target=0.3, p.true=p.true, ncohort=20, cohortsize=3,
#'                           n.earlystop=12, startdose=c(1, 1), ntrial=100)
#'
#' @section Uses:
#' This function uses \code{\link{get.boundary.comb.kb}} 
#'
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
#' @export
get.oc.comb.kb <- function(target, p.true, ncohort, cohortsize,
                           n.earlystop = 100, marginL = 0.05, marginR = 0.05,
                           startdose = c(1, 1), cutoff.eli = 0.95,
                           extrasafe = FALSE, offset = 0.05, ntrial = 1000) {

    # if (!requireNamespace("BOIN", quietly = TRUE)) {
    #     stop("Package \"BOIN\" needed for this function to work.",
    #          "Please install it.", call. = FALSE)
    # }

    JJ = nrow(p.true)
    KK = ncol(p.true)
    if (JJ > KK) {
        stop("p.true should be arranged in a way such that the number of rows",
             " is less than or equal to the number of columns (i.e. rotated).")
    }
    #if (JJ > KK) {
    #    p.true = t(p.true)
    #}
    if (target < 0.05) {
        stop("The target is too low.")
    }
    if (target > 0.6) {
        stop("The target is too high.")
    }
    if (offset >= 0.5) {
        stop("The offset is too large.")
    }
    if (n.earlystop <= 6) {
        stop("The value of n.earlystop is too low to ensure good operating ",
             "characteristics.\n  ",
             "n.earlystop = 9 to 18 is recommended.")
    }

    set.seed(6)
    ndose = length(p.true)
    npts = ncohort * cohortsize
    Y <- array(matrix(rep(0, length(p.true) * ntrial),
                      dim(p.true)[1]), dim = c(dim(p.true), ntrial))
    N <- array(matrix(rep(0, length(p.true) * ntrial),
                      dim(p.true)[1]), dim = c(dim(p.true), ntrial))
    dselect = matrix(rep(0, 2 * ntrial), ncol = 2)

    total.incoherent.dlt = rep(0,ntrial)        
    total.incoherent.no.dlt = rep(0,ntrial)
    total.incoherent = rep(0,ntrial)

    total.long.incoherent.dlt = rep(0,ntrial)       
    total.long.incoherent.no.dlt = rep(0,ntrial)
    total.long.incoherent = rep(0,ntrial)

  
    temp = get.boundary.comb.kb(target, ncohort, cohortsize, cutoff.eli=0.95)$boundary

    b.e = temp[2, ]
    b.d = temp[3, ]
    b.elim = temp[4, ]

    for (trial in 1:ntrial) {
        y <- matrix(rep(0, ndose), dim(p.true)[1], dim(p.true)[2])
        n <- matrix(rep(0, ndose), dim(p.true)[1], dim(p.true)[2])

        incoherent.dlt = incoherent.no.dlt = 0   
        long.incoherent.dlt = long.incoherent.no.dlt = 0   
        earlystop = 0
        d = startdose
        elimi = matrix(rep(0, ndose), dim(p.true)[1], dim(p.true)[2])
        for (pp in 1:ncohort) {
            y.current <- y[d[1], d[2]]    

            y[d[1], d[2]] = y[d[1], d[2]] + sum(runif(cohortsize) < p.true[d[1], d[2]])
            n[d[1], d[2]] = n[d[1], d[2]] + cohortsize
            if (n[d[1], d[2]] >= n.earlystop) {
                break
            }
            nc = n[d[1], d[2]]
            if (!is.na(b.elim[nc])) {
                if (y[d[1], d[2]] >= b.elim[nc]) {
                    for (i in min(d[1], dim(p.true)[1]):dim(p.true)[1]) {
                        for (j in min(d[2], dim(p.true)[2]):dim(p.true)[2]) {
                            elimi[i, j] = 1
                        }
                    }
                    if (d[1] == 1 && d[2] == 1) {
                        d = c(99, 99)
                        earlystop = 1
                        break
                    }
                }
                if (extrasafe) {
                    if (d[1] == 1 && d[2] == 1 && n[1, 1] >= 3) {
                        if (1 - pbeta(target, y[1, 1] + 1, n[1, 1] - y[1, 1] + 1) >
                            cutoff.eli - offset) {
                            d = c(99, 99)
                            earlystop = 1
                            break
                        }
                    }
                }
            }
            if (y[d[1], d[2]] <= b.e[nc]) {
                elevel = matrix(c(1, 0, 0, 1), 2)
                pr_H0 = rep(0, length(elevel) / 2)
                nn = pr_H0
                for (i in seq(1, length(elevel) / 2, by = 1)) {
                    if (d[1] + elevel[1, i] <= dim(p.true)[1] &&
                        d[2] + elevel[2, i] <= dim(p.true)[2]) {
                        if (elimi[d[1] + elevel[1, i], d[2] + elevel[2, i]] == 0) {
                            yn = y[d[1] + elevel[1, i], d[2] + elevel[2, i]]
                            nn[i] = n[d[1] + elevel[1, i], d[2] + elevel[2, i]]
                            pr_H0[i] <- pbeta(target+marginR, yn + 0.5, nn[i] - yn + 0.5) -
                                        pbeta(target-marginL, yn + 0.5, nn[i] - yn + 0.5)
                        }
                    }
                }
                pr_H0 = pr_H0 + nn * 5e-04
                if (max(pr_H0) == 0) {
                    d = d
                }
                else {
                    ##count incoherent: if DLT, still escalate;
                    if ((y[d[1], d[2]] - y.current) > 0) {
                        incoherent.dlt <- incoherent.dlt + 1    #Sep 27,2017
                    }
                    if (y[d[1], d[2]] / n[d[1], d[2]] > target + marginR) {
                        long.incoherent.dlt <- long.incoherent.dlt+1
                    }

                    k = which(pr_H0 == max(pr_H0))[as.integer(runif(1) * length(which(pr_H0 == max(pr_H0))) + 1)]
                    d = d + c(elevel[1, k], elevel[2, k])
                }
            }
            else if (y[d[1], d[2]] >= b.d[nc]) {
                delevel = matrix(c(-1, 0, 0, -1), 2)
                pr_H0 = rep(0, length(delevel) / 2)
                nn = pr_H0
                for (i in seq(1, length(delevel) / 2, by = 1)) {
                    if (d[1] + delevel[1, i] > 0 && d[2] + delevel[2, i] > 0) {
                        yn = y[d[1] + delevel[1, i], d[2] + delevel[2, i]]
                        nn[i] = n[d[1] + delevel[1, i], d[2] + delevel[2, i]]
                        pr_H0[i] = pbeta(target + marginR, yn + 0.5, nn[i] - yn + 0.5) -
                                   pbeta(target - marginL, yn + 0.5, nn[i] - yn + 0.5)
                    }
                }
                pr_H0 = pr_H0 + nn * 5e-04
                if (max(pr_H0) == 0) {
                    d = d
                }
                else {
                    ##count incoherent: if no-DLT, still de-escalate;
                    if ((y[d[1], d[2]] - y.current) == 0) {
                        incoherent.no.dlt <- incoherent.no.dlt + 1    #Sep 27,2017
                    }
                    if (y[d[1], d[2]] / n[d[1], d[2]] < target - marginL) {
                        long.incoherent.no.dlt <- long.incoherent.no.dlt + 1
                    }

                    k = which(pr_H0 == max(pr_H0))[as.integer(runif(1) * length(which(pr_H0 == max(pr_H0))) + 1)]
                    d = d + c(delevel[1, k], delevel[2, k])
                }
            }
            else {
                d = d
            }
        }
        Y[, , trial] = y
        N[, , trial] = n

        total.incoherent.dlt[trial] = incoherent.dlt        #Sep 27,2017
        total.incoherent.no.dlt[trial] = incoherent.no.dlt
        total.incoherent[trial] = incoherent.dlt + incoherent.no.dlt

        total.long.incoherent.dlt[trial] = long.incoherent.dlt        #Sep 27,2017
        total.long.incoherent.no.dlt[trial] = long.incoherent.no.dlt
        total.long.incoherent[trial] = long.incoherent.dlt + long.incoherent.no.dlt

        if (earlystop == 1) {
            dselect[trial, ] = c(99, 99)
        }
        else {
            # select.mtd.comb is adapted from BOIN 
 select.mtd.comb <- function (target, npts, ntox, cutoff.eli = 0.95, extrasafe = FALSE, offset = 0.05, mtd.contour = FALSE) {
    y = ntox
    n = npts
    if (nrow(n) > ncol(n) | nrow(y) > ncol(y)) {
        stop("npts and ntox should be arranged in a way (i.e., rotated) such that for each of them, the number of rows is less than or equal to the number of columns.")
    }
    elimi = matrix(0, dim(n)[1], dim(n)[2])
    if (extrasafe) {
        if (n[1, 1] >= 3) {
            if (1 - pbeta(target, y[1, 1] + 1, n[1, 1] - y[1, 
                1] + 1) > cutoff.eli - offset) {
                elimi[, ] = 1
            }
        }
    }
    for (i in 1:dim(n)[1]) {
        for (j in 1:dim(n)[2]) {
            if (n[i, j] >= 3) {
                if (1 - pbeta(target, y[i, j] + 1, n[i, j] - 
                  y[i, j] + 1) > cutoff.eli) {
                  elimi[i:dim(n)[1], j] = 1
                  elimi[i, j:dim(n)[2]] = 1
                  break
                }
            }
        }
    }
    if (elimi[1] == 1) {
        selectdose = c(99, 99)
        selectdoses = matrix(selectdose, nrow = 1)
    }
    else {
        phat = (y + 0.05)/(n + 0.1)
        phat = Iso::biviso(phat, n + 0.1, warn = TRUE)[, ]
        phat.out = phat
        phat.out[n == 0] = NA
        phat[elimi == 1] = 1.1
        phat = phat * (n != 0) + (1e-05) * (matrix(rep(1:dim(n)[1], 
            each = dim(n)[2], len = length(n)), dim(n)[1], byrow = T) + 
            matrix(rep(1:dim(n)[2], each = dim(n)[1], len = length(n)), 
                dim(n)[1]))
        phat[n == 0] = 10
        selectdose = which(abs(phat - target) == min(abs(phat - 
            target)), arr.ind = TRUE)
        if (length(selectdose) > 2) 
            selectdose = selectdose[1, ]
        aa = function(x) as.numeric(as.character(x))
        if (mtd.contour == TRUE) {
            selectdoses = cbind(row = 1:dim(n)[1], col = rep(99, 
                dim(n)[1]))
            for (k in dim(n)[1]:1) {
                kn = n[k, ]
                ky = y[k, ]
                kelimi = elimi[k, ]
                kphat = phat[k, ]
                if (kelimi[1] == 1 || sum(n[kelimi == 0]) == 
                  0) {
                  kseldose = 99
                }
                else {
                  adm.set = (kn != 0) & (kelimi == 0)
                  adm.index = which(adm.set == T)
                  y.adm = ky[adm.set]
                  n.adm = kn[adm.set]
                  selectd = sort(abs(kphat[adm.set] - target), 
                    index.return = T)$ix[1]
                  kseldose = adm.index[selectd]
                }
                selectdoses[k, 2] = ifelse(is.na(kseldose), 99, 
                  kseldose)
                if (k < dim(n)[1]) 
                  if (selectdoses[k + 1, 2] == dim(n)[2]) 
                    selectdoses[k, 2] = dim(n)[2]
                if (k < dim(n)[1]) 
                  if (aa(selectdoses[k + 1, 2]) == dim(n)[2] & 
                    aa(selectdoses[k + 1, 2]) == aa(selectdoses[k, 
                      2])) 
                    selectdoses[k, 2] = 99
            }
        }
        else {
            selectdoses = matrix(99, nrow = 1, ncol = 2)
            selectdoses[1, ] = matrix(selectdose, nrow = 1)
        }
        selectdoses = matrix(selectdoses[selectdoses[, 2] != 
            99, ], ncol = 2)
        colnames(selectdoses) = c("DoseA", "DoseB")
    }
    if (mtd.contour == FALSE) {
        if (selectdoses[1, 1] == 99 && selectdoses[1, 2] == 99) {
            message("All tested doses are overly toxic. No MTD is selected! \n")
            out = list(target = target, MTD = 99, p_est = matrix(NA, 
                nrow = dim(npts)[1], ncol = dim(npts)[2]))
        }
        else {
            out = list(target = target, MTD = selectdoses, p_est = round(phat.out, 
                2))
        }
        return(out)
    }
    else {
        if (length(selectdoses) == 0) {
            message("All tested doses are overly toxic. No MTD is selected! \n")
            out = list(target = target, MTD = 99, p_est = matrix(NA, 
                nrow = dim(npts)[1], ncol = dim(npts)[2]))
        }
        else {
            out = list(target = target, MTD = selectdoses, p_est = round(phat.out, 
                2))
        }
        return(out)
    }
}

            
            selcomb = select.mtd.comb(target, n, y, cutoff.eli, extrasafe, offset)
            dselect[trial, 1] = selcomb$MTD[1]
            dselect[trial, 2] = selcomb$MTD[2]
        }
    }
    selpercent = matrix(rep(0, ndose), dim(p.true)[1], dim(p.true)[2])
    nptsdose = apply(N, c(1, 2), mean, digits = 2, format = "f")
    ntoxdose = apply(Y, c(1, 2), mean, digits = 2, format = "f")
    for (i in 1:dim(p.true)[1]) {
        for (j in 1:dim(p.true)[2]) {
            selpercent[i, j] = sum(dselect[, 1] == i & dselect[, 2] == j) / ntrial * 100
        }
    }
    if (JJ <= KK) {
        # message("True toxicity rate of dose combinations:\n")
        # for (i in 1:dim(p.true)[2]) {
        #   message(formatC(p.true[, i], digits = 2, format = "f",
        #               width = 5), sep = "  ", "\n")
        # }
        # message("\n")
        # message("selection percentage at each dose combination (%):\n")
        # for (i in 1:dim(p.true)[2]) {
        #   message(formatC(selpercent[, i], digits = 2, format = "f",
        #               width = 5), sep = "  ", "\n")
        # }
        # message("\n")
        # message("number of patients treated at each dose combination:\n")
        # for (i in 1:dim(p.true)[2]) {
        #   message(formatC(apply(N, c(1, 2), mean)[, i], digits = 2,
        #               format = "f", width = 5), sep = "  ", "\n")
        # }
        # message("\n")
        # message("number of toxicity observed at each dose combination:\n")
        # for (i in 1:dim(p.true)[2]) {
        #   message(formatC(apply(Y, c(1, 2), mean)[, i], digits = 2,
        #               format = "f", width = 5), sep = "  ", "\n")
        # }
        # message("\n")
        # message("average number of toxicities:", formatC(t(sum(Y))/ntrial,
        #                                              digits = 1, format = "f"), "\n")
        # message("average number of patients:", formatC(t(sum(N))/ntrial,
        #                                            digits = 1, format = "f"), "\n")
        # message("selection percentage of MTD:", formatC(sum(selpercent[which(p.true ==
        #                                                                    target, arr.ind = TRUE)]), digits = 1, format = "f"),
        #     "\n")
        # message("percentage of patients treated at MTD:", formatC(sum(nptsdose[which(p.true ==
        #                                                                            target, arr.ind = TRUE)])/npts * 100, digits = 1,
        #                                                       format = "f"), "\n")
        # message("percentage of early stopping due to toxicity:",
        #     formatC(100 - sum(selpercent), digits = 2, format = "f"),
        #
        #     "\n")
        # ### Pan edit;
        # message("percentage of patients treated above MTD:", formatC(sum(nptsdose[which(p.true >
        #                                                                               target, arr.ind = TRUE)])/npts * 100, digits = 1,
        #                                                          format = "f"), "\n")
        # message("percentage of patients treated below MTD:", formatC(sum(nptsdose[which(p.true <
        #                                                                               target, arr.ind = TRUE)])/npts * 100, digits = 1,
        #                                                          format = "f"), "\n")

        over_dosing_60_num <- NULL
        for(i in 1:ntrial){
            over_dosing_60_num[i] <- sum(N[,,i][p.true > target + marginR])
        }
        over_dosing_60 = mean(over_dosing_60_num > 0.6 * npts) * 100

        over_dosing_80_num <- NULL
        for(i in 1:ntrial){
            over_dosing_80_num[i] <- sum(N[,,i][p.true > target + marginR])
        }
        over_dosing_80 = mean(over_dosing_80_num > 0.8 * npts) * 100

        under_dosing_60_num <- NULL
        for(i in 1:ntrial){
            under_dosing_60_num[i] <- sum(N[,,i][p.true < target - marginR])
        }
        under_dosing_60 = mean(under_dosing_60_num > 0.6 * npts) * 100

        under_dosing_80_num <- NULL
        for(i in 1:ntrial) {
            under_dosing_80_num[i] <- sum(N[,,i][p.true < target - marginR])
        }
        under_dosing_80 = mean(under_dosing_80_num > 0.8 * npts) * 100
        ###
    }
    else {
        # message("True toxicity rate of dose combinations:\n")
        # for (i in 1:dim(p.true)[2]) {
        #   message(formatC(p.true[, i], digits = 2, format = "f",
        #               width = 5), sep = "  ", "\n")
        # }
        # message("\n")
        # message("selection percentage at each dose combination (%):\n")
        # for (i in 1:dim(p.true)[2]) {
        #   message(formatC(selpercent[, i], digits = 2, format = "f",
        #               width = 5), sep = "  ", "\n")
        # }
        # message("\n")
        # message("number of patients treated at each dose combination:\n")
        # for (i in 1:dim(p.true)[2]) {
        #   message(formatC(apply(N, c(1, 2), mean)[, i], digits = 2,
        #               format = "f", width = 5), sep = "  ", "\n")
        # }
        # message("\n")
        # message("number of toxicity observed at each dose combination:\n")
        # for (i in 1:dim(p.true)[2]) {
        #   message(formatC(apply(Y, c(1, 2), mean)[, i], digits = 2,
        #               format = "f", width = 5), sep = "  ", "\n")
        # }
        # message("\n")
        # message("average number of toxicities:", formatC(t(sum(Y))/ntrial,
        #                                              digits = 1, format = "f"), "\n")
        # message("average number of patients:", formatC(t(sum(N))/ntrial,
        #                                            digits = 1, format = "f"), "\n")
        # message("selection percentage of MTD:", formatC(sum(selpercent[which(p.true ==
        #                                                                    target, arr.ind = TRUE)]), digits = 1, format = "f"),
        #     "\n")
        # message("percentage of patients treated at MTD:", formatC(sum(nptsdose[which(p.true ==
        #                                                                            target, arr.ind = TRUE)])/npts * 100, digits = 1,
        #                                                       format = "f"), "\n")
        # message("percentage of early stopping due to toxicity:",
        #     formatC(100 - sum(selpercent), digits = 2, format = "f"),
        #
        #     "\n")
        # ### Pan edit;
        # message("percentage of patients treated above MTD:", formatC(sum(nptsdose[which(p.true >
        #                                                                               target, arr.ind = TRUE)])/npts * 100, digits = 1,
        #                                                          format = "f"), "\n")
        # message("percentage of patients treated below MTD:", formatC(sum(nptsdose[which(p.true <
        #                                                                               target, arr.ind = TRUE)])/npts * 100, digits = 1,
        #                                                          format = "f"), "\n")
        #
        over_dosing_60_num <- NULL
        for(i in 1:ntrial) {
            over_dosing_60_num[i] <- sum(N[,,i][p.true > target + marginR])
        }
        over_dosing_60 = mean(over_dosing_60_num > 0.6 * npts) * 100

        over_dosing_80_num <- NULL
        for(i in 1:ntrial){
            over_dosing_80_num[i] <- sum(N[,,i][p.true > target + marginR])
        }
        over_dosing_80 = mean(over_dosing_80_num > 0.8 * npts) * 100

        under_dosing_60_num <- NULL
        for(i in 1:ntrial){
            under_dosing_60_num[i] <- sum(N[,,i][p.true < target - marginR])
        }
        under_dosing_60 = mean(under_dosing_60_num > 0.6 * npts) * 100

        under_dosing_80_num <- NULL
        for(i in 1:ntrial){
            under_dosing_80_num[i] <- sum(N[,,i][p.true < target - marginR])
        }
        under_dosing_80 = mean(under_dosing_80_num > 0.8 * npts) * 100
        ###
    }

    out = list(name = "get.oc.comb.kb",  ## to identify object for summary.kb() function.
               target = target,
               p.true = p.true,
               ncohort = ncohort,
               cohortsize = cohortsize,
               ntotal = ncohort * cohortsize,
               ntrial = ntrial,
               selpercent = selpercent,
               pcs = sum(selpercent[which(p.true <= target + marginR & p.true >= target - marginR, arr.ind = TRUE)]),
               # nmtd = sum(nptsdose[which(p.true <= target+marginR&p.true>=target-marginR, arr.ind = TRUE)])/npts * 100,
               nmtd = sum(nptsdose[p.true <= target + marginR & p.true >= target - marginR]) / sum(nptsdose) * 100,

               nmtd.over = sum(nptsdose[p.true > target + marginR]) / sum(nptsdose) * 100,
               nmtd.over.sd = sd( nptsdose[p.true > target + marginR] / sum(nptsdose)),
               # nmtd_ex = sum(nptsdose[p.true>target]),

               nptsdose = nptsdose,
               ntoxdose = ntoxdose,
               totaltox = round(sum(Y) / ntrial, 1), ## added 28/06/19
               totaln = round(sum(N) / ntrial, 1), ## added 28/06/19

               npercent = paste(round(sum(nptsdose[which(abs(p.true - target) == min(abs(p.true - target)), arr.ind = TRUE)])
                                      / sum(nptsdose) * 100, 1), '%', sep = ""),  ## added 28/06/19

               # low_dose=sum(nptsdose[which(p.true < target-marginR, arr.ind = TRUE)])/npts * 100,
               # high_tox = sum(nptsdose[which(p.true > target+marginR, arr.ind = TRUE)])/npts * 100,
               low_dose = sum(nptsdose[p.true < target - marginR]) / sum(nptsdose) * 100,
               high_tox = sum(nptsdose[p.true > target + marginR]) / sum(nptsdose) * 100,

               over_dosing_60 = over_dosing_60,
               over_dosing_80 = over_dosing_80,
               under_dosing_60 = under_dosing_60,
               under_dosing_80 = under_dosing_80,

               pctearlystop = sum(dselect == 99) / ntrial * 100,
               over.sel.pcs = sum(selpercent[which(p.true >= target + marginR, arr.ind = TRUE)]),
               percent.under.MTD.sel.20 = sum(selpercent[which(p.true < target - marginR, arr.ind = TRUE)]),

               total.incoherent.dlt = sum(total.incoherent.dlt) / (ncohort * ntrial) * 100, #Sep 27,2017
               total.incoherent.no.dlt = sum(total.incoherent.no.dlt) / (ncohort * ntrial) * 100,
               total.incoherent = sum(total.incoherent) / (ncohort * ntrial) * 100,

               total.long.incoherent.dlt = sum(long.incoherent.dlt) / (ncohort * ntrial) * 100,
               total.long.incoherent.no.dlt = sum(long.incoherent.no.dlt) / (ncohort * ntrial) * 100,
               total.long.incoherent = sum(total.long.incoherent) * 100
               )
    return(out)
}
