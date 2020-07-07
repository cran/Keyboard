#' Illustration of the Keyboard Design (Single-agent)
#'
#' This function serves solely as an example to illustrate how the
#' Keyboard design determines dosing decisions when conducting a single-agent
#' trial.
#'
#' @details
#' It is straightforward to visualize how the Keyboard design determines whether to
#' stay at the current dose, escalate, or de-escalate. The "target key" (shown in blue) is
#' the proper dosing interval determined for the trial by the investigators.
#' The "strongest key" (in red) is the interval with the highest posterior
#' probability (given the dosage data observed thus far) of the current dose's
#' true DLT rate. If the strongest key is to the left of the
#' target key, then we escalate the dose; if the strongest key is to the
#' right  of the target key, then we de-escalate the dose; and if the strongest
#' key is the target key, then we retain the current dose. Thus, posterior
#' probabilities are what directs dosage. Graphically, the strongest key is the
#' one with the largest area under the posterior distribution curve of the DLT
#' rate of the current dose.
#'
#' @param center The center of the target key (between 0 and 1).
#' @param half_width The width, or tolerable deviation, to each side of the
#'                   center (half of the target key's width).
#' @param s1 The lower boundary of the strongest key.
#' @param s2 The upper boundary of the strongest key.
#' @param a The alpha parameter of the beta distribution.
#' @param b The beta parameter of the beta distribution.
#'
#' @return A dose decision plot with the posterior distribution of the DLT rate
#' of the current dose and the positions of the target and the strongest keys.
#' @author Hongying Sun, Li Tang, and Haitao Pan
#' @examples
#' ## Clear all plots before switching between graphical parameters
#' opar <- par(no.readonly = TRUE)
#' on.exit(par(opar))
#' par(mfrow = c(3, 1)) # for many plots in the same screen
#' par(mar = c(5, 5, 2, 2)) # for only one plot per page
#'
#' example.kb(center = 0.19, half_width = 0.03, s1 = 0.4, s2 = 0.46, a = 3, b = 4) # de-escalation
#' example.kb(center = 0.19, half_width = 0.03, s1 = 0.04, s2 = 0.1, a = 2, b = 12) # escalation
#' example.kb(center = 0.19, half_width = 0.03, s1 = 0.16, s2 = 0.22, a = 2, b = 5) # stay
#' @references
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
example.kb <- function(center, half_width, s1, s2, a, b) {
    ## Divide the range between 0 and 1 into keys of the same width, based on
    ## the target key's width.
    ## @param center: Of the target key.
    ## @param half_width: The width to each side of the center, (half of the
    ##                    target key's width.)
    getkeys <- function(center, half_width) {
        c1 = center - half_width
        c2 = center + half_width

        delta = c2 - c1
        lkey = NULL
        rkey = NULL

        i = 1
        cutoff = c1 - (i * delta)
        while (cutoff > 0) {
            lkey = c(cutoff, lkey)
            i = i + 1
            cutoff = c1 - (i * delta)
        }
        lkey[lkey < 0] = 0

        i = 1
        cutoff = c2 + (i * delta)
        while(cutoff < 1) {
            rkey = c(rkey, cutoff)
            i = i + 1
            cutoff = c2 + (i * delta)
        }
        rkey[rkey > 1] = 1

        keys = c(lkey, c1, c2, rkey)
        return(keys)
    }

    ## Shade a given key on its beta distribution plot.
    ## @param center: of the key to be shaded.
    ## @param half_width: the width to each side of the center, or half of the key's width.
    ## @param a: alpha parameter of the beta distribution.
    ## @param b: beta parameter of the beta distribution.
    ## @param color: of the key to be shaded.
    shadekey <- function(center, half_width, a, b, color) {
        c1 = center - half_width
        c2 = center + half_width

        cord.x <- c(c1, seq(c1, c2, 0.01), c2)
        cord.y <- c(0, dbeta(seq(c1, c2, 0.01), a, b),0)
        polygon(cord.x, cord.y, col = color)
    }

    keys = getkeys(center = center, half_width = half_width)

    #.old <- par(no.readonly = TRUE)
    #par(mar = c(5, 5, 2, 2))
    plot(seq(0, 1, 0.01), dbeta(seq(0, 1, 0.01), a, b),
         type = "l", xlab = "DLT rate", ylab = "Density", cex.lab = 1.4)
    for(i in 1:length(keys)) {
        lines(c(keys[i], keys[i]), c(0, dbeta(keys[i], a, b)))
    }
    if (center + half_width == s2) {
        text(0.9, 2, "Stay", cex=1.8)
        shadekey(center, half_width, a, b, "skyblue")  #target key
        shadekey(center, half_width, a, b, "red")      #strongest key
        s_interval = sprintf("(%.2f, %.2f)", s1, s2)
        text(s2, (dbeta(s1, a, b) + dbeta(s2, a, b))/2, s_interval, cex = 1.0, pos = 4, offset = 2)
    }
    else if (center + half_width < s2) {
        text(0.9, 1.7, "De-escalate", cex = 1.8)
        shadekey(center, half_width, a, b, "skyblue")
        shadekey((s1 + s2)/2, s2 - (s1 + s2)/2, a, b, "red")
        t_interval = sprintf("(%.2f, %.2f)", center - half_width, center + half_width)
        text(center, (dbeta(center - half_width, a, b) + dbeta(center + half_width, a, b))/2,
             t_interval, cex = 1.0, pos = 2)
    }
    else {
        text(0.9, 4, "Escalate", cex=1.8)
        shadekey(center, half_width, a, b, "skyblue")
        shadekey((s1 + s2)/2, s2 - (s1 + s2)/2, a, b, "red")
        t_interval = sprintf("(%.2f, %.2f)", center - half_width, center + half_width)
        text(center, (dbeta(center - half_width, a, b) + dbeta(center + half_width, a, b))/2,
             t_interval, cex = 1.0, pos = 4)
    }
    #par(.old)
}



