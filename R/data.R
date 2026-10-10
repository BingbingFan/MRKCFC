#' Mayo Clinic primary biliary cirrhosis follow-up data
#'
#' Longitudinal follow-up of 312 patients from the Mayo Clinic trial of
#' D-penicillamine in primary biliary cirrhosis, in long format with one row
#' per visit. Included as a worked example of the sparse, irregularly observed
#' multivariate longitudinal design this package is written for: visits are
#' unequally spaced and differ in number between patients.
#'
#' \code{serBilir} and \code{albumin} are the two markers used in the
#' \code{\link{MKCFC}} example.
#'
#' @format A data frame of 1945 rows and 20 columns, 312 subjects:
#' \describe{
#'   \item{id}{subject identifier, a factor with 312 levels.}
#'   \item{years}{follow-up time to death, transplantation or censoring.}
#'   \item{status}{outcome at the end of follow-up: alive, transplanted or
#'     dead.}
#'   \item{drug}{treatment arm, placebo or D-penicillamine.}
#'   \item{age}{age at entry, in years.}
#'   \item{sex}{male or female.}
#'   \item{year}{time of this visit, in years from entry. This is the
#'     observation time for the longitudinal markers.}
#'   \item{ascites, hepatomegaly, spiders, edema}{clinical signs recorded at
#'     the visit.}
#'   \item{serBilir}{serum bilirubin, mg/dl.}
#'   \item{serChol}{serum cholesterol, mg/dl.}
#'   \item{albumin}{serum albumin, g/dl.}
#'   \item{alkaline}{alkaline phosphatase, U/l.}
#'   \item{SGOT}{aspartate aminotransferase, U/ml.}
#'   \item{platelets}{platelet count.}
#'   \item{prothrombin}{prothrombin time, seconds.}
#'   \item{histologic}{histologic stage of disease.}
#'   \item{status2}{indicator of death, 1 for dead and 0 otherwise.}
#' }
#'
#' @source Distributed with the \pkg{JM} package; originally from the Mayo
#'   Clinic trial reported by Murtaugh et al. (1994).
#'
#' @references
#' Murtaugh, P. A., Dickson, E. R., Van Dam, G. M., Malinchoc, M., Grambsch,
#' P. M., Langworthy, A. L. and Gips, C. H. (1994). Primary biliary cirrhosis:
#' prediction of short-term survival based on repeated patient visits.
#' \emph{Hepatology}, 20, 126--134.
#'
#' @examples
#' data(pbc2)
#' str(pbc2[, c("id", "year", "serBilir", "albumin")])
#' # visits per subject: sparse and irregular
#' summary(as.numeric(table(pbc2$id)))
"pbc2"
