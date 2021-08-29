#' @include utils.R
NULL

#' Pipe operator
#'
#' This package uses the pipe operator (`\%>\%`) to express nested code
#' as a series of imperative procedures.
#'
#' @param lhs,rhs An object and a function.
#' @seealso [magrittr::%>%()], [tee()].
#' @return An object.
#' @examples
#' # set seed for reproducibility
#' set.seed(500)
#'
#' # generate 100 random numbers and calculate the mean
#' mean(runif(100))
#'
#' # reset the seed
#' set.seed(500)
#'
#' # repeat the previous procedure but use the pipe operator instead of nesting
#' # function calls inside each other.
#' runif(100) %>% mean()
#'
#' @name %>%
#' @rdname pipe
#' @aliases pipe
#' @importFrom magrittr %>%
#' @export
NULL

#' Central colour repository
#' @name ibis_colours
#' @examples
#' ibis_colours[['viridis_plasma']]
#' @keywords internal
#' @noRd
ibis_colours <- list(
  sdm_colour = colorRampPalette(c('grey90','steelblue4','steelblue1','gold','red1','red4'))(100),
  viridis_orig = c("#440154FF","#482878FF","#3E4A89FF","#31688EFF","#26828EFF","#1F9E89FF","#35B779FF","#6DCD59FF","#B4DE2CFF","#FDE725FF"),
  viridis_cividis = c("#00204DFF","#00336FFF","#39486BFF","#575C6DFF","#707173FF","#8A8779FF","#A69D75FF","#C4B56CFF","#E4CF5BFF","#FFEA46FF"),
  viridis_plasma = c("#0D0887FF","#47039FFF","#7301A8FF","#9C179EFF","#BD3786FF","#D8576BFF","#ED7953FF","#FA9E3BFF","#FDC926FF","#F0F921FF"),
  distinct_random =  c("#000000", "#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941", "#006FA6", "#A30059",
                "#FFDBE5", "#7A4900", "#0000A6", "#63FFAC", "#B79762", "#004D43", "#8FB0FF", "#997D87",
                "#5A0007", "#809693", "#FEFFE6", "#1B4400", "#4FC601", "#3B5DFF", "#4A3B53", "#FF2F80",
                "#61615A", "#BA0900", "#6B7900", "#00C2A0", "#FFAA92", "#FF90C9", "#B903AA", "#D16100",
                "#DDEFFF", "#000035", "#7B4F4B", "#A1C299", "#300018", "#0AA6D8", "#013349", "#00846F",
                "#372101", "#FFB500", "#C2FFED", "#A079BF", "#CC0744", "#C0B9B2", "#C2FF99", "#001E09",
                "#00489C", "#6F0062", "#0CBD66", "#EEC3FF", "#456D75", "#B77B68", "#7A87A1", "#788D66",
                "#885578", "#FAD09F", "#FF8A9A", "#D157A0", "#BEC459", "#456648", "#0086ED", "#886F4C",
                "#34362D", "#B4A8BD", "#00A6AA", "#452C2C", "#636375", "#A3C8C9", "#FF913F", "#938A81",
                "#575329", "#00FECF", "#B05B6F", "#8CD0FF", "#3B9700", "#04F757", "#C8A1A1", "#1E6E00",
                "#7900D7", "#A77500", "#6367A9", "#A05837", "#6B002C", "#772600", "#D790FF", "#9B9700",
                "#549E79", "#FFF69F", "#201625", "#72418F", "#BC23FF", "#99ADC0", "#3A2465", "#922329",
                "#5B4534", "#FDE8DC", "#404E55", "#0089A3", "#CB7E98", "#A4E804", "#324E72", "#6A3A4C")
)
