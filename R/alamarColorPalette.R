
#' Make color palette based on Alamar colors
#'
#'@param n number of main colors to output. For `palette=1`, maximum is 13 
#' when `interpolate=FALSE`, and unlimited when `interpolate=TRUE`. For `palette=2`, 
#' maximum is 9 when `interpolate=FALSE`, and unlimited when `interpolate=TRUE`
#'@param nReps optional, number of sub-colors to output 
#' for each of the main colors. Useful for when there are technical replicates.
#' Function will output n*nReps total colors.
#' @param tint optional, used when sub-colors are created.
#' `"light"` blends main colors with white.
#' `"dark"` blends main colors with black.
#' `"lightdark"` (default) blends main colors with white & black.
#' @param interpolate logical `TRUE` or `FALSE` (default). 
#' If `TRUE`, main colors are based on interpolation of 
#' the colors denoted by the interpolate_index 
#' argument.
#' @param interpolateIndex indices of colors to interpolate. For `palette=1`, 
#' default is 1:5. For `palette=2`, default is `c(1,2,4,3,5,6)`.
#' @param palette Which set of colors to use. Integer. `"palette=1"` (default) 
#' is based on the original website color palette. `palette=2` is the new 
#' palette used in the NULISA manuscript.
#'
#' @return A vector or list of vectors of colors in hex format
#'
#' @examples
#' cols <- alamarColorPalette(n=13, nReps=3)
#' plot(1:length(unlist(cols)),
#' 1:length(unlist(cols)), col=unlist(cols),
#' pch=19, cex=4)
#' 
#' cols <- alamarColorPalette(n=10, nReps=3,
#' interpolate=TRUE, interpolateIndex=1:5)
#' plot(1:length(unlist(cols)),
#' 1:length(unlist(cols)), col=unlist(cols),
#' pch=19, cex=4)
#' 
#'
#' @export

alamarColorPalette <- function(n, 
                               nReps=1, 
                               tint='lightdark',
                               interpolate=FALSE,
                               interpolateIndex=NULL,
                               palette=1){
  if (palette==1){
    if (n > 13 & interpolate==FALSE){
      stop('number of colors cannot exceed 13 unless interpolate==TRUE')
    }
    allColors <- c('#FFBD03', # yellow  
                   '#00C28C', # light green 
                   '#4C00BF', # dark purple 
                   '#DA4E37', # orange 
                   '#46AADC', # light blue 
                   '#D51C8D', # hot pink 
                   '#5809B1', # royal purple 
                   '#1D53C2', # dark blue 
                   '#0975EF', # medium blue 
                   '#079AAF', # teal 
                   '#05C3DE', # light teal 
                   '#3B8DFD', # cornflower blue
                   '#888888') # dark grey
    if(is.null(interpolateIndex)) interpolateIndex <- 1:5
  } else if (palette==2){
    if (n > 9 & interpolate==FALSE){
      stop('number of colors cannot exceed 9 unless interpolate==TRUE')
    }
    allColors <- c('#005292', # blue (secondary color)
                   '#C91846', # red (secondary color)
                   '#99CC33', # green (secondary color)
                   '#FFCE00', # yellow (main color)
                   '#00A0DF', # light blue (main color)
                   '#12284C', # dark navy (text color)
                   '#232323', # dark grey (text color)
                   '#5A595C', # medium grey (text color)
                   '#BABCBE') # light grey (text color)
    if(is.null(interpolateIndex)) interpolateIndex <- c(1,2,4,3,5,6)
  }
  if (interpolate==FALSE){
    mainColors <- allColors[1:n]
  } else if (interpolate==TRUE){
    mainColors <- colorRampPalette(allColors[interpolateIndex])(n)
  }
  colors <- mainColors
  if (nReps > 1){
    colors <- lapply(mainColors, function(x){
      if (tint=='light'){
        subcolors <- colorRampPalette(c(x, '#FFFFFF'))(nReps + 1)
        subcolors[1:nReps]
      } else if (tint=='dark'){
        subcolors <- colorRampPalette(c('#000000', x))(nReps + 1)
        subcolors[2:(nReps+1)]
      } else if (tint=='lightdark'){
        subcolors <- colorRampPalette(c('#000000', x, '#FFFFFF'))(nReps + 2)
        subcolors[2:(nReps+1)]
      }
    })
  }
  return(colors)
}
