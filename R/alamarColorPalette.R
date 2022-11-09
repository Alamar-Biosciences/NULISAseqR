
#' Make color palette based on Alamar colors
#'
#'@param n number of main colors to output (maximum 13)
#'@param nReps optional, number of sub-colors to output 
#' for each of the main colors. Useful for when there are technical replicates.
#' Function will output n*nReps total colors.
#' @param tint optional, used when subcolors are created.
#' "light" blends main colors with white.
#' "dark" blends main colors with black.
#' "lightdark" (default) blends main colors with white & black.
#' @param interpolate logical TRUE or FALSE (default). 
#' If TRUE, main colors are based on interpolation of 
#' the colors denoted by the interpolate_index 
#' argument.
#' @param interpolateIndex indices of colors to interpolate. Default is 1:5
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
                               interpolateIndex=1:5){
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
