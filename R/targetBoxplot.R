#' Target Boxplot
#'
#' Makes boxplot of log2-transformed counts for each target. 
#' Option to show target distributions relative to the LOD (LOD is 
#' subtracted from each data point, so LOD becomes the zero) and color
#' boxplots by target detectability.
#'
#' @param data_matrix The Data matrix output from readNULISAseq.R
#' or normalized data from normalization functions (targets in columns 
#' and samples in rows).
#' @param title Title of the plot.
#' @param sortBy Order by which to sort the boxplots. 
#' Options are `'median'` for target medians, `'detect'` for
#' detectability (must provide LODs). If `subtractLOD = FALSE`, will only
#' sort by medians.
#' @param subtractLOD Should the target boxplot values be relative to LOD?
#' Default is `FALSE`. If `TRUE`, will color boxplots by detectability.
#' @param blanks A vector of column names that represent the negative controls.
#' Required when `subtractLOD = TRUE`. This is passed to the `lod()` function.
#' @param horizontal Should the boxplots be horizontal (target names along the 
#' y-axis)? Default is `FALSE` (target names along x-axis).
#' @param log2transform If TRUE (default), each value `x` will be transformed 
#' by `log2(x + 1)`. LOD values will also be transformed when 
#' `subtractLOD = TRUE`.
#' @param replace_zero_LOD If TRUE (default is FALSE), will replace LODs of zero 
#' on the unlogged scale with zeros on the log scale. Probably only makes sense 
#' to use this option when using unnormalized data, because in this case LODs of 
#' zero become negative and when subtracted will cause shift upwards for the set 
#' of targets with LOD = 0. This option is only used when log2transform is TRUE.
#' @param axis_lab_normalized If TRUE (default), will use "normalized count" for axis
#' label. If FALSE, will use "count" for axis label.
#' @param excludeTargets Vector of names of any targets to exclude from the plots 
#' (such as internal controls).
#' @param excludeSamples Vector of names of samples to exclude from plots. These 
#' correspond to the column names. 
#' @param colors Optional vector of colors that will be used to color the 
#' boxplots. Colors should be in the same order as the targets (rows) 
#' in the input data (after excluding any targets in `excludeTargets`). 
#' If provided, will override `detectHighColor` 
#' and `detectLowColor` when `subtractLOD = TRUE`. 
#' When `subtractLOD = FALSE`, default is `'grey'`
#' @param detectHighColor The color that represents 100% detectability in the 
#' gradient color scale used when `subtractLOD = TRUE`. Default is `'pink'`.
#' @param detectLowColor The color that represents 0% detectability in the 
#' gradient color scale used when `subtractLOD = TRUE`. Default is `'blue'`.
#' @param showLegend Should a legend appear in the plot? Only used when 
#' `subtractLOD = TRUE`. 
#' @param cex.targets Character expansion factor for target labels.
#' 
#' 
#'
#' @return Draws boxplots of target expression.
#'
#' 
#' @export
#' 
targetBoxplot <- function(data_matrix,
                          title=NULL,
                          sortBy = 'median',
                          subtractLOD=FALSE,
                          blanks=NULL,
                          horizontal=FALSE,
                          log2transform=TRUE,
                          replace_zero_LOD=FALSE,
                          axis_lab_normalized=TRUE,
                          excludeTargets=NULL,
                          excludeSamples=NULL,
                          colors=NULL,
                          detectHighColor='pink',
                          detectLowColor='blue',
                          showLegend=TRUE,
                          cex.targets=0.4){
  # exclude targets
  if(!is.null(excludeTargets)){
    data_matrix <- data_matrix[!(rownames(data_matrix) %in% excludeTargets),]
  }
  
  # boxplot if not subtracting LOD
  if(subtractLOD==FALSE){
    
    # make title
    if(is.null(title)){
      title <- 'Target distributions'
    }
    
    # exclude samples
    if(!is.null(excludeSamples)){
      data_matrix <- data_matrix[,!(colnames(data_matrix) %in% excludeSamples)]
    }
    
    # log2 transform
    if(log2transform==TRUE){
      data_matrix <- log2(data_matrix + 1)
      if(axis_lab_normalized==FALSE) axis_label <- 'log2(count)'
      else if (axis_lab_normalized==TRUE) axis_label <- 'log2(normalized count)'
    } else if(log2transform==FALSE){
      if(axis_lab_normalized==FALSE) axis_label <- 'count'
      else if (axis_lab_normalized==TRUE) axis_label <- 'normalized count'
    }
    
    # sort targets
    target_medians <- apply(data_matrix, 1, median)
    target_medians <- target_medians[order(target_medians, decreasing=TRUE)]
    data_matrix <- data_matrix[names(target_medians),]
    
    # use grey if colors not provided
    if(is.null(colors)){
      colors <- 'grey'
    }
    
    if(horizontal==FALSE){
      # make plot
      boxplot(t(data_matrix), 
              horizontal=FALSE,
              xlab='',
              ylab=axis_label,
              main=title,
              las=3, 
              lty=1,
              pch=16, 
              cex=0.5,
              staplelty=0,
              cex.axis=cex.targets,
              yaxt='n',
              col=colors)
      axis(side=2, las=1)
      grid(nx=NA, ny=NULL)
      
    } else if(horizontal==TRUE){
      # make plot
      boxplot(t(apply(data_matrix, 2, rev)), 
              horizontal=TRUE,
              xlab=axis_label,
              ylab='',
              main=title,
              las=1, 
              lty=1,
              pch=16, 
              cex=0.5,
              staplelty=0,
              cex.axis=cex.targets,
              xaxt='n',
              col=rev(colors))
      grid(nx=NULL, ny=NA)
      axis(side=1, las=1)
    } # end horizontal==TRUE
  } # end subtractLOD==FALSE
  
  # boxplot if subtracting LOD (detectability plot)
  if(subtractLOD==TRUE){
    
    # make title
    if(is.null(title)){
      title <- 'Target distributions relative to LOD'
    }
    
    # calculate LODs
    LOD <- lod(data_matrix,
               blanks=blanks)
    
    # exclude samples
    if(!is.null(excludeSamples)){
      data_matrix <- data_matrix[,!(colnames(data_matrix) %in% excludeSamples)]
    }
    
    # calculate detectability
    detect <- detectability(aboveLOD_matrix = LOD$aboveLOD[,!(colnames(LOD$aboveLOD) %in% excludeSamples)])$all$detectability
    # log2 transform
    if(log2transform==TRUE){
      data_matrix <- log2(data_matrix + 1)
      LOD_count_scale <- LOD$LOD
      LOD$LOD <- log2(LOD$LOD + 1)
      if(replace_zero_LOD==TRUE){
        # replace zero LOD_count_scale with zero
        LOD$LOD[LOD_count_scale == 0] <- 0
      }
      if(axis_lab_normalized==FALSE) axis_label <- 'log2(count) - LOD'
      else if (axis_lab_normalized==TRUE) axis_label <- 'log2(normalized count) - LOD'
    } else if(log2transform==FALSE){
      if(axis_lab_normalized==FALSE) axis_label <- 'count - LOD'
      else if (axis_lab_normalized==TRUE) axis_label <- 'normalized count - LOD'
    }
    
    # subtract the LODs from data
    data_matrix <- sweep(data_matrix, MARGIN=1, STATS=LOD$LOD)
    
    # sort targets
    if(sortBy=='median'){
      target_medians <- apply(data_matrix, 1, median)
      target_medians <- target_medians[order(target_medians, decreasing=TRUE)]
      detect <- detect[names(target_medians)]
      data_matrix <- data_matrix[names(detect),]
    } else if(sortBy=='detect'){
      detect <- detect[order(detect, decreasing=TRUE)]
      data_matrix <- data_matrix[names(detect),]
    }
    
    # create colors to indicate detectability
    makeColors <- colorRamp(colors=c(detectLowColor, detectHighColor))
    detectColors <- makeColors(detect/100)
    detectColors <- apply(detectColors, 1, 
                          function(x) rgb(x[1], x[2], x[3], maxColorValue=255))
    
    # if colors are provided, override the detectColors
    if(!is.null(colors)){
      detectColors <- colors
    }
    
    if(horizontal==FALSE){
      # make plot
      par(mar=c(4,4,2,0.5))
      boxplot(t(data_matrix), 
              horizontal=FALSE,
              xlab='',
              ylab=axis_label,
              main=title,
              las=3, 
              lty=1,
              pch=16, 
              cex=0.5,
              staplelty=0,
              cex.axis=cex.targets,
              yaxt='n',
              col=detectColors)
      axis(side=2, las=1)
      grid(nx=NA, ny=NULL)
      abline(h=0, col='darkred')
      axis_limits <- par('usr')
      text(x=axis_limits[1] + (axis_limits[2] - axis_limits[1])/20, y=0,
           labels='LOD', col='darkred', cex=0.8, pos=3)
      
      # make legend
      if(showLegend==TRUE){
        # add color gradient legend
        # define coordinates
        axis_limits <- par('usr')
        x_axis_length <- axis_limits[2] - axis_limits[1]
        y_axis_length <- axis_limits[4] - axis_limits[3]
        legend_length <- x_axis_length/6
        legend_right <- axis_limits[2] - x_axis_length/20
        legend_left <- legend_right - legend_length
        legend_height <- y_axis_length/20
        legend_top <- axis_limits[4] - y_axis_length/10
        legend_bottom <- legend_top - legend_height
        # create color gradient
        colfunc <- colorRampPalette(c(detectLowColor, detectHighColor))
        legend_image <- as.raster(matrix(colfunc(20), nrow=1))
        rasterImage(legend_image, legend_left, legend_bottom, legend_right, legend_top)
        # make labels
        text(x=seq(legend_left, legend_right, l=6), 
             y = legend_bottom - y_axis_length/40,
             labels = seq(0,100,l=6), 
             cex=0.8, adj=0.5)
        text(x=(legend_left+legend_right)/2, y=legend_top, 
             labels='detectability (%)', pos=3, cex=0.8)
      } # end legend
      
    } else if(horizontal==TRUE){
      # make plot
      boxplot(t(apply(data_matrix, 2, rev)), 
              horizontal=TRUE,
              xlab=axis_label,
              ylab='',
              main=title,
              las=1, 
              lty=1,
              pch=16, 
              cex=0.5,
              staplelty=0,
              cex.axis=cex.targets,
              xaxt='n',
              col=rev(detectColors))
      axis(side=1, las=1)
      grid(nx=NULL, ny=NA)
      abline(v=0, col='darkred')
      axis_limits <- par('usr')
      text(x=0, y=axis_limits[4] - (axis_limits[4] - axis_limits[3])/20,
           labels='LOD', col='darkred', cex=0.8, pos=4)
      
      # make legend
      if(showLegend==TRUE){
        # add color gradient legend
        # define coordinates
        axis_limits <- par('usr')
        x_axis_length <- axis_limits[2] - axis_limits[1]
        y_axis_length <- axis_limits[4] - axis_limits[3]
        legend_length <- x_axis_length/25
        legend_right <- axis_limits[2] - x_axis_length/9
        legend_left <- legend_right - legend_length
        legend_height <- y_axis_length/6
        legend_bottom <- axis_limits[3] + y_axis_length/20
        legend_top <- legend_bottom + legend_height
        # create color gradient
        colfunc <- colorRampPalette(c(detectHighColor, detectLowColor))
        legend_image <- as.raster(matrix(colfunc(20), ncol=1))
        rasterImage(legend_image, legend_left, legend_bottom, legend_right, legend_top)
        # make labels
        text(x = legend_left + x_axis_length/15, 
             y = seq(legend_bottom, legend_top, l=6),
             labels = seq(0,100, l=6), 
             cex=0.8, adj=0.5)
        text(x=(legend_left+legend_right)/2, y=legend_top, 
             labels='detectability (%)', pos=3, cex=0.8)
      } # end legend
    } # end horizontal==TRUE
  } # end subtractLOD==TRUE 
}