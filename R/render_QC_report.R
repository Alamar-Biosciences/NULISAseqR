#' Render NULISAseq QC HTML Report
#'
#' Generates an HTML quality control report using the NULISAseq R Markdown template skeleton.Rmd.
#'
#' @param output_filename Output HTML filename.
#' @param output_dir Output directory for the rendered report.
#' @param study_name Title of the report.
#' @param assayName Name of the assay. Default is NULL.
#' @param dataDir Directory containing XML files.
#' @param xml_files Vector of XML filenames to include.
#' @param dataRuns Optional list of preloaded data objects. Default is NULL.
#' @param Rmd_input_file Path to the R Markdown template file skeleton.Rmd. Defaults to the template in the NULISAseqR package.
#' @param plateNames Optional vector of plate names to override defaults.
#' @param assayRunInfo Optional list of lists with assay run dates and instruments. Default is NULL.
#' @param report_type Type of report: "internal" or "webApp".
#' @param advancedQC Logical; whether to include advanced QC metrics.
#' @param IC Vector of internal control target names to override default behavior.
#' @param IPC Vector of inter-plate control sample identifiers to override default behavior.
#' @param NC Vector of negative control sample identifiers to override default behavior.
#' @param SC Vector of sample control sample identifiers to override default behavior.
#' @param Bridge Vector of bridge sample identifiers to.
#' @param rowAnnotName Column name for row annotation. Default is "AUTO_WELLROW".
#' @param colAnnotName Column name for column annotation. Default is "AUTO_WELLCOL".
#' @param plateAnnotName Column name for plate ID annotations. Default is "AUTO_PLATE".
#' @param excludeSamples List of sample names to exclude per run.
#' @param excludeTargets List of target names to exclude per run.
#' @param ignoreTargetBlank List of named vectors specifying targets and associated blank samples to ignore.
#' @param minBlankNoMAD Minimum number of blanks required for MAD outlier detection. Default is 4.
#' @param plate_effect_pval_ANOVA Unadjusted p-value threshold for ANOVA plate effect detection. Default is 0.01. 
#' @param plate_effect_pval_pairwise P-value threshold for pairwise plate effect detection. P-values are Tukey-adjusted. Default is 0.01. 
#' @param plate_effect_icc_threshold ICC \% threshold for plate effect flagging of targets. Default is 10\%. 
#' @param plate_effect_sig_pct_threshold Percent threshold of significant targets needed to flag a plate. Default is 10\%.
#' @param sampleGroupCovar Column name for sample group variable used in detectability. Default is "SAMPLE_MATRIX".
#' @param outRunSummary Logical. Include Run Summary section. Default is TRUE.
#' @param outPlateLayout Logical. Include Plate Layout section. Default is TRUE.
#' @param outReadSummary Logical. Include Read Summary section. Default is TRUE.
#' @param outHeatmaps Logical. Include Heatmaps. Default is TRUE.
#' @param outQC Logical. Include QC Summary section. Default is TRUE.
#' @param outDetectability Logical. Include Detectability Analysis. Default is TRUE.
#' @param outIntraPlateNorm Logical. Include Intra-Plate CV Analysis. Default is TRUE.
#' @param outInterPlateNorm Logical. Include Inter-Plate CV Analysis. Default is TRUE.
#' @param outSampleBoxplot Logical. Include Sample Boxplots. Default is TRUE.
#' @param outSampleCorrelation Logical. Include Sample Correlation plot. Default is TRUE.
#' @param outSampleClustering Logical. Include Sample Clustering heatmap. Default is TRUE.
#' @param outSamplePCA Logical. Include PCA plot. Default is TRUE.
#' @param outPlateEffect Logical. Include Plate Effect section. Default is TRUE.
#' @param out_SC_IPC_Ratio Logical. Include SC/IPC ratio graphic. Default is TRUE.
#' @param out_SC_NC_Ratio Logical. Include SC/NC ratio graphic. Default is TRUE.
#' @param heatMapRel Logical. Whether to use relative (\% median) values in heatmaps. Default is TRUE.
#' @param outputPlots Logical. Whether to output slide-friendly plots. Default is TRUE.
#' @param outputDetectCSV Logical. Output CSV of detectability data. Default is TRUE.
#' @param outputCoefVarCSV Logical. Output CSV of coefficient of variation results. Default is TRUE.
#' @param output_SC_IPC_Ratio_CSV Logical. Output CSV of SC/IPC log2 ratios. Default is TRUE.
#' @param output_SC_NC_Ratio_CSV Logical. Output CSV of SC/NC log2 ratios. Default is TRUE.
#' @param outputPlateEffectCSV Logical. Whether to output plate effect test CSV. Default is TRUE.
#' @param outputRData Logical. Whether to output .RData file of report data. Default is TRUE.
#' @param rendered_by_shiny Logical. Indicates if the report is rendered in a Shiny app. Default is FALSE.
#' @param highlight_TAP_report_fields Logical. Highlight fields required for TAP report. Default is TRUE.
#' @param sort_by_plateNames Logical. Whether to sort runs by plate name. Default is TRUE.
#' @param ... Additional named parameters passed to the R Markdown template.
#'
#' @return Path to the rendered report (invisible).
#' @export
render_QC_report <- function(output_filename,
                             output_dir,
                             study_name = 'NULISAseq QC_Report',
                             assayName = NULL,
                             dataDir,
                             xml_files,
                             dataRuns = NULL,
                             Rmd_input_file = file.path(system.file(package = 'NULISAseqR'), 'rmarkdown/templates/nulisaseq/skeleton/skeleton.Rmd'),
                             plateNames = NULL,
                             assayRunInfo = NULL,
                             report_type = 'internal',
                             advancedQC = FALSE,
                             IC = NULL,
                             IPC = NULL,
                             NC = NULL,
                             SC = NULL,
                             Bridge = NULL,
                             rowAnnotName = "AUTO_WELLROW",
                             colAnnotName = "AUTO_WELLCOL",
                             plateAnnotName = "AUTO_PLATE",
                             excludeSamples = NULL,
                             excludeTargets = NULL,
                             ignoreTargetBlank = NULL,
                             minBlankNoMAD = 4,
                             plate_effect_pval_ANOVA = 0.01,
                             plate_effect_pval_pairwise = 0.01,
                             plate_effect_icc_threshold = 10,
                             plate_effect_sig_pct_threshold = 10,
                             sampleGroupCovar = "SAMPLE_MATRIX",
                             outRunSummary = TRUE,
                             outPlateLayout = TRUE,
                             outReadSummary = TRUE,
                             outHeatmaps = TRUE,
                             outQC = TRUE,
                             outDetectability = TRUE,
                             outIntraPlateNorm = TRUE,
                             outInterPlateNorm = TRUE,
                             outSampleBoxplot = TRUE,
                             outSampleCorrelation = TRUE,
                             outSampleClustering = TRUE,
                             outSamplePCA = TRUE,
                             outPlateEffect = TRUE,
                             heatMapRel = TRUE,
                             outputPlots = TRUE,
                             outputDetectCSV = TRUE,
                             outputCoefVarCSV = TRUE,
                             out_SC_IPC_Ratio = TRUE,
                             out_SC_NC_Ratio = TRUE,
                             output_SC_IPC_Ratio_CSV = TRUE,
                             output_SC_NC_Ratio_CSV = TRUE,
                             outputPlateEffectCSV = TRUE,
                             outputRData = TRUE,
                             rendered_by_shiny = FALSE,
                             highlight_TAP_report_fields = TRUE,
                             sort_by_plateNames = TRUE,
                             ...
) {
  # Combine explicitly passed parameters with additional ones in ...
  input_params <- modifyList(
    list(
      study_name = study_name,
      assayName = assayName,
      dataDir = dataDir,
      xmlFiles = xml_files,
      plateNames = plateNames,
      assayRunInfo = assayRunInfo,
      reportType = report_type,
      advancedQC = advancedQC,
      IC = IC,
      IPC = IPC,
      NC = NC,
      SC = SC,
      Bridge = Bridge,
      rowAnnotName = rowAnnotName,
      colAnnotName = colAnnotName,
      plateAnnotName = plateAnnotName,
      excludeSamples = excludeSamples,
      excludeTargets = excludeTargets,
      ignoreTargetBlank = ignoreTargetBlank,
      minBlankNoMAD = minBlankNoMAD,
      plate_effect_pval_ANOVA = plate_effect_pval_ANOVA,
      plate_effect_pval_pairwise = plate_effect_pval_pairwise,
      plate_effect_icc_threshold = plate_effect_icc_threshold,
      plate_effect_sig_pct_threshold = plate_effect_sig_pct_threshold,
      sampleGroupCovar = sampleGroupCovar,
      outRunSummary = outRunSummary,
      outPlateLayout = outPlateLayout,
      outReadSummary = outReadSummary,
      outHeatmaps = outHeatmaps,
      outQC = outQC,
      outDetectability = outDetectability,
      outIntraPlateNorm = outIntraPlateNorm,
      outInterPlateNorm = outInterPlateNorm,
      outSampleBoxplot = outSampleBoxplot,
      outSampleCorrelation = outSampleCorrelation,
      outSampleClustering = outSampleClustering,
      outSamplePCA = outSamplePCA,
      outPlateEffect = outPlateEffect,
      heatMapRel = heatMapRel,
      outputPlots = outputPlots,
      outputDetectCSV = outputDetectCSV,
      outputCoefVarCSV = outputCoefVarCSV,
      out_SC_IPC_Ratio = out_SC_IPC_Ratio,
      out_SC_NC_Ratio = out_SC_NC_Ratio,
      output_SC_IPC_Ratio_CSV = output_SC_IPC_Ratio_CSV,
      output_SC_NC_Ratio_CSV = output_SC_NC_Ratio_CSV,
      outputPlateEffectCSV = outputPlateEffectCSV,
      outputRData = outputRData,
      rendered_by_shiny = rendered_by_shiny,
      highlight_TAP_report_fields = highlight_TAP_report_fields,
      sort_by_plateNames = sort_by_plateNames
    ),
    list(...)
  )
  
  rmarkdown::render(
    input = Rmd_input_file,
    output_format = 'html_document',
    output_file = output_filename,
    output_dir = output_dir,
    params = input_params
  )
  
  invisible(file.path(output_dir, output_filename))
}