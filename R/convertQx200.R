################################################################################
#                  DDES QX200 Conversion                                       #
################################################################################

#' Convert QX200 QuantaSoft exports to DDES
#'
#' Converts a QX200 QuantaSoft set of head file and well amplitude files to
#' DDES main, assay, and intensity file archive.
#'
#' @param headFile (str) Path to QX200 QuantaSoft head csv file.
#' @param amplitudeFiles (str) Vector containing paths to QX200 amplitude 
#'   csv files.
#' @param outputDirectory (str) Directory to which the .ddes files are written.
#' @param experimentId (str) Experiment identifier, written in the DDES headers.
#' @param runId Run identifier written in the DDES headers. By default uses
#'   the head filename without its extension.
#' @param instrument (str) Instrument name, written in the DDES headers.
#' @param instrumentSoftwareVersion (str) Instrument software version, written
#'   DDES headers.
#' @param ddesVersion (str) DDES format version written in the DDES headers.
#' @param dateTime (str) Timestamp written to the headers and filenames. Default
#'   uses the current time, `YYYYMMDDHHMM` format.
#'
#' @return (list) containing the converted data & headers, file and DDES paths.
#'
#' @examples
#' exampleDirectory <- system.file(
#'   "extdata",
#'   "qx200_quantasoft",
#'   package = "DDES"
#' )
#' outputDirectory <- file.path(tempdir(), "DDES-QX200-example")
#'
#' conversion <- convertQx200(
#'   headFile = file.path(
#'     exampleDirectory,
#'     "20231130_JVS230_BRIP_Duplex_Normal.csv"
#'   ),
#'   amplitudeFiles = list.files(
#'     exampleDirectory,
#'     pattern = "_Amplitude[.]csv$",
#'     full.names = TRUE
#'   ),
#'   outputDirectory = outputDirectory,
#'   experimentId = "JVS230_BRIP_Duplex_Normal",
#'   instrumentSoftwareVersion = "QuantaSoft",
#'   dateTime = "202611051300"
#' )
#'
#' conversion$archive
#'
#' @export
convertQx200 <- function(
  headFile,
  amplitudeFiles,
  outputDirectory,
  experimentId = "expID",
  runId = tools::file_path_sans_ext(basename(headFile)),
  instrument = "QX200",
  instrumentSoftwareVersion = "NA",
  ddesVersion = "1.0",
  dateTime = format(Sys.time(), "%Y%m%d%H%M")
) {
  print("Reading QX200 QuantaSoft export files")

  headData <- readQx200Head(headFile)
  wellIds <- unique(headData$Well)

  # Align QuantaSoft targets with the color ids used in DDES intensity files
  headData$channel_ID <- sub("Unknown$", "", headData$TargetType)

  
  print("Creating the DDES assay file")

  wellProbeCombinations <- lapply(wellIds, function(wellId) {

    wellRows <- headData[headData$Well == wellId, , drop = FALSE]
    unique(paste(wellRows$Target, wellRows$channel_ID, sep = "."))
  })

  wellCombinationKeys <- vapply(
    wellProbeCombinations,
    paste,
    character(1),
    collapse = ";"
  )

  wellAssayIds <- vapply(wellIds, function(wellId) {

    wellTargets <- unique(
      headData[headData$Well == wellId, "Target"]
    )
    paste(wellTargets, collapse = "_")
  }, character(1))

  assayCombinations <- unique(wellCombinationKeys)
  assayCombinationRows <- match(assayCombinations, wellCombinationKeys)

  assayData <- data.frame(
    assay_ID = wellAssayIds[assayCombinationRows],
    probe_combination = assayCombinations,
    probe_concentrations_nM = vapply(
      strsplit(assayCombinations, split = ";"),
      function(probeNames) paste(rep("NA", length(probeNames)), collapse = ";"),
      character(1)
    ),
    stringsAsFactors = FALSE
  )

  assayHeader <- c(
    "#header",
    paste("#DDES_version", ddesVersion, sep = "\t"),
    paste("#datetime", dateTime, sep = "\t"),
    paste("#experiment_ID", experimentId, sep = "\t"),
    paste("#run_ID", runId, sep = "\t"),
    paste("#DDES_type", "assays", sep = "\t"),
    "#data",
    "#assay_information"
  )

  assayEditorData <- assayData
  assayEditorData$target_number <- lengths(
    strsplit(assayEditorData$probe_combination, split = ";")
  )
  assayEditorData$selector <- FALSE

  probeData <- unique(
    headData[, c("Target", "channel_ID"), drop = FALSE]
  )
  colnames(probeData)[1] <- "target_ID"
  probeData$selector <- FALSE
  rownames(probeData) <- seq_len(nrow(probeData))

  wellAssayMap <- data.frame(
    well_ID = wellIds,
    assay_ID = wellAssayIds,
    stringsAsFactors = FALSE
  )

  
  print("Creating the DDES main file")

  mainData <- headData[
    ,
    c(
      "Well",
      "Sample",
      "Target",
      "Positives",
      "Negatives",
      "Concentration",
      "Threshold"
    ),
    drop = FALSE
  ]

  mainData$assay_ID <- wellAssayMap$assay_ID[
    match(mainData$Well, wellAssayMap$well_ID)
  ]
  mainData <- mainData[
    ,
    c(
      "Well",
      "Sample",
      "assay_ID",
      "Target",
      "Positives",
      "Negatives",
      "Concentration",
      "Threshold"
    )
  ]
  colnames(mainData) <- c(
    "well_ID",
    "sample_ID",
    "assay_ID",
    "target_ID",
    "counts_positive",
    "counts_negative",
    "concentration_reaction_cp_\u00B5L",
    "threshold_if_threshold"
  )

  mainHeader <- c(
    "#header",
    paste("#DDES_version", ddesVersion, sep = "\t"),
    paste("#datetime", dateTime, sep = "\t"),
    paste("#experiment_ID", experimentId, sep = "\t"),
    paste("#run_ID", runId, sep = "\t"),
    paste("#DDES_type", "main", sep = "\t"),
    paste("#instrument", instrument, sep = "\t"),
    paste("#plate_type", "Droplets", sep = "\t"),
    paste(
      "#instrument_software_version",
      instrumentSoftwareVersion,
      sep = "\t"
    ),
    "#data"
  )

  sampleAnnotation <- unique(
    mainData[, c("well_ID", "assay_ID", "sample_ID")]
  )
  sampleAnnotation$selector <- FALSE
  rownames(sampleAnnotation) <- seq_len(nrow(sampleAnnotation))

  
  
  print("Creating the DDES intensity files")

  intensityData <- vector("list", length(wellIds))
  intensityHeaders <- vector("list", length(wellIds))
  intensityFiles <- character(length(wellIds))

  for (wellNumber in seq_along(wellIds)) {

    wellId <- wellIds[wellNumber]
    amplitudeFile <- grep(
      pattern = paste0("_", wellId, "_"),
      x = amplitudeFiles,
      value = TRUE
    )[1]
    amplitudeData <- readQx200Amplitude(amplitudeFile)
    partitionIds <- seq_len(nrow(amplitudeData))

    wellIntensity <- data.frame(
      partition_ID = rep(partitionIds, times = 2),
      channel_ID = rep(
        c("Ch1", "Ch2"),
        each = length(partitionIds)
      ),
      FI_endpoint = c(
        amplitudeData$Ch1.Amplitude,
        amplitudeData$Ch2.Amplitude
      ),
      stringsAsFactors = FALSE
    )

    intensityName <- paste0(
      "DDES_",
      experimentId,
      "_",
      runId,
      "_",
      dateTime,
      "_intensity_",
      wellId
    )

    intensityData[[wellNumber]] <- wellIntensity
    names(intensityData)[wellNumber] <- intensityName

    assayId <- wellAssayMap[
      wellAssayMap$well_ID == wellId,
      "assay_ID"
    ]

    intensityHeaders[[wellNumber]] <- c(
      "#head",
      paste("#DDES_version", ddesVersion, sep = "\t"),
      paste("#datetime", dateTime, sep = "\t"),
      paste("#experiment_ID", experimentId, sep = "\t"),
      paste("#run_ID", runId, sep = "\t"),
      paste("#DDES_type", "intensity", sep = "\t"),
      paste("#well_ID", wellId, sep = "\t"),
      paste("#assay_ID", assayId, sep = "\t"),
      "#data"
    )
    names(intensityHeaders)[wellNumber] <- intensityName
    intensityFiles[wellNumber] <- paste0(intensityName, ".tsv")
  }
  

  print("Writing the DDES files")

  dir.create(outputDirectory, recursive = TRUE, showWarnings = FALSE)

  mainFile <- file.path(
    outputDirectory,
    paste0(
      "DDES_",
      experimentId,
      "_",
      runId,
      "_",
      dateTime,
      "_main.tsv"
    )
  )
  assayFile <- file.path(
    outputDirectory,
    paste0(
      "DDES_",
      experimentId,
      "_",
      runId,
      "_",
      dateTime,
      "_assays.tsv"
    )
  )

  writeDdesTable(mainHeader, mainData, mainFile)
  writeDdesTable(assayHeader, assayData, assayFile)

  for (wellNumber in seq_along(wellIds)) {

    intensityFile <- file.path(
      outputDirectory,
      intensityFiles[wellNumber]
    )
    writeDdesTable(
      intensityHeaders[[wellNumber]],
      intensityData[[wellNumber]],
      intensityFile
    )
  }

  outputFiles <- c(
    basename(mainFile),
    basename(assayFile),
    intensityFiles
  )
  archiveFile <- file.path(
    outputDirectory,
    paste0(
      "DDES_",
      experimentId,
      "_",
      runId,
      "_",
      dateTime,
      ".ddes"
    )
  )

  zip::zip(
    zipfile = archiveFile,
    files = outputFiles,
    root = outputDirectory,
    mode = "cherry-pick"
  )

  invisible(
    list(
      mainHeader = mainHeader,
      mainData = mainData,
      assayHeader = assayHeader,
      assayData = assayEditorData,
      probeData = probeData,
      sampleAnnotation = sampleAnnotation,
      intensityHeaders = intensityHeaders,
      intensityData = intensityData,
      files = file.path(outputDirectory, outputFiles),
      archive = archiveFile
    )
  )
}


## Read a QX200 QuantaSoft head file, can hold a separator line
readQx200Head <- function(file) {
  firstLine <- readLines(file, n = 1, encoding = "UTF-8")
  firstLine <- sub("^\ufeff", "", firstLine)
  skipLines <- if (grepl("^sep=", firstLine, ignore.case = TRUE)) 1 else 0

  headData <- utils::read.csv(
    file = file,
    header = TRUE,
    skip = skipLines,
    row.names = NULL,
    stringsAsFactors = FALSE,
    check.names = TRUE
  )

  if (colnames(headData)[1] == "row.names") {
    colnames(headData) <- c(colnames(headData)[-1], "unused")
  }
  headData
}


## Read a QX200 amplitude file
readQx200Amplitude <- function(file) {
  utils::read.csv(
    file = file,
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = TRUE
  )
}
