################################################################################
#                    DDES QIAcuity Conversion                                  #
################################################################################

#' Convert QIAcuity exports to DDES
#'
#' Converts a QIAcuity set of analysis file, plate layout file, and relative
#' fluorescence unit (RFU) zip to DDES main, assay, and intensity file archive.
#'
#' @param analysisFile (str) Path to QIAcuity analysis csv file.
#' @param plateLayoutFile (str) Path to QIAcuity plate-layout csv file.
#' @param rfuArchive (str) Path to zip containing the QIAcuity RFU files.
#' @param outputDirectory (str) Directory to which the .ddes files are written.
#' @param experimentId (str) Experiment identifier, written in the DDES headers.
#' @param instrument (str) Instrument name, written in the DDES headers.
#' @param instrumentSoftwareVersion (str) Instrument software version, written
#'   in the DDES headers.
#' @param ddesVersion (str) DDES format version written in the DDES headers.
#' @param dateTime (str) Timestamp written to the headers and filenames. Default
#'   uses the current time, `YYYYMMDDHHMM` format.
#'
#' @return (list) containing the converted data & headers, file and DDES paths.
#'
#' @examples
#' exampleDirectory <- system.file("extdata", "qiacuity", package = "DDES")
#' outputDirectory <- file.path(tempdir(), "DDES-QIAcuity-example")
#'
#' conversion <- convertQiacuity(
#'   analysisFile = file.path(exampleDirectory, "current-results.csv"),
#'   plateLayoutFile = file.path(exampleDirectory, "plate-layout.csv"),
#'   rfuArchive = file.path(exampleDirectory, "rfu-export.zip"),
#'   outputDirectory = outputDirectory,
#'   experimentId = "tutorialExperiment",
#'   instrument = "QIAcuity Four",
#'   instrumentSoftwareVersion = "QIAcuity Software Suite 2.5.0",
#'   dateTime = "202506301127"
#' )
#'
#' @export
convertQiacuity <- function(
  analysisFile,
  plateLayoutFile,
  rfuArchive,
  outputDirectory,
  experimentId = "expID",
  instrument = "NA",
  instrumentSoftwareVersion = "NA",
  ddesVersion = "1.0",
  dateTime = format(Sys.time(), "%Y%m%d%H%M")
) {

  print("Reading QIAcuity files")

  # Read analysis and plate layout
  currentResults <- readQiacuityCsv(analysisFile)
  plateLayout <- readQiacuityCsv(plateLayoutFile)

  # Extract the RFU files in a temp directory, and combine
  rfuDirectory <- tempfile("qiacuity-rfu-")
  dir.create(rfuDirectory)
  utils::unzip(rfuArchive, exdir = rfuDirectory)

  rfuFiles <- list.files(
    rfuDirectory,
    pattern = "RFU.*\\.csv$",
    full.names = TRUE,
    ignore.case = TRUE
  )

  rfuData <- do.call(
    rbind,
    lapply(rfuFiles, readQiacuityRfu)
  )

  # Retrieve run metadata
  plateInfo <- currentResults[1, 1:3, drop = FALSE]
  currentResults <- currentResults[, -c(1:3), drop = FALSE]

  runId <- as.character(plateInfo$Plate.ID[1])
  plateType <- as.character(plateInfo$Plate.type[1])

  print("Creating the DDES main file")

  mainHeader <- c(
    "#header",
    paste("#DDES_version", ddesVersion, sep = "\t"),
    paste("#datetime", dateTime, sep = "\t"),
    paste("#experiment_ID", experimentId, sep = "\t"),
    paste("#run_ID", runId, sep = "\t"),
    paste("#DDES_type", "main", sep = "\t"),
    paste("#instrument", instrument, sep = "\t"),
    paste("#plate_type", plateType, sep = "\t"),
    paste(
      "#instrument_software_version",
      instrumentSoftwareVersion,
      sep = "\t"
    ),
    "#data"
  )

  concentrationColumn <- grep(
    pattern = "^Conc",
    x = colnames(currentResults),
    value = TRUE
  )[1]

  mainColumns <- c(
    "Well",
    "Sample.NTC.Control",
    "Reaction.Mix",
    "Target..Name.",
    "Partitions..Positive.",
    "Partitions..Negative.",
    concentrationColumn,
    "Threshold"
  )

  mainData <- currentResults[, mainColumns, drop = FALSE]
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
  mainData <- replace(mainData, mainData == "NA", "-")

  sampleAnnotation <- unique(
    mainData[, c("well_ID", "assay_ID", "sample_ID")]
  )
  sampleAnnotation$selector <- FALSE
  rownames(sampleAnnotation) <- seq_len(nrow(sampleAnnotation))
  wellIds <- sampleAnnotation$well_ID

  print("Creating the DDES assay file")

  assayIds <- unique(plateLayout$Reaction.mix.name)
  assayIds <- assayIds[!is.na(assayIds) & assayIds != ""]

  assayData <- do.call(
    rbind,
    lapply(assayIds, function(assayId) {

      assayRows <- plateLayout[
        plateLayout$Reaction.mix.name == assayId,
        ,
        drop = FALSE
      ]

      probeCombinations <- character()

      for (channelNumber in seq_len(5)) {

        channelName <- as.character(
          assayRows[[paste0("Channel.", channelNumber)]]
        )
        targetName <- as.character(
          assayRows[[paste0("Target.", channelNumber)]]
        )

        validProbe <- !is.na(channelName) & !is.na(targetName) &
          channelName != "" & targetName != ""

        probeCombinations <- c(
          probeCombinations,
          paste(targetName[validProbe], channelName[validProbe], sep = ".")
        )
      }

      probeCombinations <- unique(probeCombinations)

      data.frame(
        assay_ID = assayId,
        probe_combination = paste(probeCombinations, collapse = ";"),
        probe_concentrations_nM = paste(
          rep("NA", length(probeCombinations)),
          collapse = ";"
        ),
        stringsAsFactors = FALSE
      )
    })
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

  probeNames <- unique(
    unlist(strsplit(assayData$probe_combination, split = ";"))
  )
  probeParts <- strsplit(probeNames, split = "\\.")
  probeData <- data.frame(
    target_ID = vapply(probeParts, `[[`, character(1), 1),
    channel_ID = vapply(probeParts, `[[`, character(1), 2),
    selector = FALSE,
    stringsAsFactors = FALSE
  )

  print("Creating the DDES intensity files")

  intensityData <- vector("list", length(wellIds))
  intensityHeaders <- vector("list", length(wellIds))
  intensityFiles <- character(length(wellIds))

  channelNames <- c(
    C = "CRIMSON",
    G = "GREEN",
    Y = "YELLOW",
    R = "RED",
    O = "ORANGE"
  )

  for (wellNumber in seq_along(wellIds)) {

    wellId <- wellIds[wellNumber]
    assayName <- plateLayout[
      plateLayout$Well == wellId,
    ]$Reaction.mix.name[1]

    wellIntensity <- rfuData[rfuData$Well == wellId, , drop = FALSE]
    wellIntensity <- wellIntensity[, c("Partition", "Channel", "RFU")]
    colnames(wellIntensity) <- c(
      "partition_ID",
      "channel_ID",
      "FI_endpoint"
    )
    wellIntensity <- wellIntensity[!is.na(wellIntensity$FI_endpoint), ]
    wellIntensity$channel_ID <- unname(
      channelNames[wellIntensity$channel_ID]
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

    intensityHeaders[[wellNumber]] <- c(
      "#head",
      paste("#DDES_version", ddesVersion, sep = "\t"),
      paste("#datetime", dateTime, sep = "\t"),
      paste("#experiment_ID", experimentId, sep = "\t"),
      paste("#run_ID", runId, sep = "\t"),
      paste("#DDES_type", "intensity", sep = "\t"),
      paste("#well_ID", wellId, sep = "\t"),
      paste("#assay_ID", assayName, sep = "\t"),
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


## Read a QIAcuity csv file, check first for prsence of separator line
readQiacuityCsv <- function(file) {
  firstLine <- readLines(file, n = 1, encoding = "UTF-8")
  firstLine <- sub("^\ufeff", "", firstLine)
  skipLines <- if (grepl("^sep=", firstLine, ignore.case = TRUE)) 1 else 0

  utils::read.csv(
    file = file,
    header = TRUE,
    skip = skipLines,
    stringsAsFactors = FALSE,
    check.names = TRUE
  )
}


## Read QIAcuity RFU files, always skip the separator line
readQiacuityRfu <- function(file) {
  utils::read.csv(
    file = file,
    header = TRUE,
    skip = 1,
    stringsAsFactors = FALSE
  )
}


## Write a DDES header and data table
writeDdesTable <- function(header, data, file) {
  fileConnection <- file(file, open = "wt", encoding = "UTF-8")
  writeLines(header, fileConnection)
  utils::write.table(
    data,
    file = fileConnection,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
  )
  close(fileConnection)
}
