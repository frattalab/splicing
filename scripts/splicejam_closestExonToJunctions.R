
###pasteByRow

pasteByRow = function (x, sep = "_", na.rm = TRUE, condenseBlanks = TRUE,
    includeNames = FALSE, sepName = ":", blankGrep = "^[ ]*$",
    verbose = FALSE, ...)
{
    sep <- head(sep, 1)
    if (length(ncol(x)) == 0 || ncol(x) == 0) {
        return(x)
    }
    if (igrepHas("matrix", class(x))) {
        x <- as.data.frame(x)
    }
    for (iCol in seq_len(ncol(x))) {
        if (igrepHas("factor", class(x[, iCol]))) {
            x[, iCol] <- as.character(x[[iCol]])
        }
    }
    getColVals <- function(x, i, includeNames, na.rm, sepName) {
        xVals <- x[[i]]
        isNa <- (is.na(xVals))
        if (any(isNa)) {
            if (na.rm) {
                xVals[isNa] <- ""
            }
            else {
                xVals[isNa] <- "NA"
            }
        }
        if (condenseBlanks) {
            isBlank <- grep(blankGrep, xVals)
        }
        if (includeNames) {
            xVals <- paste0(colnames(x)[i], sepName, xVals)
            if (condenseBlanks && length(isBlank) > 0) {
                xVals[isBlank] <- ""
            }
        }
        else {
            if (condenseBlanks && length(isBlank) > 0) {
                xVals[isBlank] <- ""
            }
        }
        xVals
    }
    xVals <- getColVals(x, 1, includeNames, na.rm, sepName)
    if (ncol(x) > 1) {
        for (i1 in 2:ncol(x)) {
            xVals1 <- getColVals(x, i1, includeNames, na.rm,
                sepName)
            if (condenseBlanks) {
                isBlank1 <- (is.na(xVals1) | grepl(blankGrep,
                  xVals1))
                isBlank <- (is.na(xVals) | grepl(blankGrep, xVals))
                sepV <- ifelse(isBlank | isBlank1, "", sep)
            }
            else {
                sepV <- sep
            }
            xVals <- paste0(xVals, sepV, xVals1)
        }
    }
    if (!is.null(rownames(x))) {
        names(xVals) <- rownames(x)
    }
    return(xVals)
}
###function taken from https://github.com/jmw86069/splicejam - please cite them!###
splicejam_closestExonToJunctions = function (spliceGRgene, exonsGR, flipNegativeStrand = TRUE, sampleColname = "sample_id",
    reportActualCoords = FALSE, verbose = FALSE, ...)
{
    updateColnames <- c("distFrom", "distTo", "nameFrom", "nameTo",
        "genesDiffer", "genesMatch", "tooFarFrom", "tooFarTo",
        "tooFar")
    if (any(updateColnames %in% colnames(GenomicRanges::values(spliceGRgene)))) {
        GenomicRanges::values(spliceGRgene) <- GenomicRanges::values(spliceGRgene)[,
            setdiff(colnames(spliceGRgene), updateColnames),
            drop = FALSE]
    }
    if (verbose) {
        printDebug("closestExonToJunctions(): ", "Finding closest exons for splice starts.")
    }
    spliceStartExonEndD1 <- as.data.frame(GenomicRanges::distanceToNearest(GenomicRanges::resize(spliceGRgene,
        width = 1, fix = "start", ignore.strand = TRUE), GenomicRanges::resize(exonsGR,
        width = 1, fix = "end", ignore.strand = TRUE), select = "all"))
    if (verbose) {
        printDebug("closestExonToJunctions(): ", "Calculating stranded distance.")
        print(spliceStartExonEndD1)
    }
    spliceStartExonEndDactual1 <- (GenomicRanges::start(GenomicRanges::resize(spliceGRgene,
        width = 1, fix = "start", ignore.strand = TRUE)[spliceStartExonEndD1[,
        "queryHits"]]) - GenomicRanges::start(GenomicRanges::resize(exonsGR,
        width = 1, fix = "end", ignore.strand = TRUE)[spliceStartExonEndD1[,
        "subjectHits"]]))
    spliceStartExonEndDactual1end <- GenomicRanges::end(exonsGR[spliceStartExonEndD1[,
        "subjectHits"]])
    spliceStartExonEndDactual <- spliceStartExonEndDactual1 -
        sign(spliceStartExonEndDactual1)
    spliceGRgeneNeg <- as.vector(strand(spliceGRgene)) %in% "-"
    if (any(spliceGRgeneNeg)) {
        spliceStartExonEndDactual[spliceGRgeneNeg] <- spliceStartExonEndDactual[spliceGRgeneNeg] *
            -1
    }
    if (verbose) {
        printDebug("closestExonToJunctions(): ", "Finding closest exons for splice ends.")
    }
    spliceEndExonStartD1 <- as.data.frame(GenomicRanges::distanceToNearest(GenomicRanges::resize(spliceGRgene,
        width = 1, fix = "end", ignore.strand = TRUE), GenomicRanges::resize(exonsGR,
        width = 1, fix = "start", ignore.strand = TRUE), select = "all"))
    if (verbose) {
        printDebug("closestExonToJunctions(): ", "Calculating stranded distance.")
    }
    spliceEndExonStartDactual1 <- (GenomicRanges::start(GenomicRanges::resize(spliceGRgene,
        width = 1, fix = "end", ignore.strand = TRUE)[spliceEndExonStartD1[,
        "queryHits"]]) - GenomicRanges::start(GenomicRanges::resize(exonsGR,
        width = 1, fix = "start", ignore.strand = TRUE)[spliceEndExonStartD1[,
        "subjectHits"]]))
    spliceEndExonStartDactual1start <- GenomicRanges::start(exonsGR[spliceEndExonStartD1[,
        "subjectHits"]])
    spliceEndExonStartDactual <- spliceEndExonStartDactual1 -
        sign(spliceEndExonStartDactual1)
    if (any(spliceGRgeneNeg)) {
        spliceEndExonStartDactual[spliceGRgeneNeg] <- spliceEndExonStartDactual[spliceGRgeneNeg] *
            -1
    }
    spliceStartExonEndD1[, "strandedDistance"] <- spliceStartExonEndDactual
    spliceEndExonStartD1[, "strandedDistance"] <- spliceEndExonStartDactual
    GenomicRanges::values(spliceGRgene)[spliceStartExonEndD1[,
        "queryHits"], "nameFrom"] <- names(exonsGR[spliceStartExonEndD1[,
        "subjectHits"]])
    GenomicRanges::values(spliceGRgene)[spliceEndExonStartD1[,
        "queryHits"], "nameTo"] <- names(exonsGR[spliceEndExonStartD1[,
        "subjectHits"]])
    GenomicRanges::values(spliceGRgene)[, "nameFromTo"] <- pasteByRow(GenomicRanges::values(spliceGRgene)[,
        c("nameFrom", "nameTo")], sep = " ", na.rm = TRUE)
    GenomicRanges::values(spliceGRgene)[spliceStartExonEndD1[,
        "queryHits"], "distFrom"] <- spliceStartExonEndD1[, "strandedDistance"]
    GenomicRanges::values(spliceGRgene)[spliceEndExonStartD1[,
        "queryHits"], "distTo"] <- spliceEndExonStartD1[, "strandedDistance"]
    if (reportActualCoords) {
        if (verbose) {
            printDebug("closestExonToJunctions(): ", "Reporting actual coords.")
        }
        GenomicRanges::values(spliceGRgene)[spliceStartExonEndD1[,
            "queryHits"], "coordFrom"] <- spliceStartExonEndDactual1end
        GenomicRanges::values(spliceGRgene)[spliceEndExonStartD1[,
            "queryHits"], "coordTo"] <- spliceEndExonStartDactual1start
        if (flipNegativeStrand) {
        }
    }
    if (flipNegativeStrand) {
        if (verbose) {
            printDebug("closestExonToJunctions(): ", "Flipping negative strand.")
        }
        fromToCols <- paste0(rep(c("dist", "name", "coord", "tooFar"),
            each = 2), c("From", "To"))
        switchCols1 <- intersect(fromToCols, colnames(GenomicRanges::values(spliceGRgene)))
        switchCols2 <- as.vector(matrix(nrow = 2, switchCols1)[2:1,
            ])
        if (verbose) {
            printDebug("closestExonToJunctions(): ", "switchCols1:",
                switchCols1)
            printDebug("closestExonToJunctions(): ", "switchCols2:",
                switchCols2)
        }
        negStrand <- (as.vector(strand(spliceGRgene)) %in% "-")
        if (any(negStrand)) {
            GenomicRanges::values(spliceGRgene)[negStrand, switchCols1] <- GenomicRanges::values(spliceGRgene)[negStrand,
                switchCols2]
        }
    }
    retVal <- list(spliceStartExonEndD = spliceStartExonEndD1,
        spliceEndExonStartD = spliceEndExonStartD1, spliceGRgene = spliceGRgene)
    return(retVal)
}
