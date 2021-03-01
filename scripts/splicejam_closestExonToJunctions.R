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
