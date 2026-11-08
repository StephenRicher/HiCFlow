nextflow.enable.dsl = 2

include { FASTQC_RAW } from './modules/local/fastqc'
include { CUTADAPT } from './modules/local/cutadapt'
include { FASTQC_TRIMMED } from './modules/local/fastqc_trimmed'
include { HICUP_TRUNCATE } from './modules/local/hicup_truncate'
include { PREPARE_GENOME } from './modules/local/prepare_genome'
include { BOWTIE2_BUILD } from './modules/local/bowtie2_build'
include { BOWTIE2_ALIGN } from './modules/local/bowtie2_align'
include { PREPARE_BAM } from './modules/local/prepare_bam'
include { PROCESS_HIC_QC } from './modules/local/process_hic_qc'
include { PLOT_HIC_QC } from './modules/local/plot_hic_qc'
include { SAMTOOLS_STATS } from './modules/local/samtools_stats'
include { SPLIT_PAIRED_READS } from './modules/local/split_paired_reads'
include { BUILD_BASE_MATRIX } from './modules/local/build_base_matrix'
include { ADJUST_MATRIX } from './modules/local/adjust_matrix'
include { SUM_REPLICATES } from './modules/local/sum_replicates'
include { MERGE_MATRIX_BINS } from './modules/local/merge_matrix_bins'
include { MATRIX_COVERAGE } from './modules/local/matrix_coverage'
include { HICREP } from './modules/local/hicrep'
include { PLOT_HICREP } from './modules/local/plot_hicrep'
include { NORMALISE_MATRIX } from './modules/local/normalise_matrix'
include { DETECT_LOOPS } from './modules/local/detect_loops'
include { SCORE_LOOP_DIFFERENCES } from './modules/local/score_loop_differences'
include { FIND_TADS } from './modules/local/find_tads'
include { ONTAD } from './modules/local/ontad'
include { CSCORE } from './modules/local/cscore'
include { DISTANCE_NORMALISE } from './modules/local/distance_normalise'
include { PLOT_OBS_EXP_MATRIX } from './modules/local/plot_obs_exp_matrix'
include { COMPARE_MATRICES } from './modules/local/compare_matrices'
include { PLOT_COMPARE_VIEWPOINT } from './modules/local/plot_compare_viewpoint'
include { PREPARE_GENES } from './modules/local/prepare_genes'
include { PLOT_HIC } from './modules/local/plot_hic'
include { PROCESS_DIFF_TADS } from './modules/local/process_diff_tads'
include { PLOT_SUBTRACTION } from './modules/local/plot_subtraction'
include { PLOT_VIEWPOINT } from './modules/local/plot_viewpoint'
include { MULTIQC } from './modules/local/multiqc'


def createSamplesChannel(samplesheetPath) {
    def samplesheet = file(samplesheetPath, checkIfExists: true)
    def samplesheetDir = samplesheet.parent

    Channel
        .fromPath(samplesheet, checkIfExists: true)
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            def missingValues = ['experiment', 'cell_type', 'group', 'rep', 'R1', 'R2']
                .findAll { column -> !row[column] }

            if (missingValues) {
                error "Missing required values in ${samplesheet}: ${missingValues.join(', ')}"
            }

            def sampleId = "${row.group}-${row.rep}"
            def read1 = samplesheetDir.resolve(row.R1).normalize()
            def read2 = samplesheetDir.resolve(row.R2).normalize()

            if (!read1.exists()) {
                error "R1 FASTQ does not exist for ${sampleId}: ${read1}"
            }

            if (!read2.exists()) {
                error "R2 FASTQ does not exist for ${sampleId}: ${read2}"
            }

            def meta = [
                id         : sampleId,
                experiment : row.experiment,
                cell_type  : row.cell_type,
                group      : row.group,
                replicate  : row.rep as Integer
            ]

            tuple(meta, read1, read2)
        }
}

def createRegionsChannel(regionsPath) {
    Channel
        .fromPath(regionsPath, checkIfExists: true)
        .splitCsv(header: false, sep: '\t')
        .map { row ->
            if (row.size() != 4) {
                error "Expected four columns in ${regionsPath}, found: ${row}"
            }

            def region = [
                chromosome : row[0],
                start      : row[1] as Integer,
                end        : row[2] as Integer,
                name       : row[3]
            ]

            if (region.start < 0 || region.end <= region.start) {
                error "Invalid region coordinates: ${region}"
            }

            if (region.name ==~ /.*[-.\/].*/) {
                error "Region name contains a prohibited character: ${region.name}"
            }

            region
        }
}

def createPlotCoordinatesChannel(coordinatesPath) {
    Channel
        .fromPath(coordinatesPath, checkIfExists: true)
        .splitText()
        .map { line ->
            def stripped = line.trim()

            if (!stripped || stripped.startsWith('#')) {
                return null
            }

            def fields = stripped.split(/\s+/)

            if (fields.size() != 4) {
                error "Expected four columns in ${coordinatesPath}, found: ${stripped}"
            }

            def start = fields[1] as Integer
            def end = fields[2] as Integer

            [
                chromosome: fields[0],
                start     : start,
                end       : end,
                name      : fields[3],
                id        : "${fields[3]}_${fields[0]}_${start}_${end}"
            ]
        }
        .filter { it != null }
        .unique { coordinate -> coordinate.id }
}

workflow {
    if (!params.samplesheet) {
        error 'Missing required parameter: --samplesheet'
    }

    if (!params.genome) {
        error 'Missing required genome configuration'
    }

    if (!params.regions) {
        error 'Missing required parameter: --regions'
    }

    if (!params.genes?.gff3) {
        error 'Missing gene annotation: params.genes.gff3'
    }

    log.info """
    HiCFlow Nextflow
    =================
    Sample sheet: ${params.samplesheet}
    Output directory: ${params.outdir}
    """

    samples = createSamplesChannel(params.samplesheet)

    genomes = Channel.fromList(
        params.genome.collect { cellType, referencePath ->
            tuple(
                cellType,
                file(referencePath, checkIfExists: true)
            )
        }
    )

    PREPARE_GENOME(
        genomes,
        params.hicup.restrictionSequence
    )

    genomesForIndex = PREPARE_GENOME.out.genome.map {
        cellType,
        reference,
        referenceIndex,
        chromSizes,
        restrictionSites ->

        tuple(cellType, reference)
    }

    geneSettings = [
        genomeBuild: params.genomeBuild,
        typeKey    : params.genes.typeKey,
        label      : params.genes.label,
        geneID     : params.genes.geneID
    ]

    PREPARE_GENES(
        file(params.genes.gff3, checkIfExists: true),
        geneSettings
    )

    geneAnnotation = PREPARE_GENES.out.annotation

    bigWigFiles = Channel.value(
        (params.bigWig ?: [:]).collect { bigWigTitle, bigWigPath ->
            file(bigWigPath, checkIfExists: true)
        }
    )

    BOWTIE2_BUILD(genomesForIndex)

    FASTQC_RAW(samples)

    CUTADAPT(
        samples,
        params.cutadapt
    )

    FASTQC_TRIMMED(CUTADAPT.out.reads)

    HICUP_TRUNCATE(
        CUTADAPT.out.reads,
        params.hicup
    )

    readsByCellType = HICUP_TRUNCATE.out.reads.map { meta, read1, read2 ->
        tuple(meta.cell_type, meta, read1, read2)
    }

    alignmentInputs = readsByCellType
        .combine(BOWTIE2_BUILD.out.index, by: 0)
        .map { cellType, meta, read1, read2, indexFiles ->
            tuple(cellType, meta, read1, read2, indexFiles)
        }

    BOWTIE2_ALIGN(alignmentInputs)

    PREPARE_BAM(BOWTIE2_ALIGN.out.alignments)

    PROCESS_HIC_QC(
        PREPARE_BAM.out.bam,
        params.qc.samplePairs
    )

    processedQcFiles = PROCESS_HIC_QC.out.processed
        .map { meta, processedFile ->
            processedFile
        }
        .collect()

    plotQcSettings = [
        filetype: params.plot.filetype,
        dpi     : params.plot.dpi
    ]

    PLOT_HIC_QC(
        processedQcFiles,
        plotQcSettings
    )

    SPLIT_PAIRED_READS(PREPARE_BAM.out.bam)

    splitReadsByCellType = SPLIT_PAIRED_READS.out.reads.map {
        meta,
        read1Bam,
        read2Bam ->

        tuple(
            meta.cell_type,
            meta,
            read1Bam,
            read2Bam
        )
    }

    genomeResources = PREPARE_GENOME.out.genome.map {
        cellType,
        reference,
        referenceIndex,
        chromSizes,
        restrictionSites ->

        tuple(
            cellType,
            chromSizes,
            restrictionSites
        )
    }

    matrixInputs = splitReadsByCellType
        .combine(genomeResources, by: 0)
        .map {
            cellType,
            meta,
            read1Bam,
            read2Bam,
            chromSizes,
            restrictionSites ->

            tuple(
                meta,
                read1Bam,
                read2Bam,
                chromSizes,
                restrictionSites
            )
        }

    matrixSettings = [
        baseResolution       : params.resolution.base,
        restrictionSequence  : params.hicup.restrictionSequence.replace('^', ''),
        danglingSequence     : params.hicup.restrictionSequence.replace('^', ''),
        maxLibraryInsertSize : params.hic.maxLibraryInsertSize,
        minDistance          : params.hic.minDistance,
        minMappingQuality    : params.hic.minMappingQuality,
        keepSelfLigation     : params.hic.keepSelfLigation,
        keepSelfCircles      : params.hic.keepSelfCircles,
        inputBufferSize      : 400000
    ]

    BUILD_BASE_MATRIX(
        matrixInputs,
        matrixSettings
    )

    regions = createRegionsChannel(params.regions)

    regionMatrixInputs = BUILD_BASE_MATRIX.out.matrices
        .combine(regions)
        .map { meta, matrix, region ->
            tuple(meta, matrix, region)
        }

    ADJUST_MATRIX(regionMatrixInputs)

    replicateMatricesByGroup = ADJUST_MATRIX.out.matrices
        .map { meta, region, matrix ->
            def groupMeta = [
                id        : meta.group,
                experiment: meta.experiment,
                cell_type : meta.cell_type,
                group     : meta.group
            ]

            tuple(
                groupMeta,
                region,
                matrix
            )
        }
        .groupTuple(by: [0, 1])

    SUM_REPLICATES(replicateMatricesByGroup)

    if (params.hic.compartmentScore) {
        cscoreMatrices = SUM_REPLICATES.out.matrices
            .flatMap { meta, region, baseMatrix ->
                params.resolution.bins.collect { targetResolution ->
                    tuple(
                        meta.cell_type,
                        meta,
                        region,
                        targetResolution as Integer,
                        baseMatrix
                    )
                }
            }

        cscoreGenomes = PREPARE_GENOME.out.genome.map {
            cellType,
            reference,
            referenceIndex,
            chromSizes,
            restrictionSites ->

            tuple(
                cellType,
                reference,
                referenceIndex
            )
        }

        cscoreInputs = cscoreMatrices
            .combine(cscoreGenomes, by: 0)
            .map {
                cellType,
                meta,
                region,
                resolution,
                baseMatrix,
                reference,
                referenceIndex ->

                tuple(
                    meta,
                    region,
                    resolution,
                    baseMatrix,
                    reference,
                    referenceIndex
                )
            }

        CSCORE(cscoreInputs)
    }

    allBaseMatrices = ADJUST_MATRIX.out.matrices.mix(
        SUM_REPLICATES.out.matrices
    )

    rebinnedMatrixInputs = allBaseMatrices.flatMap { meta, region, matrix ->
        params.resolution.bins.collect { targetResolution ->
            def resolution = targetResolution as Integer
            def baseResolution = params.resolution.base as Integer

            if (resolution <= baseResolution) {
                error "Target resolution ${resolution} must be greater than base resolution ${baseResolution}"
            }

            if (resolution % baseResolution != 0) {
                error "Target resolution ${resolution} is not divisible by base resolution ${baseResolution}"
            }

            tuple(
                meta,
                region,
                baseResolution,
                resolution,
                matrix
            )
        }
    }

    MERGE_MATRIX_BINS(rebinnedMatrixInputs)

    coverageInputs = MERGE_MATRIX_BINS.out.matrices
        .map { meta, region, resolution, matrix ->
            tuple(
                "${meta.id}|${region.name}",
                meta,
                region,
                resolution,
                matrix
            )
        }
        .groupTuple(by: 0)
        .map { key, metadata, regions, resolutions, matrices ->
            def orderedEntries = [resolutions, matrices]
                .transpose()
                .sort { first, second ->
                    (first[0] as Integer) <=> (second[0] as Integer)
                }

            tuple(
                metadata[0],
                regions[0],
                orderedEntries.collect { entry -> entry[0] as Integer },
                orderedEntries.collect { entry -> entry[1] }
            )
        }

    coverageSettings = [
        filetype: params.plot.filetype,
        dpi     : params.plot.dpi,
        nBins   : params.qc.coverageBins,
        fontSize: params.plot.coverageFontSize
    ]

    MATRIX_COVERAGE(
        coverageInputs,
        coverageSettings
    )

    if (params.qc.runHiCRep) {
        hicRepMatrices = MERGE_MATRIX_BINS.out.matrices
            .filter { meta, region, resolution, matrix ->
                meta.containsKey('replicate') &&
                    resolution == params.qc.hicRepBin
            }
            .map { meta, region, resolution, matrix ->
                tuple(
                    meta.cell_type,
                    meta.group,
                    region.name,
                    resolution,
                    meta,
                    region,
                    matrix
                )
            }
            .groupTuple(by: [0, 1, 2, 3])

        hicRepInputs = hicRepMatrices.flatMap {
            cellType,
            group,
            regionName,
            resolution,
            metadata,
            regions,
            matrices ->

            def entries = metadata.indices.collect { index ->
                [
                    meta  : metadata[index],
                    region: regions[index],
                    matrix: matrices[index]
                ]
            }.sort { first, second ->
                first.meta.replicate <=> second.meta.replicate
            }

            def pairs = []

            for (int first = 0; first < entries.size(); first++) {
                for (int second = first + 1; second < entries.size(); second++) {
                    pairs << tuple(
                        entries[first].meta,
                        entries[second].meta,
                        entries[first].region,
                        resolution,
                        entries[first].matrix,
                        entries[second].matrix
                    )
                }
            }

            pairs
        }

        HICREP(hicRepInputs)

        hicRepPlotInputs = HICREP.out.scores
            .map { meta1, meta2, region, resolution, scoreFile ->
                tuple(
                    "${region.name}|${resolution}",
                    region,
                    resolution,
                    scoreFile
                )
            }
            .groupTuple(by: 0)
            .map { key, regions, resolutions, scoreFiles ->
                tuple(
                    regions[0],
                    resolutions[0] as Integer,
                    scoreFiles
                )
            }

        hicRepPlotSettings = [
            filetype: params.plot.filetype,
            dpi     : params.plot.dpi
        ]

        PLOT_HICREP(
            hicRepPlotInputs,
            hicRepPlotSettings
        )
    }

    NORMALISE_MATRIX(
        MERGE_MATRIX_BINS.out.matrices,
        params.hic.multiplicativeValue
    )

    if (params.loops.detectLoops) {
        DETECT_LOOPS(
            NORMALISE_MATRIX.out.corrected,
            params.loops
        )
    }

    viewpointMatrices = NORMALISE_MATRIX.out.corrected.map {
        meta,
        region,
        resolution,
        matrix ->

        tuple(
            region.name,
            meta,
            region,
            resolution,
            matrix
        )
    }

    viewpoints = createPlotCoordinatesChannel(params.plot.viewpoints)
        .map { viewpoint ->
            tuple(viewpoint.name, viewpoint)
        }

    viewpointInputs = viewpointMatrices
        .combine(viewpoints, by: 0)
        .map {
            regionName,
            meta,
            region,
            resolution,
            matrix,
            viewpoint ->

            tuple(
                meta,
                region,
                resolution,
                matrix,
                viewpoint
            )
        }

    viewpointSettings = [
        viewpointRange: params.plot.viewpointRange,
        genomeBuild   : params.genomeBuild,
        filetype      : params.plot.filetype,
        dpi           : params.plot.dpi
    ]

    PLOT_VIEWPOINT(
        viewpointInputs,
        viewpointSettings
    )

    FIND_TADS(NORMALISE_MATRIX.out.corrected)

    ONTAD(NORMALISE_MATRIX.out.corrected)

    correctedForPlots = NORMALISE_MATRIX.out.corrected.map {
        meta,
        region,
        resolution,
        matrix ->

        tuple(
            "${meta.id}|${region.name}|${resolution}",
            meta,
            region,
            resolution,
            matrix
        )
    }

    insulationForPlots = FIND_TADS.out.results.map {
        meta,
        region,
        resolution,
        boundariesBed,
        boundariesGff,
        domainsBed,
        insulation,
        zscoreMatrix ->

        tuple(
            "${meta.id}|${region.name}|${resolution}",
            insulation
        )
    }

    ontadForPlots = ONTAD.out.domains.map {
        meta,
        region,
        resolution,
        tads ->

        tuple(
            "${meta.id}|${region.name}|${resolution}",
            tads
        )
    }

    plotData = correctedForPlots
        .join(insulationForPlots)
        .join(ontadForPlots)
        .map {
            key,
            meta,
            region,
            resolution,
            matrix,
            insulation,
            tads ->

            tuple(
                region.name,
                meta,
                region,
                resolution,
                matrix,
                insulation,
                tads
            )
        }

    plotCoordinates = createPlotCoordinatesChannel(params.plot.coordinates)
        .map { coordinate ->
            tuple(coordinate.name, coordinate)
        }

    plotInputs = plotData
        .combine(plotCoordinates, by: 0)
        .map {
            regionName,
            meta,
            region,
            resolution,
            matrix,
            insulation,
            tads,
            coordinate ->

            tuple(
                meta,
                region,
                resolution,
                matrix,
                insulation,
                tads,
                coordinate
            )
        }

    DISTANCE_NORMALISE(NORMALISE_MATRIX.out.corrected)

    plotSettings = [
        filetype   : params.plot.filetype,
        dpi        : params.plot.dpi,
        genomeBuild: params.genomeBuild
    ]

    PLOT_OBS_EXP_MATRIX(
        DISTANCE_NORMALISE.out.matrices,
        plotSettings
    )

    groupMatricesByRegion = MERGE_MATRIX_BINS.out.matrices
        .filter { meta, region, resolution, matrix ->
            !meta.containsKey('replicate')
        }
        .map { meta, region, resolution, matrix ->
            tuple(
                meta.cell_type,
                region.name,
                resolution,
                meta,
                region,
                matrix
            )
        }
        .groupTuple(by: [0, 1, 2])

    matrixComparisons = groupMatricesByRegion.flatMap {
        cellType,
        regionName,
        resolution,
        metadata,
        regions,
        matrices ->

        def entries = [metadata, matrices]
            .transpose()
            .sort { first, second ->
                first[0].id <=> second[0].id
            }

        def comparisons = []

        for (int first = 0; first < entries.size(); first++) {
            for (int second = first + 1; second < entries.size(); second++) {
                comparisons << tuple(
                    cellType,
                    entries[first][0],
                    entries[second][0],
                    regions[0],
                    resolution,
                    entries[first][1],
                    entries[second][1]
                )
            }
        }

        comparisons
    }

    chromSizesByCellType = PREPARE_GENOME.out.genome.map {
        cellType,
        reference,
        referenceIndex,
        chromSizes,
        restrictionSites ->

        tuple(cellType, chromSizes)
    }

    comparisonInputs = matrixComparisons
        .combine(chromSizesByCellType, by: 0)
        .map {
            cellType,
            group1,
            group2,
            region,
            resolution,
            matrix1,
            matrix2,
            chromSizes ->

            tuple(
                group1,
                group2,
                region,
                resolution,
                matrix1,
                matrix2,
                chromSizes
            )
        }

    COMPARE_MATRICES(
        comparisonInputs,
        params.compareMatrices
    )

    if (params.loops.detectLoops) {
        groupLoops = DETECT_LOOPS.out.loops
            .filter { meta, region, resolution, loops ->
                !meta.containsKey('replicate')
            }
            .map { meta, region, resolution, loops ->
                tuple(
                    "${meta.id}|${region.name}|${resolution}",
                    loops
                )
            }

        loopComparisonMatrices = COMPARE_MATRICES.out.differences.map {
            group1,
            group2,
            region,
            resolution,
            unfilteredMatrix,
            filteredMatrix,
            subtractionMetadata ->

            tuple(
                "${group1.id}|${region.name}|${resolution}",
                group1,
                group2,
                region,
                resolution,
                filteredMatrix
            )
        }

        loopInputsWithGroup1 = loopComparisonMatrices
            .join(groupLoops)

        loopInputs = loopInputsWithGroup1
            .map {
                group1Key,
                group1,
                group2,
                region,
                resolution,
                subtractionMatrix,
                group1Loops ->

                tuple(
                    "${group2.id}|${region.name}|${resolution}",
                    group1,
                    group2,
                    region,
                    resolution,
                    subtractionMatrix,
                    group1Loops
                )
            }
            .join(groupLoops)
            .map {
                group2Key,
                group1,
                group2,
                region,
                resolution,
                subtractionMatrix,
                group1Loops,
                group2Loops ->

                tuple(
                    group1,
                    group2,
                    region,
                    resolution,
                    subtractionMatrix,
                    group1Loops,
                    group2Loops
                )
            }

        SCORE_LOOP_DIFFERENCES(
            loopInputs,
            params.loops
        )
    }

    adjustedMatricesForViewpoints = COMPARE_MATRICES.out.adjusted.flatMap {
        group1,
        group2,
        region,
        resolution,
        adjustedMatrix1,
        adjustedMatrix2 ->

        def comparison = [
            id    : "${group1.id}-vs-${group2.id}",
            group1: group1.id,
            group2: group2.id
        ]

        [
            tuple(
                region.name,
                comparison,
                'adjIF1',
                region,
                resolution,
                adjustedMatrix1
            ),
            tuple(
                region.name,
                comparison,
                'adjIF2',
                region,
                resolution,
                adjustedMatrix2
            )
        ]
    }

    compareViewpointInputs = adjustedMatricesForViewpoints
        .combine(viewpoints, by: 0)
        .map {
            regionName,
            comparison,
            adjustedSet,
            region,
            resolution,
            matrix,
            viewpoint ->

            tuple(
                comparison,
                adjustedSet,
                region,
                resolution,
                matrix,
                viewpoint
            )
        }

    compareViewpointSettings = [
        viewpointRange: params.plot.viewpointRange,
        genomeBuild   : params.genomeBuild,
        filetype      : params.plot.filetype,
        dpi           : params.plot.dpi
    ]

    PLOT_COMPARE_VIEWPOINT(
        compareViewpointInputs,
        compareViewpointSettings
    )

    comparisonDifferences = COMPARE_MATRICES.out.differences.map {
        group1,
        group2,
        region,
        resolution,
        unfilteredMatrix,
        filteredMatrix,
        subtractionMetadata ->

        def comparisonKey = "${group1.id}|${group2.id}|${region.name}|${resolution}"

        tuple(
            comparisonKey,
            group1,
            group2,
            region,
            resolution,
            filteredMatrix
        )
    }

    rawGroupMatrices = MERGE_MATRIX_BINS.out.matrices
        .filter { meta, region, resolution, matrix ->
            !meta.containsKey('replicate')
        }
        .map { meta, region, resolution, matrix ->
            tuple(
                "${meta.id}|${region.name}|${resolution}",
                matrix
            )
        }

    groupTads = ONTAD.out.domains
        .filter { meta, region, resolution, tads ->
            !meta.containsKey('replicate')
        }
        .map { meta, region, resolution, tads ->
            tuple(
                "${meta.id}|${region.name}|${resolution}",
                tads
            )
        }

    diffTadInputsWithGroup1Raw = comparisonDifferences
        .map { comparisonKey, group1, group2, region, resolution, filteredMatrix ->
            tuple(
                "${group1.id}|${region.name}|${resolution}",
                comparisonKey,
                group1,
                group2,
                region,
                resolution,
                filteredMatrix
            )
        }
        .join(rawGroupMatrices)

    diffTadInputsWithRawMatrices = diffTadInputsWithGroup1Raw
        .map {
            group1Key,
            comparisonKey,
            group1,
            group2,
            region,
            resolution,
            filteredMatrix,
            rawMatrix1 ->

            tuple(
                "${group2.id}|${region.name}|${resolution}",
                comparisonKey,
                group1,
                group2,
                region,
                resolution,
                filteredMatrix,
                rawMatrix1
            )
        }
        .join(rawGroupMatrices)

    diffTadInputsWithGroup1Tads = diffTadInputsWithRawMatrices
        .map {
            group2Key,
            comparisonKey,
            group1,
            group2,
            region,
            resolution,
            filteredMatrix,
            rawMatrix1,
            rawMatrix2 ->

            tuple(
                "${group1.id}|${region.name}|${resolution}",
                comparisonKey,
                group1,
                group2,
                region,
                resolution,
                filteredMatrix,
                rawMatrix1,
                rawMatrix2
            )
        }
        .join(groupTads)

    diffTadInputs = diffTadInputsWithGroup1Tads
        .map {
            group1Key,
            comparisonKey,
            group1,
            group2,
            region,
            resolution,
            filteredMatrix,
            rawMatrix1,
            rawMatrix2,
            tads1 ->

            tuple(
                "${group2.id}|${region.name}|${resolution}",
                group1,
                group2,
                region,
                resolution,
                filteredMatrix,
                rawMatrix1,
                rawMatrix2,
                tads1
            )
        }
        .join(groupTads)
        .map {
            group2Key,
            group1,
            group2,
            region,
            resolution,
            filteredMatrix,
            rawMatrix1,
            rawMatrix2,
            tads1,
            tads2 ->

            tuple(
                group1,
                group2,
                region,
                resolution,
                filteredMatrix,
                rawMatrix1,
                rawMatrix2,
                tads1,
                tads2
            )
        }

    PROCESS_DIFF_TADS(
        diffTadInputs,
        params.compareMatrices.zThreshold
    )

    subtractionMatricesForPlots = COMPARE_MATRICES.out.differences.map {
        group1,
        group2,
        region,
        resolution,
        unfilteredMatrix,
        filteredMatrix,
        subtractionMetadata ->

        tuple(
            "${group1.id}|${group2.id}|${region.name}|${resolution}",
            group1,
            group2,
            region,
            resolution,
            filteredMatrix
        )
    }

    diffTadsForPlots = PROCESS_DIFF_TADS.out.allTads.map {
        group1,
        group2,
        region,
        resolution,
        group1Tads,
        group2Tads ->

        tuple(
            "${group1.id}|${group2.id}|${region.name}|${resolution}",
            group1Tads
        )
    }

    subtractionPlotData = subtractionMatricesForPlots
        .join(diffTadsForPlots)
        .map {
            comparisonKey,
            group1,
            group2,
            region,
            resolution,
            subtractionMatrix,
            tads ->

            tuple(
                region.name,
                group1,
                group2,
                region,
                resolution,
                subtractionMatrix,
                tads
            )
        }

    subtractionPlotInputs = subtractionPlotData
        .combine(plotCoordinates, by: 0)
        .map {
            regionName,
            group1,
            group2,
            region,
            resolution,
            subtractionMatrix,
            tads,
            coordinate ->

            tuple(
                group1,
                group2,
                region,
                resolution,
                subtractionMatrix,
                tads,
                coordinate
            )
        }

    subtractionPlotSettings = [
        genomeBuild      : params.genomeBuild,
        colourmap        : params.compareMatrices.colourmap,
        vMax             : params.compareMatrices.vMax,
        filetype         : params.plot.filetype,
        dpi              : params.plot.dpi,
        geneLabelFontSize: params.plot.geneLabelFontSize,
        maxLabels        : params.plot.maxLabels,
        bigWigTitles     : (params.bigWig ?: [:]).keySet() as List,
        vlines           : params.plot.vlines
            ? file(params.plot.vlines, checkIfExists: true).toString()
            : null
    ]

    PLOT_SUBTRACTION(
        subtractionPlotInputs,
        geneAnnotation,
        bigWigFiles,
        subtractionPlotSettings
    )

    plotHiCSettings = [
        genomeBuild      : params.genomeBuild,
        colourmap        : params.plot.colourmap,
        filetype         : params.plot.filetype,
        dpi              : params.plot.dpi,
        geneLabelFontSize: params.plot.geneLabelFontSize,
        maxLabels        : params.plot.maxLabels,
        bigWigTitles     : (params.bigWig ?: [:]).keySet() as List,
        vlines           : params.plot.vlines
            ? file(params.plot.vlines, checkIfExists: true).toString()
            : null
    ]

    PLOT_HIC(
        plotInputs,
        geneAnnotation,
        bigWigFiles,
        plotHiCSettings
    )

    SAMTOOLS_STATS(PREPARE_BAM.out.bam)

    rawFastQcFiles = FASTQC_RAW.out.reports.flatMap { meta, htmlFiles, zipFiles ->
        [htmlFiles, zipFiles].flatten()
    }

    trimmedFastQcFiles = FASTQC_TRIMMED.out.reports.flatMap { meta, htmlFiles, zipFiles ->
        [htmlFiles, zipFiles].flatten()
    }

    cutadaptFiles = CUTADAPT.out.reports.flatMap { meta, reports ->
        [reports].flatten()
    }

    bowtieFiles = BOWTIE2_ALIGN.out.reports.flatMap { meta, read1Report, read2Report ->
        [read1Report, read2Report].flatten()
    }

    hicupFiles = HICUP_TRUNCATE.out.reports.flatMap { meta, reports ->
        [reports].flatten()
    }

    deduplicationFiles = PREPARE_BAM.out.reports.flatMap { meta, reports ->
        [reports].flatten()
    }

    samtoolsFiles = SAMTOOLS_STATS.out.reports.flatMap { meta, reports ->
        [reports].flatten()
    }

    hicExplorerQcDirectories = BUILD_BASE_MATRIX.out.qc.map { meta, qcDirectory ->
        qcDirectory
    }

    multiqcInputs = rawFastQcFiles
        .mix(trimmedFastQcFiles)
        .mix(cutadaptFiles)
        .mix(bowtieFiles)
        .mix(hicupFiles)
        .mix(deduplicationFiles)
        .mix(samtoolsFiles)
        .mix(hicExplorerQcDirectories)
        .collect()

    MULTIQC(
        multiqcInputs,
        file(params.qc.multiqcConfig, checkIfExists: true)
    )
}
