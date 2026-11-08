process COMPARE_MATRICES {
    tag "${group2.id}-minus-${group1.id}:${region.name}:${resolution}"

    label 'process_medium'

    publishDir "${params.outdir}/matrix/compare/${region.name}/${resolution}",
        mode: 'copy',
        pattern: '*.h5'

    publishDir "${params.outdir}/matrix/compare/${region.name}/${resolution}",
        mode: 'copy',
        pattern: '*-subtraction.tsv'

    publishDir "${params.outdir}/matrix/compare/${region.name}/${resolution}/text",
        mode: 'copy',
        pattern: '*.sutm|*.2d.txt'

    publishDir "${params.outdir}/qc/hiccompare/${region.name}/${resolution}",
        mode: 'copy',
        pattern: '*.png'

    publishDir "${params.outdir}/logs/compare_matrices",
        mode: 'copy',
        pattern: '*.log'

    input:
    tuple val(group1),
        val(group2),
        val(region),
        val(resolution),
        path(matrix1),
        path(matrix2),
        path(chromSizes)

    val settings

    output:
    tuple val(group1),
        val(group2),
        val(region),
        val(resolution),
        path("${group2.id}-minus-${group1.id}-${region.name}-${resolution}-LOESSdiff-noFilter-full.h5"),
        path("${group2.id}-minus-${group1.id}-${region.name}-${resolution}-LOESSdiff-medianFilter-full.h5"),
        path("${group2.id}-minus-${group1.id}-${region.name}-${resolution}-subtraction.tsv"),
        emit: differences

    tuple val(group1),
        val(group2),
        val(region),
        val(resolution),
        path("${group1.id}-vs-${group2.id}-${region.name}-${resolution}-adjIF1-full.h5"),
        path("${group1.id}-vs-${group2.id}-${region.name}-${resolution}-adjIF2-full.h5"),
        emit: adjusted

    path "*.sutm", emit: sutm
    path "*.2d.txt", emit: text
    path "*.png", emit: qc
    path "${group1.id}-vs-${group2.id}-${region.name}-${resolution}.compare_matrices.log",
        emit: logs

    script:
    def comparisonPrefix = "${group1.id}-vs-${group2.id}-${region.name}-${resolution}"
    def subtractionPrefix = "${group2.id}-minus-${group1.id}-${region.name}-${resolution}"
    def regionStart = region.start + 1

    """
    set -euo pipefail

    exec > ${comparisonPrefix}.compare_matrices.log 2>&1

    python /opt/hicflow/workflow/scripts/H5_to_SUTM.py \
        ${matrix1} \
        > ${group1.id}.sutm

    python /opt/hicflow/workflow/scripts/H5_to_SUTM.py \
        ${matrix2} \
        > ${group2.id}.sutm

    mkdir -p hiccompare qc

    Rscript /opt/hicflow/workflow/scripts/HiCcompare.R \
        hiccompare \
        qc \
        '${region.chromosome}' \
        ${regionStart} \
        ${region.end} \
        '${region.name}' \
        ${resolution} \
        full \
        ${group1.id}.sutm \
        ${group2.id}.sutm

    cp hiccompare/*.sutm .
    cp hiccompare/*.2d.txt .
    cp qc/*.png .

    grep -w '${region.chromosome}' \
        ${chromSizes} \
        > ${region.name}.chrom.sizes

    printf '%s\\t%s\\t%s\\n' \
        '${region.chromosome}' \
        '${regionStart}' \
        '${region.end}' \
        > ${region.name}.bed

    for adjusted in adjIF1 adjIF2
    do
        hicConvertFormat \
            --matrices hiccompare/${comparisonPrefix}-\${adjusted}-full.2d.txt \
            --outFileName ${comparisonPrefix}-\${adjusted}-full-allChrom.h5 \
            --inputFormat 2D-text \
            --outputFormat h5 \
            --resolutions ${resolution} \
            --chromosomeSizes ${region.name}.chrom.sizes

        hicAdjustMatrix \
            --matrix ${comparisonPrefix}-\${adjusted}-full-allChrom.h5 \
            --outFileName ${comparisonPrefix}-\${adjusted}-full.h5 \
            --regions ${region.name}.bed

        hicTransform \
            --matrix ${comparisonPrefix}-\${adjusted}-full.h5 \
            --method obs_exp \
            --outFileName ${comparisonPrefix}-\${adjusted}-obsExp-full.h5
    done

    python /opt/hicflow/workflow/scripts/compareHiC.py \
        ${comparisonPrefix}-adjIF1-obsExp-full.h5 \
        ${comparisonPrefix}-adjIF2-obsExp-full.h5 \
        --outMatrix ${subtractionPrefix}-LOESSdiff-noFilter-full.h5 \
        --outMatrixFilter ${subtractionPrefix}-LOESSdiff-medianFilter-full.h5 \
        --minSum ${settings.minSum} \
        --raw ${matrix1} ${matrix2}

    printf '%s\\t%s\\n' \
        'minuend' '${group2.id}' \
        'subtrahend' '${group1.id}' \
        'operation' '${group2.id} - ${group1.id}' \
        'positive_values' 'Higher in ${group2.id}' \
        'negative_values' 'Higher in ${group1.id}' \
        'matrix_1' '${matrix1}' \
        'matrix_2' '${matrix2}' \
        > ${subtractionPrefix}-subtraction.tsv
    """
}
