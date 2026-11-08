process PLOT_HIC {
    tag "${meta.id}:${coordinate.id}:${resolution}"

    label 'process_medium'

    publishDir "${params.outdir}/plots/${region.name}/${resolution}/pyGenomeTracks",
        mode: 'copy',
        pattern: "*.${settings.filetype}"

    publishDir "${params.outdir}/plots/${region.name}/${resolution}/pyGenomeTracks/configs",
        mode: 'copy',
        pattern: '*.ini'

    publishDir "${params.outdir}/logs/plot_hic",
        mode: 'copy',
        pattern: '*.log'

    input:
    tuple val(meta),
        val(region),
        val(resolution),
        path(matrix),
        path(insulation),
        path(tads),
        val(coordinate)

    path genes
    path bigWigs
    val settings

    output:
    tuple val(meta),
        val(region),
        val(resolution),
        val(coordinate),
        path("${meta.id}-${region.name}-${coordinate.id}-${resolution}-full.${settings.filetype}"),
        emit: plots

    path "${meta.id}-${region.name}-${coordinate.id}-${resolution}-full.ini",
        emit: configs

    path "${meta.id}-${region.name}-${coordinate.id}-${resolution}.plot_hic.log",
        emit: logs

    script:
    def prefix = "${meta.id}-${region.name}-${coordinate.id}-${resolution}"
    def depth = coordinate.end - coordinate.start
    def title = "${meta.id} : ${region.name} (${settings.genomeBuild}) at ${resolution} bin size (KR - full)"

    def bigWigArguments = settings.bigWigTitles
        .withIndex()
        .collect { bigWigTitle, index ->
            "--bigWig '${bigWigTitle},${bigWigs[index]},3'"
        }
        .join(' ')

    def vLinesArgument = settings.vlines
        ? "--vLines '${settings.vlines}'"
        : ''

    """
    set -euo pipefail

    export NUMEXPR_MAX_THREADS=1

    python /opt/hicflow/workflow/scripts/generate_config.py \
        --matrix ${matrix} \
        --log \
        --colourmap '${settings.colourmap}' \
        --depth ${depth} \
        --geneLabelFontSize ${settings.geneLabelFontSize} \
        --maxLabels ${settings.maxLabels} \
        --insulation ${insulation} \
        --tads ${tads} \
        --genes 'Genes,${genes},3' \
        ${bigWigArguments} \
        ${vLinesArgument} \
        > ${prefix}-full.ini \
        2> ${prefix}.plot_hic.log

    pyGenomeTracks \
        --tracks ${prefix}-full.ini \
        --region '${coordinate.chromosome}:${coordinate.start + 1}-${coordinate.end}' \
        --outFileName ${prefix}-full.${settings.filetype} \
        --title '${title}' \
        --dpi ${settings.dpi} \
        >> ${prefix}.plot_hic.log \
        2>&1
    """
}
