process SCIBMETRICS {
    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
        'community.wave.seqera.io/library/scanpy_python_yaml_pip_pruned:bb90eb826524be4e' }"

    input:
    tuple val(meta), path(h5ad)
    val (embeddings)

    output:
    tuple val(meta), path("*_metrics_raw.tsv"),         emit: raw_metrics_table
    tuple val(meta), path("*_metrics_raw.svg"),         emit: raw_metrics_plot
    tuple val(meta), path("*_metrics_scaled.tsv"),      emit: scaled_metrics_table
    tuple val(meta), path("*_metrics_scaled.svg"),      emit: scaled_metrics_plot
    path "versions.yml",                                emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    Run_scibmetrics.py \\
        --input_h5ad $h5ad \\
        --embeddings $embeddings \\
        --n_jobs $task.cpus \\
        --prefix $prefix \\
        $args

    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.h5ad"
    touch "versions.yml"
    """
}
