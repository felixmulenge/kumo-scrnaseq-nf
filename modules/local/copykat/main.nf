process COPYKAT {
    tag "$meta.id"
    label 'process_medium'


    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
        'biocontainers/YOUR-TOOL-HERE' }"

    input:
    tuple val(meta), path(h5ad)
    //path genome

    output:
    tuple val(meta), path("${prefix}.h5ad")     , emit: h5ad
    path "versions.yml"                         , emit: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    Run_copykat.R 
        --input_h5ad $h5ad \\
        --prefix $prefix \\
        --n_cores $task.cpus \\
        $args 

    """

   stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.h5ad"
    touch "versions.yml"
    """
}