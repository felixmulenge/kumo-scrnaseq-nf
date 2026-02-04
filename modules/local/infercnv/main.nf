process INFERCNV {
    tag "$meta.id"
    label 'process_medium'


    container "docker.io/trinityctat/infercnv:1.20.0"

    input:
    tuple val(meta), path(h5ad)
    path (gene_order_tsv)
    //path annotations_tsv
    val  ref_groups

    output:
    tuple val(meta), path("*infercnv_obj.rds"),   emit: infercnv_rds
    tuple val(meta), path("*infercnv.tsv"),       emit: metadata_table
    tuple val(meta), path("*infercnv.png"),       emit: infercnv_heatmap
    tuple val(meta), path("*seurat.rds"),         emit: seurat_rds
    path "versions.yml",                          emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
 
    """
    Run_infercnv.R \\
        --input_h5ad $h5ad \\
        --ref_groups $ref_groups \\
        --gene_order_file $gene_order_tsv \\
        --threads $task.cpus \\
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


