// TODO nf-core: If in doubt look at other nf-core/modules to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/modules/nf-core/
//               You can also ask for help via your pull request or on the #modules channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A module file SHOULD only define input and output files as command-line parameters.
//               All other parameters MUST be provided using the "task.ext" directive, see here:
//               https://www.nextflow.io/docs/latest/process.html#ext
//               where "task.ext" is a string.
//               Any parameters that need to be evaluated in the context of a particular sample
//               e.g. single-end/paired-end data MUST also be defined and evaluated appropriately.
// TODO nf-core: Software that can be piped together SHOULD be added to separate module files
//               unless there is a run-time, storage advantage in implementing in this way
//               e.g. it's ok to have a single module for bwa to output BAM instead of SAM:
//                 bwa mem | samtools view -B -T ref.fasta
// TODO nf-core: Optional inputs are not currently supported by Nextflow. However, using an empty
//               list (`[]`) instead of a file can be used to work around this issue.

process INFERCNV {
    tag "$meta.id"
    label 'process_medium'

    // TODO nf-core: See section in main README for further information regarding finding and adding container addresses to the section below.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
        'biocontainers/YOUR-TOOL-HERE' }"

    input:// TODO nf-core: Where applicable all sample-specific information e.g. "id", "single_end", "read_group"
    //               MUST be provided as an input via a Groovy Map called "meta".
    //               This information may not be required in some instances e.g. indexing reference genome files:
    //               https://github.com/nf-core/modules/blob/master/modules/nf-core/bwa/index/main.nf
    // TODO nf-core: Where applicable please provide/convert compressed files as input/output
    //               e.g. "*.fastq.gz" and NOT "*.fastq", "*.bam" and NOT "*.sam" etc.
    tuple val(meta), path(bam)

    output:
    // TODO nf-core: Named file extensions MUST be emitted for ALL output channels
    tuple val(meta), path("*.bam"), emit: bam
    // TODO nf-core: List additional required output channels/values here
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // TODO nf-core: Where possible, a command MUST be provided to obtain the version number of the software e.g. 1.10
    //               If the software is unable to output a version number on the command-line then it can be manually specified
    //               e.g. https://github.com/nf-core/modules/blob/master/modules/nf-core/homer/annotatepeaks/main.nf
    //               Each software used MUST provide the software name and version number in the YAML version file (versions.yml)
    // TODO nf-core: It MUST be possible to pass additional parameters to the tool as a command-line string via the "task.ext.args" directive
    // TODO nf-core: If the tool supports multi-threading then you MUST provide the appropriate parameter
    //               using the Nextflow "task" variable e.g. "--threads $task.cpus"
    // TODO nf-core: Please replace the example samtools command below with your module's command
    // TODO nf-core: Please indent the command appropriately (4 spaces!!) to help with readability ;)
    """
    infercnv \\
        $args \\
        -@ $task.cpus \\
        -o ${prefix}.bam \\
        $bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        infercnv: \$(infercnv --version)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // TODO nf-core: A stub section should mimic the execution of the original module as best as possible
    //               Have a look at the following examples:
    //               Simple example: https://github.com/nf-core/modules/blob/818474a292b4860ae8ff88e149fbcda68814114d/modules/nf-core/bcftools/annotate/main.nf#L47-L63
    //               Complex example: https://github.com/nf-core/modules/blob/818474a292b4860ae8ff88e149fbcda68814114d/modules/nf-core/bedtools/split/main.nf#L38-L54
    // TODO nf-core: If the module doesn't use arguments ($args), you SHOULD remove:
    //               - The definition of args `def args = task.ext.args ?: ''` above.
    //               - The use of the variable in the script `echo $args ` below.
    """
    echo $args
    
    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        infercnv: \$(infercnv --version)
    END_VERSIONS
    """
}


process INFERCNV_FROM_H5AD {
  tag "${sample_id}"

  container "docker.io/trinityctat/infercnv:latest"

  /*
   * CONTRACT / INTERFACE
   * --------------------
   * Required Outputs:
   *  - infercnv_out/ (directory)
   *      Expected to include inferCNV outputs AND (if write_feature_table=true):
   *        infercnv_out/map_metadata_from_infercnv.txt
   *  - seurat_with_infercnv_features.rds
   *
   * REQUIRED DEPENDENCIES:
   *  - R packages: infercnv, Seurat, zellkonverter, SingleCellExperiment, SummarizedExperiment, Matrix, optparse
   */

  publishDir "${params.outdir ?: 'results'}/infercnv/${sample_id}", mode: 'copy'

  input:
    tuple val(sample_id), path(h5ad)
    path annotations_tsv
    path gene_order_tsv

  output:
    tuple val(sample_id), path("infercnv_out"), emit: infercnv_dir
    tuple val(sample_id), path("seurat_with_infercnv_features.rds"), emit: seurat_rds

  script:
    def ref_groups          = params.infercnv_ref_groups ?: ""
    def cutoff              = params.infercnv_cutoff ?: 0.10
    def cluster_by_groups   = (params.infercnv_cluster_by_groups != null) ? params.infercnv_cluster_by_groups : true
    def denoise             = (params.infercnv_denoise != null) ? params.infercnv_denoise : true
    def analysis_mode       = params.infercnv_analysis_mode ?: "subclusters"
    def num_threads         = params.infercnv_num_threads ?: task.cpus
    def tumor_method        = params.infercnv_tumor_subcluster_partition_method ?: "leiden"

    def feature_top_n       = params.infercnv_feature_top_n ?: 10
    def write_feature_table = (params.infercnv_write_feature_table != null) ? params.infercnv_write_feature_table : true

    """
    Rscript ${projectDir}/bin/infercnv_from_h5ad.R \
      --h5ad ${h5ad} \
      --annotations ${annotations_tsv} \
      --gene_order ${gene_order_tsv} \
      --ref_groups '${ref_groups}' \
      --outdir infercnv_out \
      --cutoff ${cutoff} \
      --cluster_by_groups ${cluster_by_groups} \
      --denoise ${denoise} \
      --HMM true \
      --analysis_mode ${analysis_mode} \
      --num_threads ${num_threads} \
      --tumor_subcluster_partition_method ${tumor_method} \
      --extract_features true \
      --feature_top_n ${feature_top_n} \
      --write_feature_table ${write_feature_table}
    """
}







process INFERCNV_FROM_H5AD {
    tag "${meta.id}"
    label 'process_medium'

    container "docker.io/YOUR_ORG/infercnv:YOUR_TAG"

    input:
    tuple val(meta), path(h5ad)
    path annotations_tsv
    path gene_order_tsv
    val(ref_groups)
    val(cutoff)
    val(cluster_by_groups)
    val(denoise)
    val(analysis_mode)
    val(feature_top_n)
    val(write_feature_table)

    output:
    tuple val(meta), path("${prefix}.infercnv_out"), emit: infercnv_dir
    tuple val(meta), path("${prefix}.seurat_with_infercnv_features.rds"), emit: seurat_rds
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template 'infercnv_from_h5ad.R'

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}.infercnv_out
    touch ${prefix}.seurat_with_infercnv_features.rds
    cat > versions.yml <<-END
    INFERCNV_FROM_H5AD:
      r-base: "stub"
      infercnv: "stub"
      seurat: "stub"
      zellkonverter: "stub"
    END
    """
}

