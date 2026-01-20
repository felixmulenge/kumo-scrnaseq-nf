// 
//  .......  .....     ######   #
//       .  .     .    #     #  
//      .   .          #     #  #   ###    ###   ###   #   #   ###   # ###  #   #
//     .     .....     #     #  #  #      #     #   #  #   #  #   #  ##     #   #
//    .           .    #     #  #   ###   #     #   #  #   #  #####  #      #   #
//   .      .     .    #     #  #      #  #     #   #   # #   #      #       # #
//  .......  .....     ######   #   ###    ###   ###     #     ###   #        #
//                                                                           #
//                     Elucidate.  Innovate.  Accelerate.
//
//  Authors:
//  ZS Discovery
//
//  Copyright (c) ZS Discovery
//
//  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//    bayer-int/kumo-long-read-nf
//  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//    Website: www.zs.com/Discovery
//  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


include { LOAD_H5AD                            } from '../subworkflows/local/load_h5ad'
include { QUALITY_CONTROL                      } from '../subworkflows/local/quality_control'
include { CELLTYPE_ASSIGNMENT                  } from '../subworkflows/local/celltype_assignment'
include { ADATA_EXTEND as FINALIZE_QC_ANNDATAS } from '../modules/local/adata/extend'
include { COMBINE                              } from '../subworkflows/local/combine'
include { ADATA_SPLITEMBEDDINGS                } from '../modules/local/adata/splitembeddings'
include { CLUSTER                              } from '../subworkflows/local/cluster'
include { PSEUDOBULKING                        } from '../subworkflows/local/pseudobulking'
include { PER_GROUP                            } from '../subworkflows/local/per_group'
include { FINALIZE                             } from '../subworkflows/local/finalize'
include { MULTIQC                              } from '../modules/nf-core/multiqc/main'
include { softwareVersionsToYAML               } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { samplesheetToList                    } from 'plugin/nf-schema'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowScdownstream.initialise(params, log)

//
// Custom validation for pipeline parameters
//
validateInputParameters()

// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }


// Create channel from input file provided through params.input
//
ch_samplesheet = params.input ? channel.fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json")).map {
        sample -> validateInputSamplesheet(sample)
    }
    : channel.empty()

ch_base = params.base_adata ? Channel.value([[id: "base"], file(params.base_adata, checkIfExists: true)]) : Channel.value([[], []])
   


// Info required for completion email and summary
def multiqc_report = []


workflow SCDOWNSTREAM {


    ch_versions = channel.empty()
    ch_integrations = channel.empty()
    ch_obs = channel.empty()
    ch_var = channel.empty()
    ch_obsm = channel.empty()
    ch_obsp = channel.empty()
    ch_uns = channel.empty()
    ch_layers = channel.empty()
    ch_multiqc_files = channel.empty()

    if (params.input) {
        ch_obs_per_sample = channel.empty()
        ch_var_per_sample = channel.empty()
        ch_obsm_per_sample = channel.empty()
        ch_obsp_per_sample = channel.empty()
        ch_uns_per_sample = channel.empty()
        ch_layers_per_sample = channel.empty()

        //
        // Load/Convert input to h5ad
        //
        LOAD_H5AD(ch_samplesheet)
        ch_h5ad = LOAD_H5AD.out.h5ad
        ch_versions = ch_versions.mix(LOAD_H5AD.out.versions)

        //
        // Quality control per sample
        //
        QUALITY_CONTROL(
            ch_h5ad,
            params.ambient_correction,
            (!params.doublet_detection || params.doublet_detection == 'none') ? [] : params.doublet_detection.split(',').collect { it -> it.trim().toLowerCase() },
            params.doublet_detection_threshold,
            params.mito_genes,
        )
        ch_versions = ch_versions.mix(QUALITY_CONTROL.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(QUALITY_CONTROL.out.multiqc_files)
        ch_h5ad = QUALITY_CONTROL.out.h5ad

        //
        // Perform automated celltype assignment
        //
        CELLTYPE_ASSIGNMENT(ch_h5ad.map { meta, h5ad -> [meta, h5ad, meta.symbol_col] })
        ch_versions = ch_versions.mix(CELLTYPE_ASSIGNMENT.out.versions)
        ch_obs_per_sample = ch_obs_per_sample.mix(CELLTYPE_ASSIGNMENT.out.obs)

        FINALIZE_QC_ANNDATAS(
            ch_h5ad.join(ch_obs_per_sample.groupTuple(), remainder: true).join(ch_var_per_sample.groupTuple(), remainder: true).join(ch_obsm_per_sample.groupTuple(), remainder: true).join(ch_obsp_per_sample.groupTuple(), remainder: true).join(ch_uns_per_sample.groupTuple(), remainder: true).join(ch_layers_per_sample.groupTuple(), remainder: true).map { meta, h5ad, obs, var, obsm, obsp, uns, layers ->
                [meta, h5ad, obs ?: [], var ?: [], obsm ?: [], obsp ?: [], uns ?: [], layers ?: []]
            }
        )
        ch_h5ad = FINALIZE_QC_ANNDATAS.out.h5ad
        ch_versions = ch_versions.mix(FINALIZE_QC_ANNDATAS.out.versions)

        if (!params.qc_only) {
            //
            // Combine samples and perform integration
            //
            COMBINE(ch_h5ad, ch_base)
            ch_versions = ch_versions.mix(COMBINE.out.versions)
            ch_obs = ch_obs.mix(COMBINE.out.obs)
            ch_var = ch_var.mix(COMBINE.out.var)
            ch_obsm = ch_obsm.mix(COMBINE.out.obsm)
            ch_integrations = ch_integrations.mix(COMBINE.out.integrations)
            ch_finalization_base = COMBINE.out.h5ad

            ch_label_grouping = COMBINE.out.h5ad_inner
            grouping_col = "label"
        }
    }
    else {
        ch_embeddings = channel.value(params.base_embeddings.split(',').collect { it -> it.trim() })

        ADATA_SPLITEMBEDDINGS(ch_base, ch_embeddings)
        ch_versions = ch_versions.mix(ADATA_SPLITEMBEDDINGS.out.versions)
        ch_integrations = ch_integrations.mix(
            ADATA_SPLITEMBEDDINGS.out.h5ad.map { _meta, h5ads -> h5ads }.flatten().map { h5ad -> [[id: h5ad.simpleName, integration: h5ad.simpleName], h5ad] }
        )

        ch_finalization_base = ch_base
        ch_label_grouping = ch_base
        grouping_col = params.base_label_col
    }

    //
    // Perform clustering and per-cluster analysis
    //
    if (!params.qc_only) {
        CLUSTER(
            ch_integrations,
            params.cluster_per_label,
            params.cluster_global,
            params.input ? "label" : params.base_label_col,
            params.clustering_resolutions.split(','),
            "batch",
            "X_emb",
        )
        ch_versions = ch_versions.mix(CLUSTER.out.versions)
        ch_obs = ch_obs.mix(CLUSTER.out.obs)
        ch_obsm = ch_obsm.mix(CLUSTER.out.obsm)
        ch_multiqc_files = ch_multiqc_files.mix(CLUSTER.out.multiqc_files)

        if (params.pseudobulk) {
            PSEUDOBULKING(
                CLUSTER.out.h5ad_clustering,
                params.pseudobulk_groupby_labels.split(','),
                params.pseudobulk_min_num_cells,
                "X",
            )
            ch_versions = ch_versions.mix(PSEUDOBULKING.out.versions)
        }

        PER_GROUP(
            CLUSTER.out.h5ad_clustering.map { meta, h5ad -> [meta + [obs_key: "${meta.id}_leiden"], h5ad] },
            CLUSTER.out.h5ad_neighbors.map { meta, h5ad -> [meta + [obs_key: grouping_col], h5ad] },
            ch_label_grouping.map { meta, h5ad -> [meta + [obs_key: grouping_col], h5ad] },
        )
        ch_versions = ch_versions.mix(PER_GROUP.out.versions)
        ch_uns = ch_uns.mix(PER_GROUP.out.uns)
        ch_multiqc_files = ch_multiqc_files.mix(PER_GROUP.out.multiqc_files)

        FINALIZE(ch_finalization_base, ch_obs, ch_var, ch_obsm, ch_obsp, ch_uns, ch_layers)
        ch_versions = ch_versions.mix(FINALIZE.out.versions)
    }

    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_' + 'scdownstream_software_' + 'mqc_' + 'versions.yml',
            sort: true,
            newLine: true,
        )
        .set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()
ch_multiqc_logo          = params.multiqc_logo   ? Channel.fromPath(params.multiqc_logo, checkIfExists: true) : Channel.empty()


workflow_summary    = WorkflowScdownstream.paramsSummaryMultiqc(workflow, summary_params)
ch_workflow_summary = Channel.value(workflow_summary)

ch_multiqc_files = Channel.empty()
ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)


    MULTIQC(
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        [],
    )


multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
ch_versions    = ch_versions.mix(MULTIQC.out.versions)

}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)




}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
//
// Check and validate pipeline parameters
//
def validateInputParameters() {
    if (!params.input && !(params.base_adata && params.base_embeddings && params.base_label_col)) {
        throw new Exception("Either an input samplesheet or (base_adata && base_embeddings && base_label_col) must be provided")
    }

    if (params.qc_only && !params.input) {
        throw new Exception("If qc_only is set to true, an input samplesheet must be provided")
    }

    def integration_methods = params.integration_methods.split(',').collect { it.trim().toLowerCase() }
    if (params.input && params.base_adata && (integration_methods - ['scvi', 'scanvi', 'scimilarity']).size() > 0) {
        throw new Exception("Only scvi, scanvi and scimilarity integration methods are supported if base_adata is provided")
    }

    if (params.base_adata && 'scvi' in integration_methods && !params.scvi_model) {
        throw new Exception("If base_adata is provided and scvi is used as integration method, scvi_model must be provided.")
    }

    if (params.base_adata && 'scanvi' in integration_methods && !params.scanvi_model) {
        throw new Exception("If base_adata is provided and scanvi is used as integration method, scanvi_model must be provided.")
    }

    if (params.base_adata && 'scimilarity' in integration_methods && !params.scimilarity_model) {
        throw new Exception("If base_adata is provided and scimilarity is used as integration method, scimilarity_model must be provided.")
    }
}

// Validate channels from input samplesheet
//
def validateInputSamplesheet(input) {
    def (meta, filtered, unfiltered) = input
    if (!filtered && !unfiltered) {
        throw new Exception("Both filtered and unfiltered files are missing for sample ${meta.id}")
    }

    return input
}
