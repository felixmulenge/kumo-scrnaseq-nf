#!/usr/bin/env nextflow

/*
       ______   ______       ____         ____     ___    ___                 ____   ____               
           //  //          ||    \\  || //       //     //   \\   ||     || ||      ||  \\  ||   ||    
          //   //          ||     \\ || //      //     //     \\  ||     || ||      ||   || ||   ||    
         //     //___      ||     || ||  //__   ||     ||      || ||     || ||____  ||__ //  \\ //     
        //           \\    ||     || ||      \\ ||     ||      || \\    //  ||      ||  \\     ||      
       //            \\    ||     // ||      \\ \\     \\     //   \\  //   ||      ||   ||    ||      
      //_____  _____\\     ||____//  ||  ____\\  \\___  \\___//     \\//    ||____  ||   ||    ||      

                                    Elucidate.  Innovate.  Accelerate.


    Authors:
    ZS Discovery
 
    Copyright (c) ZS Discovery

 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    bayer-int/kumo-long-read-nf
 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Website: www.zs.com/Discovery
 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SCDOWNSTREAM            } from './workflows/scdownstream'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_scdownstream_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_scdownstream_pipeline'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow KUMO {

    take:
    samplesheet // channel: samplesheet read in from --input
    ch_base  // value channel: [ val(meta), path(h5ad) ]

    main:

    //
    // WORKFLOW: Run pipeline
    //
    SCDOWNSTREAM (
        samplesheet,
        ch_base
    )
    emit:
    multiqc_report = SCDOWNSTREAM.out.multiqc_report // channel: /path/to/multiqc_report.html
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input,
        params.help,
        params.help_full,
        params.show_hidden
    )

    //
    // WORKFLOW: Run main workflow
    //
    KUMO (
        PIPELINE_INITIALISATION.out.samplesheet,
        params.base_adata
            ? channel.value([[id: "base"], file(params.base_adata, checkIfExists: true)])
            : channel.value([[], []])
    )

    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
        KUMO.out.multiqc_report
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/