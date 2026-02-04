include { INFERCNV           } from '../../modules/local/infercnv'

workflow COPY_NUMBER {
    take:
    ch_h5ad // channel: [ meta, h5ad, symbol_col ]

    main:
    ch_versions = channel.empty()

    gene_order   = params.gene_order ? Channel.value(file(params.gene_order, checkIfExists: true))  : Channel.empty()
    ref_groups   = params.ref_groups ? Channel.value(file(params.ref_groups, checkIfExists: true))  : Channel.empty()

    if (params.celltypist_model) {
        INFERCNV(ch_h5ad, gene_order, ref_groups)
        ch_versions = ch_versions.mix(INFERCNV.out.versions)
    }

    emit:
    versions = ch_versions // channel: [ versions.yml ]
}