include { SCIBMETRICS           } from '../../modules/local/scibmetrics'

workflow SCIB_METRICS {
    take:
    ch_h5ad // channel: [ integration, h5ad ]
    ch_embeddings

    main:
    ch_versions = channel.empty()

    SCIBMETRICS(ch_h5ad, ch_embeddings)
    ch_versions = ch_versions.mix(SCIBMETRICS.out.versions)

    emit:
    versions = ch_versions // channel: [ versions.yml ]
}