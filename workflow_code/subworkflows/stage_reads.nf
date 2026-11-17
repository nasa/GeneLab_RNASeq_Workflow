include { COPY_READS } from '../modules/copy_reads.nf'
include { FETCH_READS } from '../modules/fetch_remote.nf'

workflow STAGE_READS {
    take:
        ch_samples
        type

    main:
        ch_split = ch_samples.branch { _meta, files ->
            remote: files && files[0]?.toString()?.contains("://")
            local: true
        }

        FETCH_READS( ch_split.remote, type )
        COPY_READS( ch_split.local, type )
        ch_reads = FETCH_READS.out.reads.mix( COPY_READS.out.reads )

        ch_reads | map { sample -> sample[1] } | collect | set { ch_all_reads }
        ch_reads | map { sample -> sample[0].id }
                        | collectFile(name: "samples.txt", sort: true, newLine: true)
                        | set { samples_txt }

    emit:
        reads = ch_reads
        ch_all_reads = ch_all_reads
        samples_txt = samples_txt
}
