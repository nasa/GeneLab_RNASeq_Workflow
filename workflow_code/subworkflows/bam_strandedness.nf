include { INFER_EXPERIMENT } from '../modules/rseqc.nf'
include { ASSESS_STRANDEDNESS } from '../modules/assess_strandedness.nf'

workflow BAM_STRANDEDNESS {
    take:
        counting_bams
        genome_bams
        genome_bed

    main:
        paired_bams = params.mode == 'microbes' ? counting_bams.map { meta, bam -> tuple(meta, bam, bam) } :
            counting_bams.map { meta, bam -> tuple(meta.id, meta, bam) }
                .join(genome_bams.map { meta, bam -> tuple(meta.id, bam) }, failOnDuplicate: true, failOnMismatch: true)
                .map { id, meta, bam, genome -> tuple(meta, bam, genome) }
        inference_bams = paired_bams.map { meta, counting, genome -> tuple(meta, genome, []) }
        INFER_EXPERIMENT(inference_bams, genome_bed)
        ASSESS_STRANDEDNESS(INFER_EXPERIMENT.out.log.map { meta, report -> report }.collect())
        strand = ASSESS_STRANDEDNESS.out.map { f -> f.text.split(':')[0] }.first()

    emit:
        counting_bams = paired_bams.map { meta, counting, genome -> tuple(meta, counting) }
        strandedness = strand
        logs = INFER_EXPERIMENT.out.log
        versions = INFER_EXPERIMENT.out.versions
}
