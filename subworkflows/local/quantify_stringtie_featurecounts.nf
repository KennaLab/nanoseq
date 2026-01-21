/*
 * Transcript Discovery and Quantification with StringTie2 and FeatureCounts
 */

include { STRINGTIE2            } from '../../modules/local/stringtie2'
include { STRINGTIE_MERGE       } from '../../modules/nf-core/stringtie/merge/main'
include { SUBREAD_FEATURECOUNTS } from '../../modules/local/subread_featurecounts'
include { MULTIQC_CUSTOM_BIOTYPE as MULTIQC_CUSTOM_BIOTYPE_GENE } from '../../modules/local/multiqc_custom_biotype'
include { MULTIQC_CUSTOM_BIOTYPE as MULTIQC_CUSTOM_BIOTYPE_TRANSCRIPT } from '../../modules/local/multiqc_custom_biotype'

// MODULES
include { SUBREAD_FEATURECOUNTS as SUBREAD_FEATURECOUNTS_GENE } from '../../modules/nf-core/subread/featurecounts/main'
include { SUBREAD_FEATURECOUNTS as SUBREAD_FEATURECOUNTS_TRANSCRIPT } from '../../modules/nf-core/subread/featurecounts/main'


workflow QUANTIFY_STRINGTIE_FEATURECOUNTS {
    take:
    ch_sample     // [ sample, barcode, fasta, gtf, is_transcripts, annotation_str ]
    ch_sortbam

    main:

    ch_sample
        .map  { it -> [ it[0], it[2], it[3] ] }
        .join ( ch_sortbam )
        .set  { ch_sample } //tuple val(meta), path(fasta), path(gtf), path(bam)

    /*
     * Novel isoform detection with StringTie
     */
    STRINGTIE2 ( ch_sample )
    ch_stringtie_gtf   = STRINGTIE2.out.stringtie_gtf
    stringtie2_version = STRINGTIE2.out.versions

    ch_sample
        .map { it -> [ it[2] ] }
        .unique()
        .set { ch_sample_gtf }

    /*
     * Merge isoforms across samples called by StringTie
     */
    STRINGTIE_MERGE ( ch_stringtie_gtf.collect(), ch_sample_gtf )
    ch_stringtie_merged_gtf = STRINGTIE_MERGE.out.gtf

    /*
     * Gene and transcript quantification with featureCounts
     */
    // ch_sample
    //     .collect { it[-1]    }
    //     .set     { ch_sample }
    // SUBREAD_FEATURECOUNTS ( ch_stringtie_merged_gtf, ch_sample )


    ch_sample
        .map {meta, fasta, gtf, bam ->
            def fmeta = [:]
            // Set meta.id
            fmeta.id = meta
            // Set meta.single_end
            fmeta.single_end = true
            return [fmeta, bam, gtf]
        }
        .set { ch_bam_gtf }

    // Run subread featurecounts with different meta_features.
    SUBREAD_FEATURECOUNTS_GENE(ch_bam_gtf)
    SUBREAD_FEATURECOUNTS_TRANSCRIPT(ch_bam_gtf)

    // Custom BIOTYPE plot from featurecounts output.
    ch_biotypes_header_multiqc   = file("$projectDir/assets/multiqc/biotypes_header.txt", checkIfExists: true)

    MULTIQC_CUSTOM_BIOTYPE_GENE ( SUBREAD_FEATURECOUNTS_GENE.out.counts, ch_biotypes_header_multiqc )
    MULTIQC_CUSTOM_BIOTYPE_TRANSCRIPT ( SUBREAD_FEATURECOUNTS_TRANSCRIPT.out.counts, ch_biotypes_header_multiqc )

    emit:
    ch_stringtie_gtf
    ch_stringtie_merged_gtf
    stringtie2_version

    ch_gene_counts                           = SUBREAD_FEATURECOUNTS_GENE.out.counts
    ch_transcript_counts                     = SUBREAD_FEATURECOUNTS_TRANSCRIPT.out.counts
    featurecounts_gene_multiqc               = SUBREAD_FEATURECOUNTS_GENE.out.summary
    featurecounts_transcript_multiqc         = SUBREAD_FEATURECOUNTS_TRANSCRIPT.out.summary
    featurecounts_gene_version               = SUBREAD_FEATURECOUNTS_GENE.out.versions
    featurecounts_transcript_version         = SUBREAD_FEATURECOUNTS_TRANSCRIPT.out.versions
    featurecounts_biotype_gene_version       = MULTIQC_CUSTOM_BIOTYPE_GENE.out.versions
    featurecounts_multiqc_biotype_gene       = MULTIQC_CUSTOM_BIOTYPE_GENE.out.tsv
    featurecounts_biotype_transcript_version = MULTIQC_CUSTOM_BIOTYPE_TRANSCRIPT.out.versions
    featurecounts_multiqc_biotype_transcript = MULTIQC_CUSTOM_BIOTYPE_TRANSCRIPT.out.tsv
}
