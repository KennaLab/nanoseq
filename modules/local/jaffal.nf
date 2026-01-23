process JAFFAL {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::jaffa=2.3.0"
    container "docker.io/wdesaint/jaffal_older_bpipe:latest"
    // container "file:///hpc/hers_en/edejong2/software/singularity_cache/ghcr.io-kennalab-jaffa-2.4.sif"
    // container "docker.io/davidsongroup/jaffa:2.4"
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/jaffa:2.3--hdfd78af_0' :
    //     'quay.io/biocontainers/jaffa:2.3--hdfd78af_0' }"

    input:
    tuple val(meta), path(fastq)
    path jaffal_ref_dir

    output:
    tuple val(meta), path("*.fasta"), emit: jaffal_fastq
    path "*.csv"                    , emit: jaffal_results
    path "*jaffa.out"               , emit: jaffal_stdout
    path "*jaffa.err"               , emit: jaffal_stderr
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    bpipe run -p refBase=$jaffal_ref_dir $jaffal_ref_dir/JAFFAL.groovy $fastq > ${prefix}.jaffa.out 2> ${prefix}.jaffa.err
    mv jaffa_results.csv ${prefix}.jaffa_results.csv
    mv jaffa_results.fasta ${prefix}.jaffa_results.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        jaffa: \$( echo 'jaffa 2.1' )
    END_VERSIONS
    """
}
