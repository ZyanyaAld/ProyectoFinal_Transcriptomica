nextflow.enable.dsl=2

/*
 * Parámetros
 */

params.input  = "${projectDir}/../data/raw/*_{1,2}.fastq.gz"
params.outdir = "${projectDir}/../results/nextflow"


/*
 * FASTQC RAW
 */

process FASTQC_RAW {

    tag "$sample_id"

    publishDir "${params.outdir}/qc/fastqc_raw", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    path "*_fastqc.html", emit: html
    path "*_fastqc.zip",  emit: zip

    script:
    """
    conda run -n base fastqc \
        ${reads[0]} \
        ${reads[1]} \
        -o .
    """
}



/*
 * FASTP
 */

process FASTP {

    tag "$sample_id"

    publishDir "${params.outdir}/fastp", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_*.clean.fastq.gz")

    script:
    """
    conda run -n fastp fastp \
        -i ${reads[0]} \
        -I ${reads[1]} \
        -o ${sample_id}_1.clean.fastq.gz \
        -O ${sample_id}_2.clean.fastq.gz \
        -h ${sample_id}.html \
        -j ${sample_id}.json
    """
}



/*
 * FASTQC CLEAN
 */

process FASTQC_CLEAN {

    tag "$sample_id"

    publishDir "${params.outdir}/qc/fastqc_clean", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    path "*_fastqc.html", emit: html
    path "*_fastqc.zip",  emit: zip

    script:
    """
    conda run -n base fastqc \
        ${reads[0]} \
        ${reads[1]} \
        -o .
    """
}



/*
 * MULTIQC RAW
 */

process MULTIQC_RAW {

    publishDir "${params.outdir}/multiqc_raw", mode: 'copy'

    input:
    path(qc_files)

    output:
    path "multiqc_raw_report.html"

    script:
    """
    mkdir raw_qc

    cp ${qc_files} raw_qc/

    conda run -n base multiqc raw_qc -o .

    mv multiqc_report.html multiqc_raw_report.html
    """
}



/*
 * MULTIQC CLEAN
 */

process MULTIQC_CLEAN {

    publishDir "${params.outdir}/multiqc_clean", mode: 'copy'

    input:
    path(qc_files)

    output:
    path "multiqc_clean_report.html"

    script:
    """
    mkdir clean_qc

    cp ${qc_files} clean_qc/

    conda run -n base multiqc clean_qc -o .

    mv multiqc_report.html multiqc_clean_report.html
    """
}



/*
 * WORKFLOW PRINCIPAL
 */

workflow {

    /*
     * Detectar pares R1/R2
     */

    Channel
        .fromFilePairs(params.input, checkIfExists: true)
        .set { fastq_files }


    /*
     * FastQC inicial
     */

    raw_qc = FASTQC_RAW(fastq_files)


    /*
     * MultiQC RAW
     */

    raw_qc.zip
        .collect()
        .set { raw_multiqc_input }

    MULTIQC_RAW(raw_multiqc_input)


    /*
     * Limpieza con fastp
     */

    cleaned_reads = FASTP(fastq_files)


    /*
     * FastQC después de limpieza
     */

    clean_qc = FASTQC_CLEAN(cleaned_reads)


    /*
     * MultiQC CLEAN
     */

    clean_qc.zip
        .collect()
        .set { clean_multiqc_input }

    MULTIQC_CLEAN(clean_multiqc_input)
}