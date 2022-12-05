params.results_dir = "/content/results"
SRA_list = params.SRA.split(",")
params.index = "/content/results/index/index.idx"
params.layout = "SINGLE"

log.info ""
log.info "  Q U A L I T Y   C O N T R O L  "
log.info "================================="
log.info "SRA number         : ${SRA_list}"
log.info "Index location     : ${params.index}"
log.info "Layout type        : ${params.layout}"
log.info "Results location   : ${params.results_dir}"

// Downloading .fastq files by the fasterq-dump tool of the sratoolkit
process DownloadFastQ {
  publishDir "${params.results_dir}"

  input:
    val sra

  output:
    path "${sra}/*"

  script:
    """
    /content/sratoolkit.3.0.0-ubuntu64/bin/fasterq-dump ${sra} -O ${sra}/
    """
}

// Creating quality control by the FastQC
process FastQC {
  publishDir "${params.results_dir}/${sra}"

  input:
    val sra
    path fastq_path

  output:
    path "fastqc/*"

  script:
    """
    mkdir fastqc
    /content/FastQC/fastqc -o fastqc $fastq_path
    """
}

// Creating quality control by the MultiQC
process MultiQC {
  publishDir "${params.results_dir}/${sra}"

  input:
    val sra
    path fastqc_data

  output:
    path "multiqc_report.html"

  script:
    """
    multiqc ${fastqc_data}
    """
}

// Kallisto quant algorithm, which quantifies abundances of transcripts
process KallistoQuant {
  publishDir "${params.results_dir}/${sra}"

  input:
    val layout_type
    path fastq_file

  output:
    path "kallisto/*"

  script:
  if( layout_type == 'SINGLE' )
    """
    /content/kallisto/build/src/kallisto quant -i ${params.index} -o kallisto/ --single -l 180 -s 20 ${fastq_file}
    """
  else if( layout_type == 'PAIRED' )
    """
    /content/kallisto/build/src/kallisto quant -i ${params.index} -o kallisto/ ${fastq_file}
    """
}

workflow {
  data = Channel.of( SRA_list )
  DownloadFastQ( data )
  FastQC( data, DownloadFastQ.out )
  MultiQC( data, FastQC.out.collect() )
  layout_type = Channel.of( params.layout )
  KallistoQuant( layout_type, DownloadFastQ.out )
}