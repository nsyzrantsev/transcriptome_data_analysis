params.results_dir = "/content/results"
SRA_list = params.SRA.split(",")
params.fasta = "ftp.ensembl.org/pub/release-108/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz"
params.layout = "SINGLE"

log.info ""
log.info "  Q U A L I T Y   C O N T R O L  "
log.info "================================="
log.info "SRA number         : ${SRA_list}"
log.info "Fasta URL          : ${params.fasta}"
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
    path "fastqc/*.html"

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

process DownloadFastaArchive {
  output:
    path "*"
  
  script:
    """
    curl -O ${params.fasta}
    """
}

// Generating index file for kallisto quant algorithm 
process KallistoIndex {
  publishDir "${params.results_dir}"

  input:
    path fasta_archive

  output:
    path "kallisto/human_transcriptome_reference.idx"

  script:
    """
    mkdir kallisto
    /content/kallisto/build/src/kallisto index -i kallisto/human_transcriptome_reference.idx ${fasta_archive}
    """
}

// Kallisto quant algorithm, which quantifies abundances of transcripts
process KallistoQuant {
  publishDir "${params.results_dir}/kallisto"

  input:
    val sra
    val layout_type
    path index_file
    path fastq_file

  output:
    path "${sra}/*"

  script:
  if( layout_type == 'SINGLE' )
    """
    mkdir ${sra}
    /content/kallisto/build/src/kallisto quant -i ${index_file} -o ${sra} --single -l 180 -s 20 ${fastq_file}
    """
  else if( layout_type == 'PAIRED' )
    """
    mkdir ${sra}
    /content/kallisto/build/src/kallisto quant -i ${index_file} -o ${sra} ${fastq_file}
    """
}

workflow {
  data = Channel.of( SRA_list )
  DownloadFastQ( data )
  FastQC( data, DownloadFastQ.out )
  MultiQC( data, FastQC.out.collect() )
  DownloadFastaArchive()
  KallistoIndex( DownloadFastaArchive.out )
  layout_type = Channel.of( params.layout )
  KallistoQuant( data, params.layout, KallistoIndex.out, DownloadFastQ.out )
}
