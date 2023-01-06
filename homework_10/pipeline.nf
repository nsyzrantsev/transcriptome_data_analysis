params.results_dir = "results"
SRA_list = params.SRA.split(",")
params.index = "/content/transcriptome.idx"
params.t2g = "/content/t2g.txt"
params.chemistry = "10xv2"
params.processing_script = "/content/batches_processing.py"

log.info ""
log.info "  Q U A L I T Y   C O N T R O L  "
log.info "================================="
log.info "SRA number         : ${SRA_list}"
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

// Kallisto | Bustools process, which evaluate gene expressions by the pseudoalignment and save expressions in .h5ad file
process KalistoBustools {
  publishDir "${params.results_dir}"

  input:
    path fastq_files
    val sra

  output:
    path "${sra}/counts_unfiltered/${sra}_adata.h5ad"

  script:
    """
    kb count --keep-tmp -i ${params.index} -g ${params.t2g} -o ${sra}/ -x ${params.chemistry} --h5ad -t 2 ${fastq_files}
    mv ${sra}/counts_unfiltered/adata.h5ad ${sra}/counts_unfiltered/${sra}_adata.h5ad
    """
}

// This process eliminates batch effect and
// integrates list of batches into one .h5ad file
process BatchesProcessing {
  publishDir "${params.results_dir}"

  input:
    path h5ad_files

  output:
    path "*"

  script:
    """
    /usr/local/bin/python ${params.processing_script} ${h5ad_files}
    """
}

workflow {
  data = Channel.of( SRA_list )
  DownloadFastQ( data )
  KalistoBustools( DownloadFastQ.out, data )
  BatchesProcessing( KalistoBustools.out.collect() )
}