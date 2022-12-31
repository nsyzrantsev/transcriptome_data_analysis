params.results_dir = "results/"
SRA_list = params.SRA.split(",")
params.index = "/content/transcriptome.idx"
params.r_script = "/content/dropletutils_script.R"

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
    path x

  output:
    path "kalisto_bustools/*"

  script:
    """
    mkdir kalisto_bustools
    /content/kallisto/build/src/kb count -i ${params.index} -o kalisto_bustools -x 10xv3 --h5ad $x
    """
}

// Delete droplets from .h5ad files
process DropletUtils {
  publishDir "${params.results_dir}"

  input:
    path x

  output:
    path "*.h5ad"

  script:
    """
    /usr/local/bin/Rscript ${params.r_script} ${x}
    """
}

workflow {
  data = Channel.of( SRA_list )
  DownloadFastQ( data )
  KalistoBustools( DownloadFastQ.out )
  DropletUtils( KalistoBus.out )
}
