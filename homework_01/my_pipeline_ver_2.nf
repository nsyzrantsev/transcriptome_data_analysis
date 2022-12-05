params.SRA = "SRR000000"
params.results_dir = "results/"

log.info ""
log.info "  Q U A L I T Y   C O N T R O L  "
log.info "================================="
log.info "SRA number         : ${params.SRA}"
log.info "Results location   : ${params.results_dir}"

process DownloadFastQ {
  publishDir "${params.results_dir}"

  output:
    path "reads/*"

  script:
    """
    /content/sratoolkit.3.0.0-ubuntu64/bin/fasterq-dump ${params.SRA} -O reads/
    """
}

process QC {
  publishDir "${params.results_dir}"

  input:
    path x

  output:
    path "qc/*.html"

  script:
    """
    mkdir qc
    /content/FastQC/fastqc -o qc $x
    """
}

workflow {
  DownloadFastQ()
  QC( DownloadFastQ.out.collect() )
}