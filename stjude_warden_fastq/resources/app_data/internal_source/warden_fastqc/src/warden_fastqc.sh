#!/bin/bash
# shellcheck disable=SC2154
set -e -x -o pipefail

main() {
    wget -nv https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip
    unzip -q fastqc_v0.11.9.zip
    fastqc=FastQC/fastqc
    chmod +x $fastqc

    dx-download-all-inputs --parallel

    out_dir=out/out_files
    mkdir -p $out_dir
    basename="${fastq_input_name%.fa*}"
    basename="${basename%.fq*}"

    fastqc_comamnd="$fastqc --noextract  $fastq_input_path"
    eval "$fastqc_comamnd"

    html_file=$(dx upload /home/dnanexus/in/fastq_input/*fastqc.html --brief)
    dx tag "$html_file" sjcp-result-file
    dx-jobutil-add-output html_file "$html_file" --class=file

    zip_file=$(dx upload /home/dnanexus/in/fastq_input/*fastqc.zip --brief)

    dx tag "$zip_file" sjcp-result-file
    dx-jobutil-add-output zip_file "$zip_file" --class=file
}
