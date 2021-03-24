#!/bin/bash
# shellcheck disable=SC2154
set -e -o -x

main() {
    wget -nv http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig
    chmod +x bedGraphToBigWig
    mv bedGraphToBigWig /usr/bin/

    echo "Value of is_sorted: '$is_sorted'"

    dx-download-all-inputs --parallel
    base_name=$(basename "$bedgraph_file_path")
    out_big_wig="$base_name.bw"
    bg2bw_command="bedGraphToBigWig $bedgraph_file_path $genome_sizes_file_path $out_big_wig"
    eval "$bg2bw_command"

    bigwig=$(dx upload "$out_big_wig" --brief)
    dx tag "$bigwig" sjcp-result-file
    dx-jobutil-add-output bigwig "$bigwig" --class=file
}
