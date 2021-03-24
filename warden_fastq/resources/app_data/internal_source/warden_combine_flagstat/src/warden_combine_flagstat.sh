#!/bin/bash
# shellcheck disable=SC2154

main() {
    echo "Start combine flagstat"
    dx-download-all-inputs --parallel
    mkdir FLAGSTAT
    mv in/flagstat_files/*/* FLAGSTAT
    process_flagstat.py "$sample_list_path"
    combined_flagstat=$(dx upload alignment_statistics.txt --brief)
    dx tag "$combined_flagstat" sjcp-result-file
    dx-jobutil-add-output combined_flagstat "$combined_flagstat" --class=file
}
