#!/bin/bash
# shellcheck disable=SC2154
set -e -x -o

system_stats() {
    while true; do
        echo "###### MEMORY & SPACE ######"
        free -m
        df
        sleep 60
    done
}
main() {
    wget -nv "https://github.com/samtools/samtools/releases/download/1.11/samtools-1.11.tar.bz2" -O - | tar xj
    cd samtools-1.11
    make -s
    make -s install
    cd ..

    dx download "$input_bam" -o input_bam

    base_name=$(basename "$input_bam_name" .bam)
    output_name="${base_name}.name_sorted.bam"

    sort_command="samtools sort -n -T temp -o $output_name input_bam"
    eval "$sort_command"

    output_bam=$(dx upload "$output_name" --brief)
    dx-jobutil-add-output output_bam "$output_bam" --class=file
}
