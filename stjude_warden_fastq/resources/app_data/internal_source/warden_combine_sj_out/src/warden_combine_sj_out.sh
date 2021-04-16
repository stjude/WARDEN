#!/bin/bash
# shellcheck disable=SC2154
# warden_combine_sj_out 1.0.2

main() {
    wget -nv "https://github.com/arq5x/bedtools2/releases/download/v2.29.2/bedtools.static.binary"
    mv bedtools.static.binary /usr/bin/bedtools
    chmod +x /usr/bin/bedtools

    for i in "${!sj_out_files[@]}"; do
        dx download "${sj_out_files[$i]}" -o sj_out_files-$i
        cat sj_out_files-$i > all_files.bed
    done
    sort -k1,1 -k2,2n all_files.bed > all_files.sorted.bed
    bedtools merge -i all_files.sorted.bed > sj_out.merged.bed

    combined_sj_out=$(dx upload sj_out.merged.bed --brief)
    dx-jobutil-add-output combined_sj_out "$combined_sj_out" --class=file
}
