#!/bin/bash
# shellcheck disable=SC2154
set -e -x -o pipefail

system_stats() {
    while true; do
        echo "###### MEMORY & SPACE ######"
        free -m
        df
        echo "############################"
        sleep 60
    done
}

main() {
    wget -nv "https://github.com/arq5x/bedtools2/releases/download/v2.29.2/bedtools.static.binary"
    mv bedtools.static.binary /usr/bin/bedtools
    chmod +x /usr/bin/bedtools
    wget -nv "https://github.com/samtools/samtools/releases/download/1.11/samtools-1.11.tar.bz2" -O - | tar xj
    cd samtools-1.11
    make -s
    make -s install
    cd ..

    dx mkdir -p /tmp

    dx-download-all-inputs --parallel

    # system_stats&

    sorted_bam="sorted.bam"
    if [ "$sorted" = "false" ]; then
        samtools sort "$input_bam_path" > "$sorted_bam"
    else
        sorted_bam=$input_bam_path
    fi

    samtools_command="samtools flagstat $sorted_bam > unsplit.flagstat.txt"
    eval "$samtools_command"

    counts=$(get_mapped_counts.py unsplit.flagstat.txt)

    default_scale=0
    if [ "$scale" == "$default_scale" ]; then
        scale=$(python -c "print float(10**7)/float($counts)")
    else
        scale=1
    fi

    genome_coverage_command="bedtools genomecov $output_options -scale $scale -split -ibam $sorted_bam -g $genome_sizes_file_path > coverage_file.bed"
    eval "$genome_coverage_command"

    sort_command="LC_COLLATE=C sort -k1,1 -k2,2n coverage_file.bed > $output_prefix.coverage_file.bed"
    eval "$sort_command"

    all_out=$(dx upload --path /tmp/ "$output_prefix.coverage_file.bed" --brief)
    dx-jobutil-add-output all_coverage_file "$all_out" --class=file

    if [ "$strandedness" != "no" ]; then
        samtools view -h -f 80 "$sorted_bam" > s_r1.bam
        samtools view -h -f 160 "$sorted_bam" > s_f2.bam
        samtools view -h -f 96 "$sorted_bam" > as_f1.bam
        samtools view -h -f 144 "$sorted_bam" > as_r2.bam
        if [ "$strandedness" = "yes" ]; then
            samtools merge -f sense.bam as_f1.bam as_r2.bam
            samtools merge -f antisense.bam s_r1.bam s_f2.bam
        else
            samtools merge -f sense.bam s_r1.bam s_f2.bam
            samtools merge -f antisense.bam as_f1.bam as_r2.bam
        fi

        gcb_sense_command="bedtools genomecov $output_options -scale $scale -split -ibam sense.bam -g $genome_sizes_file_path > pos_coverage_file.bed"
        eval "$gcb_sense_command"
        gcb_antisense_command="bedtools genomecov $output_options -scale $scale -split -ibam antisense.bam -g $genome_sizes_file_path > tmp_neg_coverage_file.bed"
        eval "$gcb_antisense_command"

        SORT1COMMAND="LC_COLLATE=C sort -k1,1 -k2,2n pos_coverage_file.bed > $output_prefix.pos_coverage_file.bed"
        eval "$SORT1COMMAND"
        SORT2COMMAND="LC_COLLATE=C sort -k1,1 -k2,2n tmp_neg_coverage_file.bed > tmp_sorted_neg_coverage_file.bed"
        eval "$SORT2COMMAND"

        awk 'BEGIN{FS="\t";OFS="\t";} {$4=-$4 ; print $0}' tmp_sorted_neg_coverage_file.bed > "$output_prefix".neg_coverage_file.bed

        # FLAGSTATCOMMAND1="samtools flagstat sense.bam > sense.flagstat.txt"
        # eval "$FLAGSTATCOMMAND1"
        # FLAGSTATCOMMAND2="samtools flagstat antisense.bam > antisense.flagstat.txt"
        # eval "$FLAGSTATCOMMAND2"

        # sense_counts=$(get_mapped_counts.py sense.flagstat.txt)
        # anti_sense_counts=$(get_mapped_counts.py antisense.flagstat.txt)
        # final_counts=$(($sense_counts + $anti_sense_counts))
        # if [ "$final_counts" -ne "$counts" ]
        # then
        #     dx-jobutil-report-error "Counts before and after splitting do not match Before=$counts after=$final_counts"
        # fi

        pos_out=$(dx upload --path /tmp/ "$output_prefix.pos_coverage_file.bed" --brief)
        dx-jobutil-add-output pos_coverage_file "$pos_out" --class=file
        neg_out=$(dx upload --path /tmp/ "$output_prefix.neg_coverage_file.bed" --brief)
        dx-jobutil-add-output neg_coverage_file "$neg_out" --class=file
    fi
}
