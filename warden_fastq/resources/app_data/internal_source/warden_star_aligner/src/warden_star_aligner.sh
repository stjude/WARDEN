#!/bin/bash
# shellcheck disable=SC2154
set -e -x -o pipefail

system_stats() {
    while true; do
        echo "###### MEMORY & SPACE ######"
        free -m
        df
        sleep 60
    done
}

main() {
    wget -nv https://github.com/broadinstitute/picard/releases/download/2.23.8/picard.jar
    mv picard.jar /usr/bin/
    chmod +x /usr/bin/picard.jar
    wget -nv https://github.com/alexdobin/STAR/raw/2.5.3a/bin/Linux_x86_64/STAR
    mv STAR /usr/bin
    chmod +x /usr/bin/STAR
    wget -nv https://github.com/samtools/samtools/releases/download/1.11/samtools-1.11.tar.bz2 -O - | tar -xj
    cd samtools-1.11
    make -s
    make -s install
    cd ..
    wget -nv https://github.com/lh3/seqtk/archive/v1.3.tar.gz -O - | tar -xz
    cd seqtk-1.3
    make
    make install
    cd ..
    wget -nv -O jq-1.6 https://github.com/stedolan/jq/releases/download/jq-1.6/jq-linux64
    chmod +x ./jq-1.6

    ###SETUP#####
    STARINDEXPATH="/home/dnanexus/in/star_index_archive/STARINDEX"
    mkdir -p $STARINDEXPATH
    #download and unzip index at same time
    dx cat "$star_index_archive" | pigz -d - | tar -xf - --owner root --group root --no-same-owner -C $STARINDEXPATH --strip-components=1 &

    file1_ext="${read_file1_name##*.}"
    read_file1_name="$output_prefix.r1.fastq.$file1_ext"
    dx download "$read_file1" -o "$read_file1_name" &
    if [ -n "$read_file2_path" ]; then
        file2_ext="${read_file2_name##*.}"
        read_file2_name="$output_prefix.r2.fastq.$file2_ext"
        dx download "$read_file2" -o "$read_file2_name" &
    fi
    sjdbFileChrStartEnd_line=""
    if [ -n "$sjdbFileChrStartEnd" ]; then
        dx download "$sjdbFileChrStartEnd" -o sjdbFileChrStartEnd.bed &
        sjdbFileChrStartEnd_line="--sjdbFileChrStartEnd sjdbFileChrStartEnd.bed"
    fi
    annotation_line=""
    overhang_line=""
    #if  gtf file is provided, (otherwise assume index has been indexed with annotations)
    if [ -n "$transcriptome_gtf" ]; then
        dx download "$transcriptome_gtf" -o annotation.gtf &
        annotation_line="--sjdbGTFfile annotation.gtf"
        overhang_line="--sjdbOverhang $sjdbOverhang"
    fi
    wait

    quantMode=""
    if [ "$generate_transcriptome_BAM" = "true" ]; then
        quantMode="--quantMode TranscriptomeSAM"
    fi

    read_file_names="$read_file1_name"
    readFilesCommand=""
    if [ "$read_file2_path" ]; then
        if [ "$file1_ext" != "$file2_ext" ]; then
            echo 2> "File extensions do not match"
            exit 1
        fi
        read_file_names="$read_file1_name $read_file2_name"
    fi
    if [ "$file1_ext" = "gz" ]; then
        readFilesCommand="--readFilesCommand zcat"
    elif [ "$file1_ext" = "bz2" ]; then
        readFilesCommand="--readFilesCommand bzip2 -c"
    fi
    output_prefix="$output_prefix."

    outSAMtype=""
    if [ "$first_pass" = "true" ]; then
        outSAMtype="--outSAMtype None --outSAMmode None"
    else
        outSAMtype="--outSAMtype BAM Unsorted"
    fi

    if [ "$subsample_target" -ge 1 ]; then
        subsample_read1_name="$(basename "$read_file1_name" "$file1_ext")"subsampled.fq
        seqtk sample -s100 "$read_file1_name" "$subsample_target" > "$subsample_read1_name"
        read_file_names="$subsample_read1_name"
        if [ "$read_file2_path" ]; then
            subsample_read2_name="$(basename "$read_file2_name" "$file2_ext")"subsampled.fq
            seqtk sample -s100 "$read_file2_name" "$subsample_target" > "$subsample_read2_name"
            read_file_names="$subsample_read1_name $subsample_read2_name"
        fi
        readFilesCommand=""
    fi

    COMMAND="STAR \
        --runMode alignReads \
        --genomeDir $STARINDEXPATH \
        --readFilesIn $read_file_names \
        --outFileNamePrefix $output_prefix \
        --outSAMunmapped $outSAMunmapped \
        --runThreadN $(nproc) \
        $annotation_line \
        $quantMode \
        $overhang_line \
        --outSAMattributes $outSAMattributes \
        --outFilterMultimapNmax $outFilterMultimapNmax \
        --outFilterMismatchNmax $outFilterMismatchNmax \
        --alignIntronMax $alignIntronMax \
        $sjdbFileChrStartEnd_line \
        --chimSegmentMin $chimSegmentMin \
        --chimJunctionOverhangMin $chimSegmentMin \
        --outSAMstrandField $outSAMstrandField \
        --outBAMcompression $outBAMcompression \
        $outSAMtype \
        $readFilesCommand"
    eval "$COMMAND" || cat ./*Log*

    if [ "$first_pass" = "false" ]; then
        ####REMOVE SECONDARY##########
        out_unsorted_bam="${output_prefix}Aligned.out.bam"
        out_sorted_coord_bam_path="${output_prefix}Aligned.sortedByCoord.out.bam"

        samtools sort --threads "$(nproc)" "$out_unsorted_bam" -o "$out_sorted_coord_bam_path"

        remove_secondary="${output_prefix}Aligned.sorted_by_coord.rm2.out.bam"

        remove_secondary_command="samtools view -b -h -o $remove_secondary -F 256 $out_sorted_coord_bam_path"
        eval "$remove_secondary_command"

        out_final_bam="${output_prefix}Aligned.bam"
        if [[ "$mark_duplicates" = "true" ]]; then
            echo "running mark_duplicates"

            metrics_file="${output_prefix}mark_duplicates_metrics.txt"

            MARKDUPLICATESCOMMAND="java -jar /usr/bin/picard.jar MarkDuplicates \
                I=$remove_secondary \
                o=$out_final_bam \
                m=$metrics_file \
                VERBOSITY=WARNING"
            eval "$MARKDUPLICATESCOMMAND"
        else
            mv "$remove_secondary" "$out_final_bam"
        fi

        #### FLAGSTAT ##########
        flagstat_file="${output_prefix}Aligned.flagstat.txt"

        FLAGSTATCOMMAND="samtools flagstat $out_final_bam > $flagstat_file"
        eval "$FLAGSTATCOMMAND"

        WORKFLOW_ID=$(./jq-1.6 --raw-output ".analysis" <<< "$(dx describe --json "$DX_JOB_ID")")
        PARENT_JOB=$(./jq-1.6 --raw-output ".workflow.createdBy.job" <<< "$(dx describe --json "$WORKFLOW_ID")")
        PARENT_FOLDER=$(./jq-1.6 --raw-output ".folder" <<< "$(dx describe --json "$PARENT_JOB")")
        APPFOLDER=$(./jq-1.6 --raw-output ".folder" <<< "$(dx describe --json "$DX_JOB_ID")")

        OUTPATH="$DX_PROJECT_CONTEXT_ID:$PARENT_FOLDER/$APPFOLDER/"

        dx mkdir -p "$OUTPATH"
        dx mkdir -p "$OUTPATH"/LOGS
        dx mkdir -p "$OUTPATH"/FLAGSTATS
        out_coord_bam_file=$(dx upload --path "$OUTPATH" "$out_final_bam" --brief)
        dx tag "$DX_PROJECT_CONTEXT_ID:$out_coord_bam_file" sjcp-result-file
        dx-jobutil-add-output sorted_by_coord_bam "$out_coord_bam_file" --class=file

        log_final_out_file="${output_prefix}Log.final.out"
        log_final_out=$(dx upload --path "$OUTPATH"/LOGS/ "$log_final_out_file" --brief)
        dx-jobutil-add-output log_final_out "$log_final_out" --class=file

        flagstat_out=$(dx upload --path "$OUTPATH"/FLAGSTATS/ "$flagstat_file" --brief)
        dx-jobutil-add-output flagstat_out "$flagstat_out" --class=file

        if [ "$generate_transcriptome_BAM" = "true" ]; then
            dx_transcriptome_dir="TRANSCRIPTOME/"
            dx mkdir -p $dx_transcriptome_dir

            transcriptome_file="${output_prefix}Aligned.toTranscriptome.out.bam"
            out_transcriptome_file="${output_prefix}Aligned.to_transcriptome.bam"
            mv "$transcriptome_file" "$out_transcriptome_file"
            transcriptome_bam_out=$(dx upload --path $dx_transcriptome_dir "$out_transcriptome_file" --brief)
            dx tag "$transcriptome_bam_out" sjcp-result-file
            dx-jobutil-add-output to_transcriptome_bam "$transcriptome_bam_out" --class=file
        fi

        if [ "$chimSegmentMin" -gt 0 ]; then
            dx_chimeric_dir="CHIMERIC/"
            dx mkdir -p $dx_chimeric_dir

            chimeric_sam_file="${output_prefix}Chimeric.out.sam"
            chimeric_bam_file="${output_prefix}Chimeric.out.bam"
            samtools view -b "$chimeric_sam_file" > "$chimeric_bam_file"
            chimeric_bam=$(dx upload --path $dx_chimeric_dir "$chimeric_bam_file" --brief)
            dx tag "$chimeric_bam" sjcp-result-file
            dx-jobutil-add-output chimeric_bam "$chimeric_bam" --class=file

            chimeric_junction_file="${output_prefix}Chimeric.out.junction"
            chimeric_junction=$(dx upload --path $dx_chimeric_dir "$chimeric_junction_file" --brief)
            dx tag "$chimeric_junction" sjcp-result-file
            dx-jobutil-add-output chimeric_junction "$chimeric_junction" --class=file
        fi
    fi

    dx mkdir -p TABS/
    sj_tab_file="${output_prefix}SJ.out.tab"
    sj_tab_out=$(dx upload --path TABS/ "$sj_tab_file" --brief)
    dx tag "$sj_tab_out" sjcp-result-file
    dx-jobutil-add-output sj_tab_out "$sj_tab_out" --class=file
}
