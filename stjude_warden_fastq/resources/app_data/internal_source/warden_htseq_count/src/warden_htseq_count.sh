#!/bin/bash
# shellcheck disable=SC2154
set -e -x -o pipefail

main() {
    python3 -m pip install --no-cache-dir HTSeq
    wget -nv -O jq-1.6 https://github.com/stedolan/jq/releases/download/jq-1.6/jq-linux64
    chmod +x ./jq-1.6

    echo "Value of feature: $feature_type"
    echo "Value of gff_feature: $id_attribute"
    echo "Value of strand: $strand"
    min_qual="10"
    dx-download-all-inputs --parallel
    out_name=$prefix.htseq_counts.txt

    echo -e "$id_attribute\t$prefix" > "$out_name"
    htseq_command="htseq-count -f bam -r ${order} -i ${id_attribute} -s ${strand} \
        -t ${feature_type} -m ${mode} --nonunique=${nonunique} \
        --secondary-alignments=${secondary_alignments} \
        --supplementary-alignments=${supplementary_alignments} \
        -a $min_qual \
        $input_bam_path $annotation_file_path \
        >> $out_name"
    eval "$htseq_command"

    if [ -n "$gene_length_file" ] && [ "$id_attribute" = "gene_name" ]; then
        fpkm_file=$prefix.fpkm.txt
        fpkmlog2File=$prefix.fpkm.log2.txt
        python3 /usr/bin/calc_fpkm.py "$out_name" "$gene_length_file_path" "$fpkm_file" "$fpkmlog2File"

        dx mkdir -p FPKM
        fpkm=$(dx upload --path FPKM/ "$fpkm_file" --brief)
        dx tag "$fpkm" sjcp-result-file
        dx-jobutil-add-output fpkm "$fpkm" --class=file

        fpkmlog2=$(dx upload --path FPKM/ "$fpkmlog2File" --brief)
        dx tag "$fpkmlog2" sjcp-result-file
        dx-jobutil-add-output fpkm_log2 "$fpkmlog2" --class=file
    fi

    WORKFLOW_ID=$(./jq-1.6 --raw-output ".analysis" <<< "$(dx describe --json "$DX_JOB_ID")")
    PARENT_JOB=$(./jq-1.6 --raw-output ".workflow.createdBy.job" <<< "$(dx describe --json "$WORKFLOW_ID")")
    PARENT_FOLDER=$(./jq-1.6 --raw-output ".folder" <<< "$(dx describe --json "$PARENT_JOB")")
    APPFOLDER=$(./jq-1.6 --raw-output ".folder" <<< "$(dx describe --json "$DX_JOB_ID")")

    OUTPATH="$DX_PROJECT_CONTEXT_ID:$PARENT_FOLDER/$APPFOLDER/"

    dx mkdir -p "$OUTPATH"
    htseq_counts=$(dx upload --path "$OUTPATH" "$out_name" --brief)
    dx tag "$DX_PROJECT_CONTEXT_ID:$htseq_counts" sjcp-result-file
    dx-jobutil-add-output htseq_counts "$htseq_counts" --class=file
}
