#!/bin/bash
# shellcheck disable=SC2154
set -e -x -o
main() {
    wget -nv -O jq-1.6 https://github.com/stedolan/jq/releases/download/jq-1.6/jq-linux64
    chmod +x ./jq-1.6

    dx-download-all-inputs --parallel

    printf '%s\n' "${sample_files_path[@]}" > sample_files_list.txt
    printf '%s\n' "${count_files_path[@]}" > count_files_list.txt

    echo "Starting process_count_files.py"
    process_count_files.py sample_files_list.txt count_files_list.txt "combined_counts.$name_value.txt" "combined_sample_file.$name_value.txt"
    echo "finished process_count_files.py"
    err_length=$(stat -c%s errors.txt)
    if [ "$err_length" -gt 0 ]; then
        errors=$(cat errors.txt)
        echo "Error message: $error"
        dx-jobutil-report-error "$errors"
    fi

    WORKFLOW_ID=$(./jq-1.6 --raw-output ".analysis" <<< "$(dx describe --json "$DX_JOB_ID")")
    PARENT_JOB=$(./jq-1.6 --raw-output ".workflow.createdBy.job" <<< "$(dx describe --json "$WORKFLOW_ID")")
    PARENT_FOLDER=$(./jq-1.6 --raw-output ".folder" <<< "$(dx describe --json "$PARENT_JOB")")
    APPFOLDER=$(./jq-1.6 --raw-output ".folder" <<< "$(dx describe --json "$DX_JOB_ID")")

    OUTPATH="$DX_PROJECT_CONTEXT_ID:$PARENT_FOLDER/$APPFOLDER/"

    dx mkdir -p "$OUTPATH"
    count_file=$(dx upload --path "$OUTPATH" "combined_counts.$name_value.txt" --brief)
    sample_file=$(dx upload "combined_sample_file.$name_value.txt" --brief)
    dx tag "$DX_PROJECT_CONTEXT_ID:$count_file" sjcp-result-file
    dx tag "$sample_file" sjcp-result-file
    dx-jobutil-add-output count_file "$count_file" --class=file
    dx-jobutil-add-output sample_file "$sample_file" --class=file
}
