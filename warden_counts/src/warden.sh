#!/bin/bash
# shellcheck disable=SC2154

wget -nv -O jq-1.6 https://github.com/stedolan/jq/releases/download/jq-1.6/jq-linux64
chmod +x ./jq-1.6

echo "########################START APP################################"
echo "JOBID: $DX_JOB_ID"
echo "Job workspace: $DX_WORKSPACE_ID"

APPFOLDER=$(./jq-1.6 --raw-output ".folder" <<< "$(dx describe --json "$DX_JOB_ID")")
echo "APPFOLDER: $APPFOLDER"
FOLDER_PROVIDED=""
if [ "$APPFOLDER" = "" ] || [ "$APPFOLDER" = "/" ]; then
    APPFOLDER="/WARDEN_$(echo "$DX_JOB_ID" | cut -d '-' -f 2)"
    echo "No output folder given."
    echo "Creating and using $APPFOLDER"
    FOLDER_PROVIDED="false"
else
    FOLDER_PROVIDED="true"
fi

echo "Building applet dependencies..."
for app in /app_data/internal_source/*; do
    dx build -f "$app"
done
echo "Applets built"
echo ""
if [ -n "$Genome" ]; then
    Genome=$(echo "$Genome" | awk '{print $1}')
    genome_json=/app_data/genome_data.json
    limma_DE_viewer=$(./jq-1.6 --raw-output ".$Genome.viewers.LIMMA_DifEx_Viewer" $genome_json)
    echo "DIFEXVIEWER: $limma_DE_viewer"
    echo ""
fi

############################INPUT FILES#############################################################
sample_list_extension="${sample_list_path##*.}"
sample_list_extension="${sample_list_extension,,}"

main() {
    dx download "$sample_list" -o "sample_list.$sample_list_extension"
    if [ "$sample_list_extension" == "txt" ]; then
        dos2unix -q sample_list.txt
    elif [ "$sample_list_extension" == "xlsx" ]; then
        python /usr/bin/parse_excel_sample_list.py sample_list.xlsx > sample_list.txt
    else
        dx-jobutil-report-error "Improper Sample List Extension. This should be a .txt or .xlsx file" appError
    fi

    printf '%s\n' "${COUNT_FILES_path[@]}" > count_list.txt
    for count in "${COUNT_FILES[@]}"; do
        dx describe "$count" --delim ',' | grep ID | cut -f 2 -d ',' >> count_link_list.txt
    done

    echo "###PROCESSING FILES###"
    process_files.py "true" sample_list.txt count_list.txt count_link_list.txt

    PROCESSFILE_ERR=$(cat process_files_errors.txt)
    IS_PROCESSFILE_ERR=${#PROCESSFILE_ERR} #get size
    if [ "$IS_PROCESSFILE_ERR" -gt 0 ]; then
        echo "Error: $PROCESSFILE_ERR"
        dx-jobutil-report-error "$PROCESSFILE_ERR" appError
    fi

    final_sample_list_id=$(dx upload --brief cleaned_sample_list.txt)
    comparisons_limma_id=$(dx upload --brief comparisons_limma.txt)
    comparisons_all_id=$(dx upload --brief comparisons_simple_DE.txt)

    echo "count_file_list.txt:"
    cat count_file_list.txt
    echo "cleaned_sample_list.txt:"
    cat cleaned_sample_list.txt
    echo ""

    contrasts_limma=$(cat comparisons_limma.txt)
    echo "contrastslimma:$contrasts_limma"
    contrasts_simple=$(cat comparisons_simple_DE.txt)
    echo "contrastssimple:$contrasts_simple"

    if [ "$(wc -m < comparisons_limma.txt)" -gt 3 ]; then
        limma_runnable="true"
    else
        limma_runnable="false"
    fi

    num_samples=$(($(wc -l < cleaned_sample_list.txt) - 1))
    echo "Number of samples: $num_samples"
    echo ""
    if [ "$num_samples" -gt 64 ]; then
        echo "Error: Number of samples greater than 64.  The app limits samples to 64"
        dx-jobutil-report-error "Number of samples greater than 64.  The app limits samples to 64" appError
    fi
    ###############
    {
        echo -e "project\t$DX_PROJECT_CONTEXT_ID"
        echo -e "Output\t$APPFOLDER"
        echo -e "folder_provided\t$FOLDER_PROVIDED"
        echo -e "app_version\t$appversion"
        echo -e "calc_norm_factors_method\t$calc_norm_factors_method"

        echo -e "filter_count_type\t$filter_count_type"
        echo -e "filter_count\t$filter_count"
        echo -e "p_value_adjust\t$p_value_adjust"
        echo -e "contrasts\t$contrasts_limma"
        echo -e "limma_runnable\t$limma_runnable"
        echo -e "limma_DE_viewer\t$limma_DE_viewer"
    } > workflow_parameter_file.txt

    dx mkdir -p "$DX_PROJECT_CONTEXT_ID:$APPFOLDER/AUXILIARY"
    echo ""

    echo "RUNNING create workflow"
    wf_id=$(create_workflow.py workflow_parameter_file.txt count_file_list.txt "$final_sample_list_id" "$comparisons_limma_id" "$comparisons_all_id")

    dx describe --json "$DX_JOB_ID" | ./jq-1.6 '.input' > WARDEN_parameters.json
    params_id=$(dx upload --path "$DX_PROJECT_CONTEXT_ID:$APPFOLDER/AUXILIARY/" WARDEN_parameters.json --brief)

    {
        echo -e "{\n"
        echo "\"parameters\":{\"\$dnanexus_link\":\"$params_id\"},"
    } > job_output.json

    echo "LAUNCHING WORKFLOW"
    analysis_id=$(dx run -y "$wf_id" --brief --extra-args '{"executionPolicy": {"onNonRestartableFailure": "failStage"}}')
    echo "Waiting for $analysis_id to complete..."
    dx wait "$analysis_id"

    dx describe --json "$analysis_id" | ./jq-1.6 '.output' > analysis.json
    {
        echo "$(./jq-1.6 'with_entries(if (.key | test("combined_counts")) then ( {key: .key, value: .value } ) else empty end) | flatten | { "combined_counts": .[0] }' analysis.json | head -n -1 | tail -n +2)",
        echo "$(./jq-1.6 'with_entries(if (.key | test("limma_outfiles")) then ( {key: .key, value: .value } ) else empty end) | flatten | { "limma_outfiles": . }' analysis.json | head -n -1 | tail -n +2)",
        echo "$(./jq-1.6 'with_entries(if (.key | test("limma_viewer")) then ( {key: .key, value: .value } ) else empty end) | flatten | { "limma_viewer": .[0] }' analysis.json | head -n -1 | tail -n +2)",
        echo "$(./jq-1.6 'with_entries(if (.key | test("simple_DE_outfiles")) then ( {key: .key, value: .value } ) else empty end) | flatten | { "simple_DE_outfiles": . }' analysis.json | head -n -1 | tail -n +2)",
        echo "$(./jq-1.6 'with_entries(if (.key | test("simple_DE_viewer")) then ( {key: .key, value: .value } ) else empty end) | flatten | { "simple_DE_viewer": .[0] }' analysis.json | head -n -1 | tail -n +2)"'}'
    } >> job_output.json
}
