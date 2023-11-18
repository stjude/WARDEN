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

#################Setting up genome annotations to use################
Genome=$(echo "$Genome" | awk '{print $1}')
annotation_description=""
index_description=""
genome_json=/app_data/genome_data.json
index_file=$(./jq-1.6 --raw-output ".$Genome.StarIndex.IndexFile" $genome_json)
gtf_file=$(./jq-1.6 --raw-output ".$Genome.Annotations.GTFFile" $genome_json)
genome_sizes_file=$(./jq-1.6 --raw-output ".$Genome.GenomeSizes.SizeFile" $genome_json)
indexed_with_gtf=$(./jq-1.6 --raw-output ".$Genome.StarIndex.indexedWithGTF" $genome_json)
index_description=$(./jq-1.6 --raw-output ".$Genome.StarIndex.indexDescription" $genome_json)
annotation_description=$(./jq-1.6 --raw-output ".$Genome.Annotations.annotationDescription" $genome_json)
gene_length_file=$(./jq-1.6 --raw-output ".$Genome.Annotations.geneLengthFile" $genome_json)

limma_DE_viewer=$(./jq-1.6 --raw-output ".$Genome.viewers.LIMMA_DifEx_Viewer" $genome_json)
BW_VIEWER=$(./jq-1.6 --raw-output ".$Genome.viewers.BW_VIEWER" $genome_json)

raw_sequencing_strandedness=$(echo "$sequencing_strandedness" | awk '{print $1}')
if [ "$raw_sequencing_strandedness" = "Unstranded" ]; then
    sequencing_strandedness="no"
elif [ "$raw_sequencing_strandedness" = "First" ]; then
    sequencing_strandedness="reverse"
elif [ "$raw_sequencing_strandedness" = "Second" ]; then
    sequencing_strandedness="yes"
fi

echo "### Genome Parameters ###"
echo "Sequencing Strandedness: $sequencing_strandedness"
echo "GENOME: $Genome"
echo "index_file: $index_file"
echo "gtf_file: $gtf_file"
echo "genome_sizes_file: $genome_sizes_file"
echo "indexed_with_gtf: $indexed_with_gtf"
echo "gene_length_file: $gene_length_file"
echo ""
echo "DIFEXVIEWER: $limma_DE_viewer"
echo "BW_VIEWER: $BW_VIEWER"
echo ""

############################INPUT FILES#############################################################
sample_list_extension="${sample_list_path##*.}"
sample_list_extension="${sample_list_extension,,}"

main() {
    dx download "$sample_list" -o "sample_list.$sample_list_extension"
    if [ "$sample_list_extension" == "txt" ]; then
        dos2unix -q sample_list.txt
    elif [ "$sample_list_extension" == "xlsx" ]; then
        python3 /usr/bin/parse_excel_sample_list.py sample_list.xlsx > sample_list.txt
    else
        dx-jobutil-report-error "Improper Sample List Extension. This should be a .txt or .xlsx file" AppError
    fi

    printf '%s\n' "${BAM_FILES_path[@]}" > bam_list.txt
    for bam in "${BAM_FILES[@]}"; do
        dx describe "$bam" --delim ',' | grep ID | cut -f 2 -d ',' >> bam_link_list.txt
    done

    echo "###PROCESSING FILES###"
    process_files.py "true" sample_list.txt bam_list.txt bam_link_list.txt

    PROCESSFILE_ERR=$(cat process_files_errors.txt)
    IS_PROCESSFILE_ERR=${#PROCESSFILE_ERR} #get size
    if [ "$IS_PROCESSFILE_ERR" -gt 0 ]; then
        echo "Error: $PROCESSFILE_ERR"
        dx-jobutil-report-error "$PROCESSFILE_ERR" AppError
    fi

    final_sample_list_id=$(dx upload --brief cleaned_sample_list.txt)
    comparisons_limma_id=$(dx upload --brief comparisons_limma.txt)
    comparisons_all_id=$(dx upload --brief comparisons_simple_DE.txt)

    echo "bam_file_list.txt:"
    cat bam_file_list.txt
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
    ###############
    {
        echo -e "project\t$DX_PROJECT_CONTEXT_ID"
        echo -e "Output\t$APPFOLDER"
        echo -e "folder_provided\t$FOLDER_PROVIDED"
        echo -e "index_file\t$index_file"
        echo -e "gtf_file\t$gtf_file"
        echo -e "genome_sizes_file\t$genome_sizes_file"
        echo -e "indexed_with_gtf\t$indexed_with_gtf"
        echo -e "index_description\t$index_description"
        echo -e "annotation_description\t$annotation_description"
        echo -e "gene_length_file\t$gene_length_file"
        echo -e "app_version\t$appversion"
        echo -e "strandedness\t$sequencing_strandedness"
        echo -e "calcNormFactors_method\t$calcNormFactors_method"

        echo -e "filter_count_type\t$filter_count_type"
        echo -e "filter_count\t$filter_count"
        echo -e "p_value_adjust\t$p_value_adjust"
        echo -e "contrasts\t$contrasts_limma"
        echo -e "limma_runnable\t$limma_runnable"
        echo -e "limma_DE_viewer\t$limma_DE_viewer"
        echo -e "BW_VIEWER\t$BW_VIEWER"
        echo -e "run_coverage\t$run_coverage"
        echo -e "run_limma\t$run_limma"
        echo -e "run_simple_dif_ex\t$run_simple_dif_ex"

        echo -e "sort_order\t$sort_order"
        echo -e "feature_type\t$feature_type"
        echo -e "id_attribute\t$id_attribute"
        echo -e "mode\t$mode"
        echo -e "nonunique\t$nonunique"
        echo -e "secondary_alignments\t$secondary_alignments"
        echo -e "supplementary_alignments\t$supplementary_alignments"
        echo -e "htseq_instance\t$htseq_instance"
        echo -e "combine_counts_instance\t$combine_counts_instance"
    } > workflow_parameter_file.txt

    dx mkdir -p "$DX_PROJECT_CONTEXT_ID:$APPFOLDER/AUXILIARY"
    echo ""

    echo "RUNNING create workflow"
    wf_id=$(create_workflow.py workflow_parameter_file.txt bam_file_list.txt "$final_sample_list_id" "$comparisons_limma_id" "$comparisons_all_id")

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
        echo "$(./jq-1.6 'with_entries(if (.key | test("^(?!combined_).*_htseqcounts"; "")) then ( {key: .key, value: .value } ) else empty end) | flatten | { "htseqcounts": . }' analysis.json | head -n -1 | tail -n +2)",
        echo "$(./jq-1.6 'with_entries(if (.key | test("^(?!combined_).*_fpkm$"; "")) then ( {key: .key, value: .value } ) else empty end) | flatten | { "fpkms": . }' analysis.json | head -n -1 | tail -n +2)",
        echo "$(./jq-1.6 'with_entries(if (.key | test("^(?!combined_).*_fpkm_log2"; "")) then ( {key: .key, value: .value } ) else empty end) | flatten | { "fpkm_log2s": . }' analysis.json | head -n -1 | tail -n +2)",
        echo "$(./jq-1.6 'with_entries(if (.key | test("_all_bigwig")) then ( {key: .key, value: .value } ) else empty end) | flatten | { "all_bigwigs": . }' analysis.json | head -n -1 | tail -n +2)",
        echo "$(./jq-1.6 'with_entries(if (.key | test("_pos_bigwig")) then ( {key: .key, value: .value } ) else empty end) | flatten | { "pos_bigwigs": . }' analysis.json | head -n -1 | tail -n +2)",
        echo "$(./jq-1.6 'with_entries(if (.key | test("_neg_bigwig")) then ( {key: .key, value: .value } ) else empty end) | flatten | { "neg_bigwigs": . }' analysis.json | head -n -1 | tail -n +2)",
        echo "$(./jq-1.6 'with_entries(if (.key | test("combined_counts")) then ( {key: .key, value: .value } ) else empty end) | flatten | { "combined_counts": .[0] }' analysis.json | head -n -1 | tail -n +2)",
        echo "$(./jq-1.6 'with_entries(if (.key | test("^combined_fpkm$"; "")) then ( {key: .key, value: .value } ) else empty end) | flatten | { "combined_fpkm": .[0] }' analysis.json | head -n -1 | tail -n +2)",
        echo "$(./jq-1.6 'with_entries(if (.key | test("^combined_fpkm_log2$"; "")) then ( {key: .key, value: .value } ) else empty end) | flatten | { "combined_fpkm_log2": .[0] }' analysis.json | head -n -1 | tail -n +2)",
        echo "$(./jq-1.6 'with_entries(if (.key | test("limma_outfiles")) then ( {key: .key, value: .value } ) else empty end) | flatten | { "limma_outfiles": . }' analysis.json | head -n -1 | tail -n +2)",
        echo "$(./jq-1.6 'with_entries(if (.key | test("limma_viewer")) then ( {key: .key, value: .value } ) else empty end) | flatten | { "limma_viewer": .[0] }' analysis.json | head -n -1 | tail -n +2)",
        echo "$(./jq-1.6 'with_entries(if (.key | test("simple_DE_outfiles")) then ( {key: .key, value: .value } ) else empty end) | flatten | { "simple_DE_outfiles": . }' analysis.json | head -n -1 | tail -n +2)",
        echo "$(./jq-1.6 'with_entries(if (.key | test("simple_DE_viewer")) then ( {key: .key, value: .value } ) else empty end) | flatten | { "simple_DE_viewer": .[0] }' analysis.json | head -n -1 | tail -n +2)",
        echo "$(./jq-1.6 'with_entries(if (.key | test("bw_viewer")) then ( {key: .key, value: .value } ) else empty end) | flatten | { "bw_viewer": .[0] }' analysis.json | head -n -1 | tail -n +2)"'}'
    } >> job_output.json
}
