#!/bin/bash
# warden_voom_limma 1.0.0
# shellcheck disable=SC2154
set -e -x -o

mkdir -p out/out_files
mkdir -p out/viewer_files
mkdir out/sessionInfo

# upload_file() {
#     upfile=$1
#     output_name=$2
#     is_array=$3
#     read -r -a path <<< "$4"
#     upload_id=$(dx upload "${path[@]}" "$upfile" --brief)
#     dx-jobutil-add-output ${is_array:+"--array"} "$output_name" "$upload_id"
# }

# upload_file_tag() {
#     upfile=$1
#     output_name=$2
#     is_array=$3
#     read -r -a path <<< "$4"
#     upload_id=$(dx upload "${path[@]}" "$upfile" --brief)
#     dx tag "$upload_id" sjcp-result-file
#     dx-jobutil-add-output ${is_array:+"--array"} "$output_name" "$upload_id"
# }

main() {
    wget -nv https://cran.r-project.org/src/base/R-3/R-3.5.3.tar.gz
    tar -xzf ./R-3.5.3.tar.gz
    cd R-3.5.3/
    ./configure \
        --with-readline=no \
        --with-x=no \
        --with-cairo=yes \
        --with-libpng=yes \
        --with-recommended-packages=no \
        --enable-R-profiling=no \
        > /dev/null
    make -s -j "$(nproc)" > /dev/null
    make -s -j "$(nproc)" install > /dev/null

    cd /R_voom_limma
    R -e 'install.packages("packrat",repos="http://cran.us.r-project.org", quiet=TRUE)'
    R -e 'packrat::restore()'
    cd

    dx mkdir -p AUXILIARY
    dx mkdir -p DIFFERENTIAL_EXPRESSION

    dx download "$input_count_file" -o /R_voom_limma/count_file.txt
    dx download "$sample_list_file" -o /R_voom_limma/limma_sample_list.txt

    echo "contrastfile:$contrasts_file"
    if [ -n "$contrasts_file" ]; then
        dx download "$contrasts_file" -o /R_voom_limma/contrasts_file.txt
        contrasts=$(cat /R_voom_limma/contrasts_file.txt)
        contrasts_final=$(parse_contrasts.py "$contrasts")
    elif [ -n "$contrasts" ]; then
        contrasts_final=$(parse_contrasts.py "$contrasts")
    else
        contrasts_final=""
    fi

    {
        echo -e "calcNormFactors_method\t$calcNormFactors_method"
        echo -e "filter_count_type\t$filter_count_type"
        echo -e "filter_count\t$filter_count"
        echo -e "p_value_adjust\t$p_value_adjust"
        echo -e "contrasts\t$contrasts_final"
    } > /R_voom_limma/R_parameters.txt
    cat /R_voom_limma/R_parameters.txt

    cd /R_voom_limma
    R --slave -f voom_limma.R
    cd

    upload_id=$(dx upload --path mds_plot.limma.png /R_voom_limma/mds_plot.png --brief)
    dx tag "$upload_id" sjcp-result-file
    dx-jobutil-add-output --array out_files "$upload_id"
    upload_id=$(dx upload --path mean_variance.png /R_voom_limma/mean_variance.png --brief)
    dx tag "$upload_id" sjcp-result-file
    dx-jobutil-add-output --array out_files "$upload_id"

    if [ -z "$difex_viewer" ] || [ "$difex_viewer" == "None" ]; then
        echo "No viewer set. Not creating shortcuts."
    else
        viewer_id=$(jq --raw-output ".id" <<< "$(dx describe --json "$difex_viewer")")
        echo "Viewer id:$viewer_id"
    fi

    while read -r contrast; do
        (
            upload_id=$(dx upload --path AUXILIARY/ /R_voom_limma/GSEA.tStat."$contrast".txt --brief)
            dx tag "$upload_id" sjcp-result-file
            dx-jobutil-add-output --array out_files "$upload_id"
            upload_id=$(dx upload --path AUXILIARY/ /R_voom_limma/GSEA.input."$contrast".txt --brief)
            dx tag "$upload_id" sjcp-result-file
            dx-jobutil-add-output --array out_files "$upload_id"
            upload_id=$(dx upload --path DIFFERENTIAL_EXPRESSION/MA_plot.limma."$contrast".png /R_voom_limma/MA_plot."$contrast".png --brief)
            dx tag "$upload_id" sjcp-result-file
            dx-jobutil-add-output --array out_files "$upload_id"
            upload_id=$(dx upload --path DIFFERENTIAL_EXPRESSION/volcano_plot.limma."$contrast".png /R_voom_limma/volcano_plot."$contrast".png --brief)
            dx tag "$upload_id" sjcp-result-file
            dx-jobutil-add-output --array out_files "$upload_id"
            if [ -z "$difex_viewer" ] || [ "$difex_viewer" == "None" ]; then
                upload_id=$(dx upload --path DIFFERENTIAL_EXPRESSION/results.limma."$contrast".txt /R_voom_limma/results."$contrast".txt --brief)
                dx tag "$upload_id" sjcp-result-file
                dx-jobutil-add-output --array out_files "$upload_id"
            else
                (
                    upload_id=$(dx upload --path DIFFERENTIAL_EXPRESSION/results.limma."$contrast".txt /R_voom_limma/results."$contrast".txt --brief)
                    dx tag "$upload_id" sjcp-result-file
                    dx-jobutil-add-output --array out_files "$upload_id"
                    details="{\"fileViewer\":{\"id\":\"$viewer_id\"},\"preselectedIDs\":[\"$upload_id\"]}"

                    bookmark_id=$(dx new record limma."$contrast".viewer \
                        --details "$details" \
                        --type ViewerShortcut \
                        --type SJCloudVisualization \
                        --brief \
                        --close)
                    dx-jobutil-add-output viewer_bookmark "$bookmark_id" --class=record
                )
            fi
        )
    done < /R_voom_limma/contrast_files.txt
}
