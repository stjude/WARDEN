#!/bin/bash
# warden_simple_differential_expression 1.0.4
# shellcheck disable=SC2154
set -e -x -o

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
    ./configure --with-readline=no --with-x=no --with-cairo=yes --with-libpng=yes --with-recommended-packages=no --enable-R-profiling=no > /dev/null
    make -s -j "$(nproc)" > /dev/null
    make -s -j "$(nproc)" install > /dev/null
    cd ..

    cd /R_simple_DE
    R -e 'install.packages("packrat",repos="http://cran.us.r-project.org", quiet=TRUE)'
    R -e 'packrat::restore()'
    cd

    dx mkdir -p DIFFERENTIAL_EXPRESSION
    dx mkdir -p AUXILIARY
    mkdir -p out/out_files
    mkdir -p out/gz_out_files

    dx download "$input_count_file" -o /R_simple_DE/count_file.txt
    dx download "$sample_list_file" -o /R_simple_DE/sample_list.txt

    if [ -n "$contrasts_file" ]; then
        dx download "$contrasts_file" -o /R_simple_DE/contrasts_file.txt
        contrasts=$(cat /R_simple_DE/contrasts_file.txt)
        contrasts_final=$(parse_contrasts.py "$contrasts")
    elif [ -n "$contrasts" ]; then
        contrasts_final=$(parse_contrasts.py "$contrasts")
    else
        contrasts_final=""
    fi
    echo -e "contrasts\t$contrasts_final" >> /R_simple_DE/R_parameters_simple.txt

    cd /R_simple_DE
    R --slave -f simple_DE.R
    cd

    mkdir -p out/out_files

    if [ -z "$difex_viewer" ] || [ "$difex_viewer" == "None" ]; then
        echo "No viewer set. Not creating shortcuts."
    else
        viewer_id=$(jq --raw-output ".id" <<< "$(dx describe --json "$difex_viewer")")
        echo "Viewer id:$viewer_id"
    fi

    upload_id=$(dx upload --path mds_plot.norm_cpm.png /R_simple_DE/mds_plot.norm_CPM.png --brief)
    dx tag "$upload_id" sjcp-result-file
    dx-jobutil-add-output --array out_files "$upload_id"
    upload_id=$(dx upload --path DIFFERENTIAL_EXPRESSION/ /R_simple_DE/CPM_normalized_values.txt --brief)
    dx tag "$upload_id" sjcp-result-file
    dx-jobutil-add-output --array out_files "$upload_id"

    while read -r contrast; do
        (
            upload_id=$(dx upload --path DIFFERENTIAL_EXPRESSION/MA_plot.simple_DE."$contrast".png /R_simple_DE/MA_plot.simple_DE."$contrast".png --brief)
            dx tag "$upload_id" sjcp-result-file
            dx-jobutil-add-output --array out_files "$upload_id"

            if [ -z "$difex_viewer" ] || [ "$difex_viewer" == "None" ]; then
                upload_id=$(dx upload --path DIFFERENTIAL_EXPRESSION/results.simple_DE."$contrast".txt /R_simple_DE/simple_DE."$contrast".txt --brief)
                dx tag "$upload_id" sjcp-result-file
                dx-jobutil-add-output --array out_files "$upload_id"
            else
                (
                    upload_id=$(dx upload --path DIFFERENTIAL_EXPRESSION/results.simple_DE."$contrast".txt /R_simple_DE/simple_DE."$contrast".txt --brief)
                    dx tag "$upload_id" sjcp-result-file
                    dx-jobutil-add-output --array out_files "$upload_id"
                    details="{\"fileViewer\":{\"id\":\"$viewer_id\"},\"preselectedIDs\":[\"$upload_id\"]}"

                    bookmark_id=$(dx new record simple_DE."$contrast".viewer \
                        --details "$details" \
                        --type ViewerShortcut \
                        --type SJCloudVisualization \
                        --brief \
                        --close)
                    dx-jobutil-add-output viewer_bookmark "$bookmark_id" --class=record
                )
            fi
        )
    done < /R_simple_DE/contrast_files.txt
}
