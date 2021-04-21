#!/bin/bash
# shellcheck disable=SC2154
# warden_bw_to_viewer 1.0.3
set -e -o -x

main() {

    echo "Value of viewer: '$viewer'"
    echo "Value of bigwig_files: '${bigwig_files[*]}'"
    viewer_id=$(jq --raw-output ".id" <<< "$(dx describe --json "$viewer")")
    viewer_project_id=$(jq --raw-output ".project" <<< "$(dx describe --json "$viewer")")
    echo "Viewer id:$viewer_id"
    echo "Viewer project id:$viewer_project_id"
    details="{\"fileViewer\":{\"id\":\"$viewer_id\"},\"preselectedIDs\":["

    for i in "${!bigwig_files[@]}"; do
        file_id=$(jq --raw-output ".id" <<< "$(dx describe --json "${bigwig_files[$i]}")")
        details+="\"$file_id\","
    done
    details="${details::-1}"
    details+="]}"
    echo "$details"

    viewer_bookmark=$(dx new record bigwig_viewer \
        --details "$details" \
        --type ViewerShortcut \
        --type SJCloudVisualization \
        --brief \
        --close)

    dx-jobutil-add-output viewer_bookmark "$viewer_bookmark" --class=record
}
