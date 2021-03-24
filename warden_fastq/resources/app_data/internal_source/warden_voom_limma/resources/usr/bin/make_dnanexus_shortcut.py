#!/usr/bin/env python3

from __future__ import print_function

import os
import sys
import json
import logging
import argparse

from shutil import which
from subprocess import run, PIPE


_log = logging.getLogger(__name__)


def shell(statement):
    "Run a simple shell shell, returning whether the command completed successfully."
    _log.debug("Running command: '%s'", statement)
    return run(statement, shell=True, stdout=PIPE, stderr=PIPE)


def exit(statement, code):
    "Write an error to stderr then exit with code"
    sys.stderr.write(statement)
    sys.stderr.write("\n")
    sys.stderr.flush()
    sys.exit(code)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Creates a viewer shortcut object on DNAnexus platform "
                                                 "using the dx-toolkit.")
    parser.add_argument("-v", "--viewer", help="Name of the viewer in DNAnexus", required=True)
    parser.add_argument("-f", "--files", nargs="+", help="Files to include in the shortcut", required=True)
    parser.add_argument("-o", "--output", help="Name of bookmark output object", required=True)
    parser.add_argument("-p", "--project", help="Project to work within (default: the selected project in dx-toolkit).")
    parser.add_argument("--verbose", help="Verbosity for debugging", default=False, action="store_true")
    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(level=logging.DEBUG, format="[LOG] %(message)s")

    if shell("which dx").returncode != 0:
        exit("dx-toolkit not on PATH!", 1)
    _log.debug("dx-toolkit located.")

    if args.project:
        if shell("dx select %s" % (args.project)).returncode != 0:
            exit("Could not select project: %s!" % (args.project), 2)
        _log.debug("Selected project '%s'.", args.project)

    # get project id
    p = shell("dx env| grep 'project-' |cut -f2")
    if p.returncode != 0:
        exit("Could not retrieve project id!", 3)

    project_id = p.stdout.decode("utf-8").strip()
    _log.debug("Project ID: %s", project_id)

    # get viewer id
    p = shell("dx ls --brief %s" % (args.viewer))
    if p.returncode != 0:
        exit("Invalid fileviewer: %s" % (args.viewer), 4)

    viewer_id = p.stdout.decode("utf-8").strip()
    _log.debug("Viewer ID: %s", viewer_id)

    # get file ids
    file_ids = []
    for _file in args.files:
        p = shell("dx ls --brief %s" % (_file))
        if p.returncode != 0:
            exit("Invalid file: %s" % (_file), 5)

        this_file_id = p.stdout.decode("utf-8").strip()
        _log.debug("File ID for '%s': %s", _file, viewer_id)
        file_ids.append(project_id + ":" + this_file_id)

    mapper = {
        "preselectedIDs": file_ids,
        "fileViewer": {
            "project": project_id,
            "id": viewer_id
        }
    }

    mapper_json_str = json.dumps(mapper)
    _log.debug("JSON object: {}".format(mapper_json_str))

    # create a new record object
    p = shell("dx new record %s "
              "--details '%s' "
              "--type ViewerShortcut "
              "--type SJCloudVisualization "
              "--brief "
              "--close" % (args.output, mapper_json_str))

    if p.returncode != 0:
        exit("Could not create record object with name '%s'! " % (args.output), 6)
    _log.debug("Created record object.")

    record_id = p.stdout.decode("utf-8").strip()
    _log.debug("Record ID: %s", record_id)


# vim: ft=python:
