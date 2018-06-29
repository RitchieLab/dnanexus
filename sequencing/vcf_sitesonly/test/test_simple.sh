#!/bin/bash

##################################################################################
# Strategy:
#
# 1) Build App
#    Assumes directory structure is
#    ├── App Dir
#    │   ├── Test
#    │       ├── src/< app script >
#    │       ├── test
#    │           ├── < This Test Script >
#      
# 2) Run test using base_input dictionary
#    - Target a specific folder for all App output
#
# 3) Wait for job completion
#
# 4) Remove all contents in target folder
#
# 5) Exit 0 or non-0 code depening on Job status
#
# Todo:
#   # UNQ_APP_NAME not even close gaurantee to be unique, is the a guid/uuid pre-implemented in bash?
##################################################################################

# Parameters to fill in

TARGET_PROJECT="project-BzQf6k80V3bJk7x0yv6z82j7"  # DNAnexus Regression Testing Project AWS US east

# Some setup
set -x

SCRIPTPATH="$( cd "$(dirname "$0")" ; pwd -P )"
APP_DIR=$(dirname "${SCRIPTPATH}")

APP_NAME=$(jq '.name' <"${APP_DIR}/dxapp.json" | tr -d '"')
UNQ_APP_NAME="${APP_NAME}_${RANDOM}"

TARGET_PATH="${TARGET_PROJECT}:/${APP_NAME}/${APP_NAME}_${RANDOM}"


#  You need to modify this function to construct the correct inputs
construct-inputs () {
	inputs=()
	inputs+=(-ivariants_vcfgzs="file-F6PYf480BZxbgYQxK80FGqJz")
	#inputs+=(-imore_inputs="I am an Inputs")
}

build-applet () {
	echo "Building appdir: ${APP_DIR}" 1>&2
	dx mkdir -p "${TARGET_PATH}"
	retval=$(dx build "${APP_DIR}" -a --destination "${TARGET_PATH}/${UNQ_APP_NAME}" | jq '.id' | tr -d '"')
	echo "App ID: ${retval}" 1>&2
	echo "${retval}"
}

execute-applet () {
	app_id="${1}"
	run_name="${APP_NAME}_simple_run"
	construct-inputs  # Constructing input array, since it can't be passed as input to bash func

	retval=$(dx run "${app_id}" "${inputs[@]}" --name "${run_name}" -y --brief --destination "${TARGET_PATH}")
	echo "Running ${run_name}: ${retval}" 1>&2
	echo "${retval}"
}

verify-result () {
	job_id=$1
	dx wait "${job_id}" 1>&2
	retval=$?
	echo "${retval}"
}

clean-up () {
	dx rm -r "${TARGET_PATH}"
}

main () {
	app_id=$(build-applet)

	job_id=$(execute-applet "${app_id}")

	exit_code=$(verify-result "${job_id}")
	[[ $exit_code == 0 ]] && echo "${APP_NAME} SUCCESS" || echo "${APP_NAME} FAIL"

	clean-up
}

main
