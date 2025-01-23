#!/bin/bash
#
# write-config.sh extracts configuration information from GSM and writes it to a set of files
# in a directory. This simplifies access to the secrets from other scripts and applications.
#
# We want to use this in a gradle task, so it takes arguments both as command line
# options and as envvars. For automations, like GHA, the command line can be specified.
# For developer use, we can set our favorite envvars and let gradle ensure the directory is properly populated.
#
# The environment passed in is used to configure several other parameters, including the target for running
# the integration tests.
#
# For personal environments, we assume that the target name is the same as the personal namespace name.
# The output directory includes the following files:
#   ---------------------------+-------------------------------------------------------------------------
#   target.txt                 | the target that generated this set of config files. Allows the script
#                              | to skip regenerating the environment on a rerun.
#   ---------------------------+-------------------------------------------------------------------------
#   teaspoons-sa.json         | SA for running Teaspoons - this is taken from GSM for the provided target
#   ---------------------------+-------------------------------------------------------------------------

function usage {
  cat <<EOF
Usage: $0 [<target>] [<outputdir>] "

  <target> can be:
    local - for testing against a local server (bootRun)
    dev - uses secrets from the dev environment
    qa - uses secrets from the qa environment
    help or ? - print this help
    clean - removes all files from the output directory
    * - anything else is assumed to be a personal environment using the terra-kernel-k8s
  If <target> is not specified, then use the envvar TEASPOONS_WRITE_CONFIG
  If TEASPOONS_WRITE_CONFIG is not specified, then use local

  <outputdir> defaults to "../config/" relative to the script. When run from the gradle rootdir, it will be
  in the expected place for automation.

EOF
 exit 1
}

# Get the inputs with defaulting
script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )/.." &> /dev/null  && pwd )"
default_outputdir="${script_dir}/config"
default_target=${TEASPOONS_WRITE_CONFIG:-local}
target=${1:-$default_target}
outputdir=${2:-$default_outputdir}

case $target in
    help | ?)
        usage
        ;;

    clean)
        rm "${outputdir}"/* &> /dev/null
        exit 0
        ;;

    local)
        # for local development we will use the dev environment configuration because our app is currently set up to work
        # with the dev environment by default
        fcenv=dev
        ;;

    dev)
        fcenv=dev
        ;;

    qa)
        fcenv=qa
        ;;


    *) # personal env
        k8senv=integration
        namespace=$target
        fcenv=dev
        ;;
esac

# Create the output directory if it doesn't already exist
mkdir -p "${outputdir}"

# If there is a config and it matches, don't regenerate
if [ -e "${outputdir}/target.txt" ]; then
    oldtarget=$(<"${outputdir}/target.txt")
    if [ "$oldtarget" = "$target" ]; then
        echo "Config for $target already written"
        exit 0
    fi
fi

# read a secret out of google secrets manager
read_secret_gsm() {
  gcloud secrets versions access latest --project="$1" --secret="$2"
}

# grab teaspoons service account json from vault
read_secret_gsm "broad-dsde-${fcenv}" "teaspoons-sa-secret" > "${outputdir}/teaspoons-sa.json"

# We made it to the end, so record the target and avoid redos
echo "$target" > "${outputdir}/target.txt"
