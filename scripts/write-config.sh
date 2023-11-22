#!/bin/bash
#
# write-config.sh extracts configuration information from vault and writes it to a set of files
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
#   db-connection-name.txt     | Connection string for CloudSQL proxy
#   ---------------------------+-------------------------------------------------------------------------
#   db-name.txt                | Database name
#   ---------------------------+-------------------------------------------------------------------------
#   db-password.txt            | Database password
#   ---------------------------+-------------------------------------------------------------------------
#   db-username.txt            | Database username
#   ---------------------------+-------------------------------------------------------------------------
#   sqlproxy-sa.json           | SA of the CloudSQL proxy
#   ---------------------------+-------------------------------------------------------------------------
#   target.txt                 | the target that generated this set of config files. Allows the script
#                              | to skip regenerating the environment on a rerun.
#   ---------------------------+-------------------------------------------------------------------------
#   testrunner-sa.json         | SA for running TestRunner - this is always taken from integration/common
#   ---------------------------+-------------------------------------------------------------------------
#   testrunner-k8s-sa-token.txt| Credentials for TestRunner to manipulate the Kubernetes cluster under
#   testrunner-k8s-sa-key.txt  | test. Not all environments have this SA configured. If the k8env is 
#                              | integration and there is no configured SA, then the wsmtest one will be
#                              | retrieved. It won't work for all test runner tests.
#   ---------------------------+-------------------------------------------------------------------------
#   user-delegated-sa.json     | Firecloud SA used to masquerade as test users
#   ---------------------------+-------------------------------------------------------------------------
#   policy-sa.json             | SA that the policy service runs as
#   ---------------------------+-------------------------------------------------------------------------

function usage {
  cat <<EOF
Usage: $0 [<target>] [<vaulttoken>] [<outputdir>] [<vaultenv>]"

  <target> can be:
    local - for testing against a local server (bootRun)
    dev - uses secrets from the dev environment
    alpha - alpha test environment
    staging - release staging environment
    help or ? - print this help
    clean - removes all files from the output directory
    * - anything else is assumed to be a personal environment using the terra-kernel-k8s
  If <target> is not specified, then use the envvar TSPS_WRITE_CONFIG
  If TSPS_WRITE_CONFIG is not specified, then use local

  <vaulttoken> defaults to the token found in ~/.vault-token.

  <outputdir> defaults to "../config/" relative to the script. When run from the gradle rootdir, it will be
  in the expected place for automation.

  <vaultenv> can be:
    docker - run a docker image containing vault
    local  - run the vault installed locally
  If <vaultenv> is not specified, then use the envvar TSPS_VAULT_ENV
  If TSPS_VAULT_ENV is not specified, then we use docker
EOF
 exit 1
}

# Get the inputs with defaulting
script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )/.." &> /dev/null  && pwd )"
default_outputdir="${script_dir}/config"
default_target=${TSPS_WRITE_CONFIG:-local}
target=${1:-$default_target}
vaulttoken=${2:-$(cat "$HOME"/.vault-token)}
outputdir=${3:-$default_outputdir}
default_vaultenv=${TSPS_VAULT_ENV:-docker}
vaultenv=${4:-$default_vaultenv}

# The vault paths are irregular, so we map the target into separate variables:
# fcenv     - the firecloud delegated service account environment: dev, qa

case $target in
    help | ?)
        usage
        ;;

    clean)
        rm "${outputdir}"/* &> /dev/null
        exit 0
        ;;

    local)
        # for local development we will use the QA environemnt stuff for now to mimic our BEEs
        fcenv=qa
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

# Run vault either using a docker container or using the installed vault
function dovault {
    local dovaultpath=$1
    local dofilename=$2
    case $vaultenv in
        docker)
            docker run --rm -e VAULT_TOKEN="${vaulttoken}" broadinstitute/dsde-toolbox:dev \
                   vault read -format=json "${dovaultpath}" > "${dofilename}"
            ;;

        local)
            VAULT_TOKEN="${vaulttoken}" VAULT_ADDR="https://clotho.broadinstitute.org:8200" \
                   vault read -format=json "${dovaultpath}" > "${dofilename}"
            ;;
    esac
}

# Read a vault path into an output file, decoding from base64
# To detect missing tokens, we need to capture the docker result before
# doing the rest of the pipeline.
function vaultgetb64 {
    vaultpath=$1
    filename=$2
    fntmpfile=$(mktemp)
    dovault "${vaultpath}" "${fntmpfile}"
    result=$?
    if [ $result -ne 0 ]; then return $result; fi
    jq -r .data.key "${fntmpfile}" | base64 -d > "${filename}"
}

# Read a vault path into an output file
function vaultget {
    vaultpath=$1
    filename=$2
    fntmpfile=$(mktemp)
    dovault "${vaultpath}" "${fntmpfile}"
    result=$?
    if [ $result -ne 0 ]; then return $result; fi
    jq -r .data "${fntmpfile}" > "${filename}"
}

# Read database data from a vault path into a set of files
function vaultgetdb {
    vaultpath=$1
    fileprefix=$2
    fntmpfile=$(mktemp)
    dovault "${vaultpath}" "${fntmpfile}"
    result=$?
    if [ $result -ne 0 ]; then return $result; fi
    datafile=$(mktemp)
    jq -r .data "${fntmpfile}" > "${datafile}"
    jq -r '.db' "${datafile}" > "${outputdir}/${fileprefix}-name.txt"
    jq -r '.password' "${datafile}" > "${outputdir}/${fileprefix}-password.txt"
    jq -r '.username' "${datafile}" > "${outputdir}/${fileprefix}-username.txt"
}

# grab tsps service account json from vault
vaultget "secret/dsde/firecloud/${fcenv}/tsps/tsps-account.json" "${outputdir}/tsps-sa.json"

# We made it to the end, so record the target and avoid redos
echo "$target" > "${outputdir}/target.txt"
