
#----------------------------------------------------------------------
# run_app.sh is called to launch the server in either 
# interactive (developer) or daemon (production) mode
# ---------------------------------------------------------------------
# this script should be sourced by a parent script
# that collected these variables from user:
#   MODE            either interactive|dev|daemon|prod
#   ACTIONS_PATH    where this script can be found
#   DATA_SCRIPTS_PATH  where q master and master strain table are found
#   PIPELINE_NAME   the name of this pipeline
#   DATA_PATH       where relevant data are found (projects, samples...)
#   GENOMES_DIR     where genome annotation information can be found
#   DEV_PORT        http port using in development
#   PROD_PORT       http port using in production
#----------------------------------------------------------------------

# set a default mode to developer
if [ "$MODE" = "" ]; then export MODE=dev; fi

# check required variables
if [ "$ACTIONS_PATH" = "" ]; then
    echo "missing variable: ACTIONS_PATH"
    exit 1
fi
if [ "$DATA_SCRIPTS_PATH" = "" ]; then
    echo "missing variable: DATA_SCRIPTS_PATH"
    exit 1
fi
if [ "$PIPELINE_NAME" = "" ]; then
    echo "missing variable: PIPELINE_NAME"
    exit 1
fi
if [ "$DATA_PATH" = "" ]; then
    echo "missing variable: DATA_PATH"
    exit 1
fi
if [ "$GENOMES_DIR" = "" ]; then
    echo "missing variable: GENOMES_DIR"
    exit 1
fi

# launch app with output written to interactive shell
# used when developing app
if [[ "$MODE" == "interactive" || "$MODE" == "dev" ]]; then
    Rscript $ACTIONS_PATH/app.R

# launch app running as low priority background service
# used in production
elif [[ "$MODE" == "daemon" || "$MODE" == "prod" ]]; then
    rm -f nohup.out
    nohup nice Rscript $ACTIONS_PATH/app.R &
    echo $! > daemon.pid # for use by kill_daemon.sh

# report errors 
else
    echo
    echo "missing or unknown server mode"
    echo
    echo "USAGE: run_app.sh <interactive|dev|daemon|prod>"
    echo
    exit 1
fi

