#!/bin/bash

# Usage in systemd service file:
# ExecStart=/bin/bash -c eval "$(conda shell.bash hook)"; conda activate /opt/locus-focus/locusfocus_env && /opt/locus-focus/setup/start-celery.sh $ENV_FILE

# Get path to env file from first argument
ENV_FILE=$1

if [ -z "$ENV_FILE" ]; then
    echo "ERROR: no env file specified"
    exit 1
fi

if [ ! -e "$ENV_FILE" || ! -f "$ENV_FILE" ]; then
    echo "ERROR: env file $ENV_FILE does not exist"
    exit 1
fi

export $(grep -v '^#' $ENV_FILE | xargs)

# Start celery worker
celery -A make_celery.celery_app worker --loglevel=INFO
