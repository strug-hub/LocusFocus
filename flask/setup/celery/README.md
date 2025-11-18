# Celery

Celery is a job scheduler for Python that's used by LocusFocus to run colocalization analyses asynchronously.

Celery depends on Redis to handle message passing.

This directory contains scripts to start and stop the Redis server and the Celery worker.

## Requirements

- Docker (optional for Redis)

## Start Redis

To start the Redis server, run the `start-redis-docker.sh` script:

```bash
./start-redis-docker.sh
```

## Start Celery

To start the Celery worker, run the `start-celery.sh` script:

```bash
./start-celery.sh .env
```

The `.env` file is a file that contains environment variables that are used by Celery.
