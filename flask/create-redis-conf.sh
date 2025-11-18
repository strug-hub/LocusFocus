#!/bin/bash

# This script creates a redis.conf file with the password from the .env file
# And also binds to localhost, so that only servers on the same machine can connect

if [ -f .env ]; then
    if grep -q "REDIS_PASSWORD=" .env; then
        echo "Redis password found in .env file"
    else
        echo "No REDIS_PASSWORD found in .env file. Please add it."
        exit 1
    fi
    REDIS_PASSWORD=$(grep -oP '(?<=REDIS_PASSWORD=).*' .env)
    echo "requirepass $REDIS_PASSWORD" > redis/redis.conf
    echo "bind 127.0.0.1" >> redis/redis.conf
    echo "Password saved to redis/redis.conf"
else
    echo "No .env file found. Please create one and add REDIS_PASSWORD to it."
fi
