#!/bin/bash

docker run -d -p 6379:6379 --name locusfocus-redis --restart always redis
