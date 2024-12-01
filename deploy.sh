#!/bin/bash

# method name
METHOD_NAME=absolute

# dockerhub id
DOCKERHUB_ID=phylyc
VERSION=1.6

# build and push together
docker buildx build --platform linux/amd64 -t ${DOCKERHUB_ID}/${METHOD_NAME}:${VERSION} --push .
