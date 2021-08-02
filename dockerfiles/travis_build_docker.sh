#!/bin/bash

# Log in to docker hub
echo "$DOCKER_PASSWORD" | docker login -u "$DOCKER_USERNAME" --password-stdin &&

# Build all images
BUILD_ARRAY=("" "-mvapich-ib" "-sandybridge-mvapich-ib")
for BUILD_NAME in ${BUILD_ARRAY[@]}; do
  if [[ ${#TRAVIS_TAG} -gt 0 ]];
  then
    DOCKER_TAG="-t $DOCKER_USERNAME/pyopatra$BUILD_NAME:$TRAVIS_TAG -t $DOCKER_USERNAME/pyopatra$BUILD_NAME:latest"
  else
    DOCKER_TAG="-t $DOCKER_USERNAME/pyopatra$BUILD_NAME:unstable"
  fi
  docker build $DOCKER_TAG -f dockerfiles/Dockerfile"$BUILD_NAME" . \
  && docker push "$DOCKER_USERNAME"/pyopatra"$BUILD_NAME":"$DOCKER_TAG"
done