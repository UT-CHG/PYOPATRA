#!/bin/bash

# Log in to docker hub
echo "$DOCKER_PASSWORD" | docker login -u "$DOCKER_USERNAME" --password-stdin &&

# Build all images
BUILD_ARRAY=("" "-mvapich-ib" "-sandybridge-mvapich-ib")
for BUILD_NAME in "${BUILD_ARRAY[@]}"; do
  if [[ ${#TRAVIS_TAG} -gt 0 ]];
  then
#    DOCKER_TAG="-t $DOCKER_USERNAME/pyopatra$BUILD_NAME:$TRAVIS_TAG -t $DOCKER_USERNAME/pyopatra$BUILD_NAME:latest"
    TAG_ARRAY=("$TRAVIS_TAG" "latest")
  else
    TAG_ARRAY=("unstable")
#    DOCKER_TAG="-t $DOCKER_USERNAME/pyopatra$BUILD_NAME:unstable"
  fi

  DOCKER_TAG=""
  for TAG_NAME in ${TAG_ARRAY[@]}; do
    DOCKER_TAG="$DOCKER_TAG -t $DOCKER_USERNAME/pyopatra$BUILD_NAME:$TAG_NAME"
  done

  docker build $DOCKER_TAG -f dockerfiles/Dockerfile"$BUILD_NAME" .

  for TAG_NAME in ${TAG_ARRAY[@]}; do
    docker push "$DOCKER_USERNAME"/pyopatra"$BUILD_NAME":"$TAG_NAME"
  done
done