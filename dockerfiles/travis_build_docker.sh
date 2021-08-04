#!/bin/bash

# Only run if main or tagged
if [[ "$TRAVIS_TAG" == "$TRAVIS_BRANCH" ]] || [[ "$TRAVIS_BRANCH" == "main" ]];
then
  # Log in to docker hub
  echo "$DOCKER_PASSWORD" | docker login -u "$DOCKER_USERNAME" --password-stdin &&

  if [[ ${#TRAVIS_TAG} -gt 0 ]];
  then
    TAG_ARRAY=("$TRAVIS_TAG" "latest")
  else
    TAG_ARRAY=("unstable")
  fi

  DOCKER_TAG=""
  for TAG_NAME in "${TAG_ARRAY[@]}"; do
    docker pull "$DOCKER_USERNAME"/pyopatra"$DOCKER_BUILD_NAME":"$TAG_NAME" || echo "No existing image for pyopatra$DOCKER_BUILD_NAME:$TAG_NAME"
    DOCKER_TAG="$DOCKER_TAG -t $DOCKER_USERNAME/pyopatra$DOCKER_BUILD_NAME:$TAG_NAME"
  done

  docker build $DOCKER_TAG -f dockerfiles/Dockerfile"$DOCKER_BUILD_NAME" .

  for TAG_NAME in ${TAG_ARRAY[@]}; do
    docker push "$DOCKER_USERNAME"/pyopatra"$DOCKER_BUILD_NAME":"$TAG_NAME"
  done
fi