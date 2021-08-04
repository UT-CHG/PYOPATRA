#!/bin/bash

if [[ "${{ env.GITHUB_TAG }}" =~ v[0-9]+\.[0-9]+\.[0-9]+ ]] || [[ "${{ env.GITHUB_TAG }}" == "main" ]];
then
  # Log in to docker hub
  echo "${{ secrets.DOCKER_PASSWORD }}" | docker login -u "${{ secrets.DOCKER_USER }}" --password-stdin &&

  if [[ "${{ env.GITHUB_TAG }}" =~ v[0-9]+\.[0-9]+\.[0-9]+ ]]
  then
    TAG_ARRAY=("${{ env.GITHUB_TAG }}" "latest")
  else
    TAG_ARRAY=("unstable")
  fi

  DOCKER_TAG=""
  for TAG_NAME in "${TAG_ARRAY[@]}"; do
    DOCKER_TAG="$DOCKER_TAG -t ${{ secrets.DOCKER_USER }}/pyopatra${{ matrix.build_name }}:$TAG_NAME"
  done

  docker build $DOCKER_TAG -f dockerfiles/Dockerfile"${{ matrix.build_name }}" .

  for TAG_NAME in ${TAG_ARRAY[@]}; do
    docker push "${{ secrets.DOCKER_USER }}"/pyopatra"${{ matrix.build_name }}":"${{ env.GITHUB_TAG }}"
  done
fi