#!/bin/bash
sudo rm -rf export/*
sudo rm -rf export/.distribution_config
sudo cp post-start-actions.sh  export/
sudo cp tfwelcome.html export/
docker build -t quay.io/fubar2/toolfactory-galaxy-docker  .
./starttoolfactory.sh

