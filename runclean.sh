#!/bin/bash
sudo rm -rf export/*
sudo rm -rf export/.distribution_config
sudo cp files/post-start-actions.sh  export/
sudo cp files/tfwelcome.html export/
docker build -t quay.io/fubar2/toolfactory-galaxy-docker  .
./starttoolfactory.sh

