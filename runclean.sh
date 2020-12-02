#!/bin/bash
sudo rm -rf export/*
sudo rm -rf export/.distribution_config
sudo cp files/post-start-actions.sh  export/
sudo cp files/tfwelcome.html export/
# comment next line to use the image rather than building it
docker build -t quay.io/fubar2/toolfactory-galaxy-docker  .
./start.sh

