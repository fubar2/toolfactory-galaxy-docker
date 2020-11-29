#!/bin/bash
sudo rm -rf /home/ross/rossgit/toolfactory_docker/export/*
sudo rm -rf /home/ross/rossgit/toolfactory_docker/export/.distribution_config
sudo cp post-start-actions.sh  /home/ross/rossgit/toolfactory_docker/export/
sudo cp tfwelcome.html export/
#sudo chmod a+x /home/ross/rossgit/toolfactory_docker/export/
docker build -t toolfactory_docker  .
#sudo cp post-start-actions.sh /home/ross/rossgit/toolfactory_docker/export/
#Cannot get that to work.
#sudo chmod a+x /home/ross/rossgit/toolfactory_docker/export/post-start-actions.sh
./starttoolfactory.sh

