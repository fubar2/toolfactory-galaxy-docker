dockergid=`getent group docker | cut -d: -f3`
sudo docker run -d -p 8080:80 -p 8021:21 -p 8800:8800 -p 9009:9009 \
    --privileged=true -e GALAXY_DOCKER_ENABLED=True -e HOST_DOCKER_LEGACY=True \
    -v /home/ross/rossgit/toolfactory_docker/export:/export/ -e HOST_DOCKER_GID=$dockergid\
    -v /var/run/docker.sock:/var/run/docker.sock \
    toolfactory_docker
