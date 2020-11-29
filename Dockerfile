# ToolFactory docker based on Bjoern's galaxy docker image
# galaxy needs sudo
FROM quay.io/bgruening/galaxy:20.09
MAINTAINER Ross Lazarus ross.lazarus@gmail.com
# most of the action moved to post-start-actions.sh
# MUST be copied to export root and made executable
ENV GALAXY_CONFIG_BRAND ToolFactory Docker
ENV HOST_DOCKER_GID $HOST_DOCKER_GID
# need to pass in HOST_DOCKER_GID using (e.g.) dockergid=`getent group docker | cut -d: -f3` from the run script
# the docker group inside the container is changed to match the host so members of the docker group
# can use the passed docker.sock
COPY files/hackadmin.py files/install-history.py files/install-deps.py files/galaxy_wait.py /usr/local/bin/
COPY files/startup.sh /usr/bin/startup
# uses [ -f post-start-actions.sh ] - original has an executable test and it fails even with shebang and ax+
COPY files/tfwelcome.html /etc/galaxy/web/welcome.html
# that seems to work.
COPY files/galaxy.yml files/tool_shed.yml files/tool_sheds_conf.xml  /etc/galaxy/
COPY files/galaxy.yml files/tool_shed.yml files/tool_sheds_conf.xml $GALAXY_ROOT/config/

RUN chmod a+x /usr/bin/startup \
  && groupadd docker \
  && usermod -G docker galaxy \
  && mkdir -p $GALAXY_ROOT/database/dependencies \
  && mkdir -p $GALAXY_ROOT/config/workflows \
  && mkdir -p $GALAXY_ROOT/config/tools \
  && mkdir -p $GALAXY_ROOT/config/histories \
  && mkdir -p $GALAXY_ROOT/database/dependencies \
  && mkdir -p $GALAXY_ROOT/tools/toolfactory \
  && /galaxy_venv/bin/python3 -m pip install --upgrade pip \
  && chown -R galaxy /home/galaxy /galaxy-central \
  && echo "galaxy ALL=(ALL:ALL) NOPASSWD: SETENV: /usr/bin/docker\n" >> /etc/sudoers \
  && apt update -y && apt install -y wget python3-venv python3-pip python3-dev gcc build-essential \
  && apt-get clean && apt-get purge \
  &&  rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
COPY files/TFhistory.tar.gz $GALAXY_ROOT/config/histories/TFhistory.tar.gz
ADD files/TFsample.ga $GALAXY_ROOT/config/workflows/tf.ga
COPY files/TFtools.yml $GALAXY_ROOT/config/tools/TFtools.yml


