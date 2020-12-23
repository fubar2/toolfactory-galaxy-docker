# ToolFactory docker based on Bjoern's galaxy docker image
# quay.io/fubar2/toolfactory-galaxy-docker
FROM quay.io/bgruening/galaxy:20.09
MAINTAINER Ross Lazarus ross.lazarus@gmail.com
# most of the action moved to post-start-actions.sh
# MUST be copied to export root and made executable - needed to tweak the startup script test for that file from -x (executable) to -f
ENV GALAXY_CONFIG_BRAND ToolFactory Docker
# need to pass in HOST_DOCKER_GID using (e.g.) dockergid=`getent group docker | cut -d: -f3` from the run script
# so the host docker GID is passed in
# otherwise the host doesn't seem to be able to launch the planemo-biodocker container
ENV HOST_DOCKER_GID $HOST_DOCKER_GID
COPY files/hackadmin.py files/install-history.py files/install-deps.py files/galaxy_wait.py files/restartall.sh /usr/local/bin/
COPY files/startup.sh /usr/bin/startup
### starts toolshed and uses -f post-start-actions.sh because Bjoern's one has an executable test and it fails for some reason
COPY files/whooshcrontab /galaxy-central
COPY files/tfwelcome.html /etc/galaxy/web/welcome.html
COPY files/galaxy.yml files/tool_shed.yml files/tool_sheds_conf.xml  /etc/galaxy/
COPY files/galaxy.yml files/tool_shed.yml files/tool_sheds_conf.xml $GALAXY_ROOT/config/

RUN chmod -R a+x /usr/bin/startup \
  && groupadd docker \
  && usermod -G docker galaxy \
  && mkdir -p $GALAXY_ROOT/database/dependencies \
  && mkdir -p $GALAXY_ROOT/config/workflows \
  && mkdir -p $GALAXY_ROOT/config/tools \
  && mkdir -p $GALAXY_ROOT/config/histories \
  && mkdir -p $GALAXY_ROOT/database/dependencies \
  && mkdir -p $GALAXY_ROOT/tools/toolfactory \
  && /galaxy_venv/bin/python3 -m pip install --upgrade pip \
  && echo "galaxy ALL=(ALL:ALL) NOPASSWD: SETENV: /usr/bin/docker\n" >> /etc/sudoers \
  && apt update -y && apt upgrade -y && apt install -y wget python3-venv python3-pip python3-dev gcc fail2ban build-essential sqlite3 \
  && apt-get clean && apt-get purge \
  &&  rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* \
  && crontab -u root /galaxy-central/whooshcrontab
COPY files/TFhistory.tar.gz $GALAXY_ROOT/config/histories/TFhistory.tar.gz
ADD files/TFsample.ga $GALAXY_ROOT/config/workflows/tf.ga
COPY files/TFtools.yml $GALAXY_ROOT/config/tools/TFtools.yml


