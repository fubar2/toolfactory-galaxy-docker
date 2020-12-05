#!/bin/bash
rm -rf /export/tfvm
OLDPATH=$PATH
PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
/galaxy_venv/bin/python3 -m venv /export/tfvm && /export/tfvm/bin/python3 -m pip install ephemeris bioblend requests planemo galaxyxml
PATH=/export/tfvm/bin:$PATH
/export/tfvm/bin/workflow-install --workflow_path $GALAXY_ROOT/config/workflows/ -g http://localhost:8080 -a $GALAXY_DEFAULT_ADMIN_KEY
/export/tfvm/bin/python3 /usr/local/bin/install-history.py -a $GALAXY_DEFAULT_ADMIN_KEY -i $GALAXY_ROOT/config/histories/TFhistory.tar.gz -g http://localhost:8080
/export/tfvm/bin/shed-tools install -g http://localhost:8080 -a $GALAXY_DEFAULT_ADMIN_KEY -t $GALAXY_ROOT/config/tools/TFtools.yml
/export/tfvm/bin/python3 /usr/local/bin/hackadmin.py -d $GALAXY_ROOT/database/community.sqlite
echo "Loaded ToolFactory workflow and history for admin"
PATH=$OLDPATH
/galaxy_venv/bin/python3 $GALAXY_ROOT/scripts/tool_shed/build_ts_whoosh_index.py -c $GALAXY_ROOT/config/tool_shed.yml --config-section tool_shed
chown -R galaxy /home/galaxy/*
chown -R galaxy /galaxy-central/*
chown -R galaxy /galaxy-central/database/
# needed because it's a soft link....
mv /export/post-start-actions.sh /export/post-start-actions-has-been-run.sh
