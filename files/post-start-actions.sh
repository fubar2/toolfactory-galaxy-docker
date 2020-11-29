#!/bin/bash
rm -rf /export/tfvm
OLDPATH=$PATH
PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
/galaxy_venv/bin/python3 -m venv /export/tfvm && /export/tfvm/bin/python3 -m pip install ephemeris bioblend requests planemo galaxyxml
PATH=/export/tfvm/bin:$PATH
/export/tfvm/bin/workflow-install --workflow_path $GALAXY_ROOT/config/workflows/ -g http://localhost -a $GALAXY_DEFAULT_ADMIN_KEY
/export/tfvm/bin/python3 /usr/local/bin/install-history.py -a $GALAXY_DEFAULT_ADMIN_KEY -i $GALAXY_ROOT/config/histories/TFhistory.tar.gz -g http://localhost
/export/tfvm/bin/shed-tools install -g http://localhost -a $GALAXY_DEFAULT_ADMIN_KEY -t $GALAXY_ROOT/config/tools/TFtools.yml
/export/tfvm/bin/python3 /usr/local/bin/hackadmin.py -a $GALAXY_DEFAULT_ADMIN_KEY -g http://localhost -d $GALAXY_ROOT/database/community.sqlite
/export/tfvm/bin/python3 /usr/local/bin/install-deps.py -a $GALAXY_DEFAULT_ADMIN_KEY -t rgtf2 -g http://localhost
echo "Loaded ToolFactory workflow and history for admin"
PATH=$OLDPATH
/galaxy_venv/bin/python3 $GALAXY_ROOT/scripts/tool_shed/build_ts_whoosh_index.py -c config/tool_shed.yml --config-section tool_shed

