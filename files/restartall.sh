#!/bin/sh
supervisorctl restart galaxy:
sh /galaxy-central/run_tool_shed.sh --stop-daemon
sh /galaxy-central/run_tool_shed.sh --daemon
