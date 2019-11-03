#!/bin/bash
/Users/sean/Desktop/kb_concoct/test_local/run_docker.sh run --rm -v /Users/sean/Desktop/kb_concoct/test_local/subjobs/$1/workdir:/kb/module/work -v /Users/sean/Desktop/kb_concoct/test_local/workdir/tmp:/kb/module/work/tmp $4 -e "SDK_CALLBACK_URL=$3" $2 async
