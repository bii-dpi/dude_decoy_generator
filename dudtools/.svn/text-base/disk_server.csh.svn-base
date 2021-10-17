#!/bin/csh

set base=$0:h
if ( "$base" == "$0" ) then
    set base="."
endif
set dud = "$base"

set port = `grep 'DISK_PORT =' $dud/dudtools.py | awk '{print $3}'`
set today = `date +'%Y%m%d'`
$dud/disk_server.py 0.0.0.0 $port >&! disk_${today}_${port}.log
