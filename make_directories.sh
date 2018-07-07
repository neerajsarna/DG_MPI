#!/bin/sh

BASEDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
BASEDIR="$BASEDIR/$1"
echo "Base Directory of the present script: "
echo $BASEDIR

for i in {3..20}
do
		dirname="$BASEDIR/M$i"
		echo $dirname
 	        mkdir -p -- "$dirname"
done
