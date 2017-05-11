#!/bin/bash
rm -f $2
mkfifo $2
(zcat $1 > $2; rm $2) &
exit 0
