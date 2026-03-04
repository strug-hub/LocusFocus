#!/bin/bash
set -e

wget https://yanglab.westlake.edu.cn/software/smr/download/smr-1.3.2-linux-x86_64.zip
unzip smr-1.3.2-linux-x86_64.zip
mv smr-1.3.2-linux-x86_64/smr misc/
rm smr-1.3.2-linux-x86_64.zip
