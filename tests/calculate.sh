#!/bin/bash
set -o nounset
set -o errexit
set -o pipefail

cd $1

find . -regex 'chord_pred.txt$' -exec wc -l {} \;
cat chord_pred.txt | cut -f 1-7 | md5sum
