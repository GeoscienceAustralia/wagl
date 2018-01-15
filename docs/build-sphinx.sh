#!/bin/bash

export MY_PATH=$(readlink -f ${0%/*})
export MY_SOURCE_PATH=$MY_PATH/source
echo $MY_SOURCE_PATH
sphinx-apidoc -f -o $MY_SOURCE_PATH $MY_PATH/../wagl
rm $MY_SOURCE_PATH/wagl.test*
find . -name "*.rst" | dos2unix
