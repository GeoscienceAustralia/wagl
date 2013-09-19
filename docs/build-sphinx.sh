#!/bin/bash

export MY_PATH=$(readlink -f ${0%/*})
export MY_SOURCE_PATH=$MY_PATH/source
echo $MY_SOURCE_PATH
sphinx-apidoc -f -o $MY_SOURCE_PATH $MY_PATH/../ULA3
sphinx-apidoc -f -o $MY_SOURCE_PATH $MY_PATH/../image_processor
rm $MY_SOURCE_PATH/ULA3.test*
find . -name "*.rst" | dos2unix
