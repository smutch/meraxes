#!/usr/bin/env bash

echo "If you haven't already, please be sure and run"
echo "  spack compiler find && spack external find"
echo

while true; do
    read -p "Do you wish to setup the spack env with all requirements? [Y/n] " yn
    case $yn in
        [Yy]* ) break;;
        ''    ) break;;
        [Nn]* ) exit;;
        * ) echo "Please answer yes or no.";;
    esac
done

cd ..
spack env create -d ./spack-env spack.yaml
eval `spack env activate --sh ./spack-env`
spack concretize
spack install
