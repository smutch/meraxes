#!/usr/bin/env bash

echo "If you haven't already, please be sure and run"
echo "  spack compiler find && spack external find"
echo
echo "This bootstrap script will create a new spack environment in ../.meraxes-spack_env"
echo

while true; do
    read -p "Do you wish to setup the spack env (../.meraxes-spack_env) with all requirements? [Y/n] " yn
    case $yn in
        [Yy]* ) break;;
        ''    ) break;;
        [Nn]* ) exit;;
        * ) echo "Please answer yes or no.";;
    esac
done

spack env create -d ../.meraxes-spack_env spack.yaml
eval `spack env activate --sh ../.meraxes-spack_env`
spack concretize
spack install
