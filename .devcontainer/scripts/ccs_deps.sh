#!/usr/bin/env bash

# Setting up ownership
sudo chown vscode:vscode /workspaces
cd /workspaces || exit

# Cloning ccs-dependencies
git clone -b python https://github.com/asimovpp/ccs-dependencies.git
cd ccs-dependencies || exit

# Building and install deps
poetry update
poetry install
poetry run ccs-install --env gnu_devcontainer

# fix up
eval $(poetry run ccs-setup --env gnu_devcontainer)
mv /workspaces/ccs-libs/rcmf90-gnu ${RCM_F90}