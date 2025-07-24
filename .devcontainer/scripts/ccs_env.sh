#!/usr/bin/env bash

export PATH=/workspaces/ccs-libs/makedepf90-gnu/bin:${PATH}
pushd /workspaces/ccs-dependencies || exit
poetry run ccs-setup > env.sh
source env.sh
popd

export