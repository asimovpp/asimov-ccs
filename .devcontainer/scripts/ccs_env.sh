#!/usr/bin/env bash

export PATH=/workspaces/ccs-libs/makedepf90-gnu/bin:${PATH}
pushd /workspaces/ccs-dependencies || exit
source env.sh
popd
