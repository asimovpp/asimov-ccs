Contributions to ASiMoV-CCS are welcomed via pull requests.

When preparing a contribution we ask that you read the developers guide (see `dev_guide/`) to check for design and coding style intentions.
As described in the build tools documentation (`build_tools/build_system_readme.md`) the `flint` program is used to apply linting rules to the ASiMoV-CCS source.

We request that contributions are based on, and target, the `upstreaming` branch - this will help us
to upstream the contributions.
The `upstreaming` branch is kept in sync with the `develop` branch - ASiMoV-CCS developers should
ensure the `upstream` branch is up to date when reviewing any pull requests, and that the proposed
changes are compatible with this.
Once a contribution is accepted, the ASiMoV-CCS developer should merge the `upstreaming` branch into
`develop` on the main repository and this will be reflected in the public repository.

Before a contribution is merged it will be reviewed by one of the ASiMoV-CCS developers, in particular to check the coding standard is followed and to verify the full test suite must still pass.
