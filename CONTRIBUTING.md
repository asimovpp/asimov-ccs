Contributions to ASiMoV-CCS are welcomed via pull requests.

When preparing a contribution we ask that you read the developers guide (see `dev_guide/`) to check for design and coding style intentions.
As described in the build tools documentation (`build_tools/build_system_readme.md`) the `flint` program is used to apply linting rules to the ASiMoV-CCS source.

Before a contribution is merged it will be reviewed by one of the ASiMoV-CCS developers, in particular to check the coding standard is followed and to verify the full test suite must still pass.

# Instructions for contributing on GitHub

If you have not done so already, create a fork of the ASiMoV-CCS repository where you can prepare your contributions.

We request that contributions are based on, and target, the `upstreaming` branch - this will help us
to upstream the contributions.
The `upstreaming` branch should be kept up to date with the `develop` branch to include the latest developments.

## Before starting a contribution

We appreciate people taking the time and effort to contribute to ASiMoV-CCS.
To prevent potential duplication of effort, we ask that before embarking on making code changes that you open an issue
on the GitHub repository outlining the proposed change and why it is necessary, it may be that someone is already
working on something similar or that there would be a conflict with some other planned developments, this will give the
opportunity to flag potential issues in advance.

## Preparing a contribution

Before beginning, make sure your `develop` and `upstreaming` branches are synchronised with the main repository, this
will ensure your work does not deviate from the rest of the code, easing the merge process.
Merge the `develop` branch into `upstreaming` to make sure it is up to date, and start a new branch from `upstreaming`,
commit and push work to your fork on your working branch.
We recommend pulling the `develop` branch and merging into your working branch somewhat regularly if the development of
your contribution continues for some time.

## Submitting a contribution

Once you are satisfied with your contribution, open a pull request against `upstreaming` on the main repository, linking
it to the issue where you proposed the change.
This will begin a process of further interaction with the ASiMoV-CCS developers, who may request changes before the
contribution is accepted.

Note that you don't need to have a "complete" contribution to open a pull request.
Opening a draft pull request early in the development process to begin discussions on code changes can also be beneficial.

# Instructions for core developers accepting contributions

The process for accepting a contribution into the upstream GitLab repository is outlined below.
Currently, a manual process is required to push the contribution to the upstream repository.

1. When a pull request is opened, check it is targetting the `upstreaming` branch, change the target branch if necessary.
2. Ensure that the `upstreaming` branch is up to date by merging `develop` into `upstreaming` and pushing to GitHub if 
necessary, this might raise merge conflicts with the pull request if it has drifted from upstream developments, work
with the contributor to resolve these.
3. Once satisfied with a contribution, accept the pull request into `upstreaming` on GitHub.
4. Pull `upstreaming` locally and merge into your local `develop` branch, if `upstreaming` was up to date there should
be no issues here.
5. Push `develop` to the upstream repository (this will then be mirrored back to GitHub)
6. Merge `develop` back into `upstreaming` and push back to GitHub, ensuring the branch is up to date for users - do
not push `upstreaming` to the upstream GitLab repository.
