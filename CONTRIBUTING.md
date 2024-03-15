Contributions to ASiMoV-CCS are welcomed via pull requests.

When preparing a contribution we ask that you read the developers guide (see `dev_guide/`) to check
for design and coding style intentions.
As described in the build tools documentation (`build_tools/build_system_readme.md`) the `flint`
program is used to apply linting rules to the ASiMoV-CCS source.

Before a contribution is merged it will be reviewed by one of the ASiMoV-CCS developers, in
particular to check the coding standard is followed and to verify the full test suite still
passes.

# Instructions for contributing on GitHub

If you have not done so already, create a fork of the ASiMoV-CCS repository where you can prepare
your contributions.

The `develop` branch contains the latest version of the code, new contributions should therefore
start from this branch.

## Before starting a contribution

We appreciate people taking the time and effort to contribute to ASiMoV-CCS.
To prevent potential duplication of effort, we ask that before embarking on making code changes that
you open an issue on the GitHub repository outlining the proposed change and why it is necessary, it
may be that someone is already working on something similar or that there would be a conflict with
some other planned developments, this will give the opportunity to flag potential issues in advance.

## Create a new branch

Hopefully there are no clashes with the proposed contribution and ongoing work and the proposed work
can proceed.
In the issue proposing the contribution click the `Create a branch` link, this will open a popup
where you can name your branch and select your forked repository as the branch destination.
At this point you will be able to checkout your branch locally and begin working.

### Before you start

If you have already worked on CCS then you should ensure your local repository view of `develop` is
up to date before beginning new work.
To do so we recommend using `git rebase` to keep the git history clean, assuming that the main
ASiMoV-CCS repository is `origin` and your branch is called `my-branch`:
```
git checkout develop
git pull origin/develop
git checkout my-branch
git rebase develop
```

You are now ready to begin.

### Keeping up to date

We recommend periodically pulling the `origin/develop` branch and rebasing your working branch onto
`develop`.
If the development of your contribution continues for some time, this will ensure that it does not
stray too far from the eventual merge target.

## Submitting a contribution

Once you are satisfied with your contribution, open a pull request against `upstreaming` on the main
repository, linking it to the issue where you proposed the change - if you followed the above steps
to create a new branch this linking should occur automatically, except you may need to change the
target branch of the pull request on GitHub.
This will begin a process of further interaction with the ASiMoV-CCS developers, who may request
changes before the contribution is accepted.

Note that you don't need to have a "complete" contribution to open a pull request.
Opening a draft pull request early in the development process to begin discussions on code changes
can also be beneficial.

# Instructions for core developers accepting contributions

The process for accepting a contribution into the upstream GitLab repository is outlined below.
Currently, a manual process is required to push the contribution to the upstream repository.

0. When a potential contributor opens an issue regarding making a contribution:
- is the proposed contribution well-described? Ask for clarification if not. **TODO:** develop a
  template for bug reports and proposed changes.
- check whether the proposed contribution conflicts or overlaps with any ongoing or planned work -
  try to discuss this with the contributor, perhaps they could collaborate on any preexisting work?
1. When a pull request is opened check it targets the `upstreaming` branch, change the target
   branch if necessary.
2. Ensure that the `upstreaming` branch is up to date by rebasing and pushing to GitHub if
   necessary, this might raise merge conflicts with the pull request if it has drifted from upstream
   developments, work with the contributor to resolve these.
3. Once satisfied with a contribution, accept the pull request into `upstreaming` on GitHub. Note that
   additional care should be taken when contributions affect tests or the running of tests, as a policy
   two developers should approve such contributions.
4. Pull `upstreaming` locally and merge into your local `develop` branch, if `upstreaming` was up to 
   date there should be no issues here.
5. Push `develop` to the upstream repository (this will then be mirrored back to GitHub)
6. Rebase `upstreaming` onto the new `develop` and push back to GitHub, ensuring the branch is up to
   date for users - do not push `upstreaming` to the upstream GitLab repository.
