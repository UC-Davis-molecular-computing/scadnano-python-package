# Contributing to the scadnano Python package 
First off, thanks for taking the time to contribute!

The following is a set of guidelines for contributing to scadnano.
Feel free to propose changes to this document in a pull request, 
or post questions as issues on the [issues page](https://github.com/UC-Davis-molecular-computing/scadnano/issues).



## Table of contents

* [What should I know before I get started?](#what-should-i-know-before-i-get-started)
* [Cloning](#cloning)
* [Pushing to the repository dev branch and documenting changes (done on all updates)](#pushing-to-the-repository-dev-branch-and-documenting-changes-done-on-all-updates)
* [Pushing to the repository main branch and documenting changes (done less frequently)](#pushing-to-the-repository-main-branch-and-documenting-changes-done-less-frequently)
* [Styleguide](#styleguide)



## What should I know before I get started?

### Python
First, read the [README](README.md) to familiarize yourself with the package from a user's perspective.

The scadnano Python package requires at least Python 3.7. See the [README for instructions](README.md#installation) for using Python 3.6 if you cannot upgrade to 3.7 or later for some reason. 

### What to install

Follow the [installation instructions](README.md#installation) to install the correct version of Python if you don't have it already.

It is actually unnecessary for you to install scadnano via pip, so you can skip that step. In developing, you will have a local version of the package that you run and modify.

I suggest using a powerful IDE such as [PyCharm](https://www.jetbrains.com/pycharm/download/download-thanks.html). [Visual Studio Code](https://code.visualstudio.com/) is also good with the right plugins. The scadnano Python package uses type hints, and these tools are very helpful in giving static analysis warnings about the code that may represent errors that will manifest at run time.

### Keeping the scadnano package simple for users to install

One goal is to make the package as easy to install as possible, even for users who have trouble installing scadnano via pip. For this reason, we have two self-imposed constraints:

1. There are minimal package dependencies. scadnano can be run in most circumstances with a standard Python 3.7 (or above) installation. (One exception is the package [xlwt](https://pypi.org/project/xlwt/), which is required to call the method [`Design.write_idt_plate_excel_file()`](https://scadnano-python-package.readthedocs.io/#scadnano.Design.write_idt_plate_excel_file).)

2. All the required code is in a single file, [scadnano.py](scadnano/scadnano.py). This is one reason an IDE will help, because navigating a large source code file is easier in an IDE.

These two constraints imply that a user who has trouble installing via pip can simply copy the file scadnano.py into their working directory (or in some directory on their `PYTHONPATH`) and import it as normal.


### git

We use [git](https://git-scm.com/docs/gittutorial) and [GitHub](https://guides.github.com/activities/hello-world/). You can use the command-line git, or a GUI such as [GitHub desktop](https://desktop.github.com/), which is very easy to use and supports the most common git commands, but it is not fully-featured, so you may want another [git client](https://www.google.com/search?q=git+client) if you prefer not to use the command line.















## Cloning

The first step is cloning the repository so you have it available locally.

```
git clone https://github.com/UC-Davis-molecular-computing/scadnano-python-package.git
```

Changes to the scadnano package should be pushed to the
[`dev`](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/tree/dev) branch.
`main` is the default branch, but it holds released versions only: every push to `main` publishes a
new release to PyPI (see
[the section on pushing to main](#pushing-to-the-repository-main-branch-and-documenting-changes-done-less-frequently)).
So after cloning, switch to the `dev` branch:

```
git checkout dev
```











## Pushing to the repository dev branch and documenting changes (done on all updates)

Minor changes, such as updating README, adding example files, etc., can be committed directly to the `dev` branch. (Note: currently this option is only available to administrators; other contributors should follow the instructions below.)

For any more significant change that is made (e.g., closing an issue, adding a new feature), follow these steps:

1. If there is not already a GitHub issue describing the desired change, make one. Make sure that its title is a self-contained description, and that it describes the change we would like to make to the software. For example, *"problem with importing gridless design"* is a bad title. A better title is *"fix problem where importing gridless design with negative x coordinates throws exception"*.

2. Make a new branch specifically for the issue, **based on `dev`**.

    The easiest way is the helper script in the root of the repository, which takes the issue number:

    ```
    .\branch-for-issue.ps1 123      # Windows PowerShell
    ./branch-for-issue.sh 123       # macOS / Linux / Git Bash
    ```

    It names the branch from the issue number and title, creates it from `dev`, links it to the
    issue (so it appears in the issue's "Development" section), and checks it out. It needs the
    [GitHub CLI](https://cli.github.com/): `winget install --id GitHub.cli -e` on Windows or
    `brew install gh` on macOS, then `gh auth login`.

    **Do not use the "Create a branch" link on the issue page.** It always bases the new branch on
    the default branch, `main`, with no option to choose `dev`, and there is no equivalent
    issue-aware command in GitHub Desktop. The script exists precisely to work around this.

    If you create the branch by hand instead, make sure `dev` is checked out first, since both
    `git checkout -b` and GitHub Desktop base a new branch on whatever is currently checked out.
    The title of the issue (with appropriate hyphenation) is a good name for the branch.

3. If it is about fixing a bug, *first* add tests to reproduce the bug before working on fixing it. (This is so-called [test-driven development](https://www.google.com/search?q=test-driven+development))

4. If it is about implementing a feature, first add tests to test the feature. For instance, if you are adding a new method, this involves writing code that calls the method and tests various combinations of example inputs and expected output.

5. Work entirely in that branch to fix the issue.

6. Run unit tests and ensure they pass.

7. Commit the changes. In the commit message, reference the issue using the phrase "fixes #123" or "closes #123" (see [here](https://docs.github.com/en/enterprise/2.16/user/github/managing-your-work-on-github/closing-issues-using-keywords)). Also, in the commit message, describe the issue that was fixed (one easy way is to copy the title of the issue); this message will show up in automatically generated release notes, so this is part of the official documentation of what changed.

8. Create a pull request (PR) into `dev`. GitHub will propose `main` as the base, because `main` is
   the default branch — but **you do not need to remember to change it**: a workflow retargets any
   newly opened PR from `main` to `dev` automatically and leaves a comment saying so. (If you ever
   genuinely mean to target `main`, just change the base back; the workflow only runs when a PR is
   first opened.) Branches created with `branch-for-issue` also record `dev` as their base, so
   `gh pr create` targets `dev` directly.

9. Wait for all checks to complete (see next section), and then merge the changes from the new branch into `dev`. This will typically require someone else to review the code first and possibly request changes.

10. After merging, it will say that the branch you just merged from can be safely deleted. Delete the branch.

11. Locally, remember to switch back to the `dev` branch and pull it. (Although you added those changes locally, they revert back once you switch to your local `dev` branch, which needs to be synced with the remote repo for you to see the changes that were just merged from the now-deleted temporary branch.)

### What happens to the issue when you merge to dev

When a commit saying "fixes #123" reaches `dev`, issue #123 gets the label
[`closed in dev`](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/labels/closed%20in%20dev)
and **stays open**. That is deliberate: the fix exists, but users install from PyPI, which is
published only from `main`, so the issue is not really resolved for anyone yet. The label is applied
automatically; you no longer need to remember to add it.

GitHub's own "fixes #123" auto-closing only acts on commits reaching the **default** branch, which
is `main`. That is the main reason `main` is kept as the default: it means merging a fix to `dev`
does not prematurely tell users the issue is done.

The issue is closed for good once the fix is actually released. GitHub does that itself: the
release merges those same commits into `main`, the default branch, where its "fixes #123" handling
applies. The release workflow then adds a comment linking the release and the PyPI package, and
removes the `closed in dev` label. Issues you close by hand are not touched by any of this.









## Pushing to the repository main branch and documenting changes (done less frequently)

Less frequently, pull requests (abbreviated PR) can be made from `dev` to `main`, but make sure that `dev` is working before merging to `main` as all changes to `main` are automatically built and deployed to [PyPI](https://pypi.org/project/scadnano/), which is the site hosting the pip installation package, and [readthedocs](https://scadnano-python-package.readthedocs.io/en/latest/), which is the site hosting the API documentation. That is, changes to main immediately affect users installing via pip or reading online documentation, so it is critical that these work.

**WARNING:** Always wait for the checks to complete. This is important to ensure that unit tests pass. 

The release PR is the one PR that really does target `main`. Since `main` is the default branch,
GitHub proposes it as the base, so there is nothing to change. The workflow that retargets PRs to
`dev` deliberately skips this one, recognizing it by its `dev` head branch *in this repository*
(a contributor's fork may also have a `dev` branch, and those PRs are retargeted normally).

**Every push to `main` is a release.** The release workflow reads `__version__` from
[scadnano/scadnano.py](scadnano/scadnano.py), creates the tag `v{version}`, creates a GitHub release,
and publishes to PyPI. So you **must bump `__version__` on `dev` before merging `dev` into `main`**.
If you forget, the release workflow fails with a red X and an explanatory message, and nothing is
tagged or published; bump the version on `dev` and merge again to recover. (Practically, this means
there is no such thing as a casual push to `main` — even a README typo fix rides along with a version
bump, or waits for the next release.)

We have an automated release system (through a GitHub action) that automatically creates the version
tag and release notes when changes are merged into the main branch, publishes to PyPI, and closes the
issues fixed by the release.

Although the GitHub web interface abbreviates long commit messages, the full commit message is included for each commit in a PR.

However, commit descriptions are not shown in the release notes. In GitHub desktop these are two separate fields; on the command line they appear to be indicated by two separate usages of the `-m` flag: https://stackoverflow.com/questions/16122234/how-to-commit-a-change-with-both-message-and-description-from-the-command-li.

So make sure that everything people should see in the automatically generated release notes is included in the commit message. (If not, then more manual editing of the release notes is required.) GitHub lets you [automatically close](https://docs.github.com/en/enterprise/2.16/user/github/managing-your-work-on-github/closing-issues-using-keywords) an issue by putting a phrase such as "closes #14". The release notes link to the issue that was closed, but do not describe it in any other way. So it is important, for the sake of having readable release notes, to describe briefly the issue that was closed in the commit message.

One simple way to do this is to copy/paste the title of the issue into the commit message. For this reason, issue titles should be stated in terms of what change should happen to handle an issue. For example, instead of the title being *"3D position is improperly calculated from grid position"*, a better issue title is *"calculate 3D position correctly from grid position"*. That way, when the issue is fixed in a commit, that title can simply be copied and pasted as the description of what was done for the commit message. (But you should still add "fixes #<issue_number>" in the commit message, e.g., the full commit message could be *"fixes #101; calculate 3D position correctly from grid position"* .)

Users can read the description by clicking on the link to the commit or the pull request, but anything is put there, then the commit message should say something like "click on commit/PR for more details".

Breaking changes should be announced explicitly, perhaps in the commit message, but ideally also manually added at the top of the release notes, indicating what users need to do to deal with the change.

See here for an example: https://github.com/UC-Davis-molecular-computing/scadnano-python-package/releases/tag/v0.10.0

So the steps for committing to the main branch are:

1. If necessary, follow the instructions above to merge changes from a temporary branch to the `dev` branch. There will typically be several of these. Despite GitHub's suggestions to keep commit messages short and put longer text in descriptions, because only the commit message is included in the release notes, it's okay to put more detail in the message (but very long stuff should go in the description, or possibly documentation such as the README.md file).

2. Bump the version number on `dev`. This is **required** before every merge into `main`; the release
   workflow fails if the version was not bumped. We follow [semantic versioning](https://semver.org/):
   a string of the form `"MAJOR.MINOR.PATCH"`, e.g., `"0.9.3"`. Bump PATCH for bug fixes only, and
   MINOR for backwards-compatible feature additions.
    - For the web interface repo scadnano, this is located at the top of the file https://github.com/UC-Davis-molecular-computing/scadnano/blob/main/lib/src/constants.dart
    - For the Python library repo scadnano-python-package, there is a single source of truth: the
      `__version__` line near the top of the file
      [scadnano/scadnano.py](scadnano/scadnano.py) (as `__version__ = "0.9.3"` or something similar).
      Keep the trailing `# version line; WARNING: ...` comment intact — `setup.py` finds this line by
      searching for that comment, and the release workflow reads the same line.

    The PATCH version numbers are not always synced between the two repos, but, they should stay synced on MAJOR and MINOR versions. **Note:** right now this isn't quite true since MINOR versions deal with backwards-compatible feature additions, and some features are supported on one but not the other; e.g., modifications can be made in the Python package but not the web interface, and calculating helix rolls/positions from crossovers can be done in the web interface but not the Python package. But post-version-1.0.0, the major and minor versions of the  should be enforced.

3. Ensure all unit tests pass.

4. In the Python repo, ensure that the documentation is generated without errors. First, run `pip install sphinx sphinx_rtd_theme`. This installs [Sphinx](https://www.sphinx-doc.org/en/main/), which is the most well-supported documentation generator for Python. (It's not very friendly, the syntax for things like links in docstrings is awkward, but it's well supported, so we use it.) Then, from within the subfolder `doc`, run the command `make html` (or `make.bat html` on Windows), ensure there are no errors, and inspect the documentation it generates in the folder `_build`.

5. Create a PR to merge changes from dev into main. `main` is the default base, so there is nothing
   to change here.

6. Once the PR is reviewed and approved, do the merge.

7. The merge triggers the release workflow, which automatically:
    - creates the tag `v{version}` at the merge commit — **no manual tagging or retagging is needed
      any more**;
    - creates a release here: https://github.com/UC-Davis-molecular-computing/scadnano-python-package/releases,
      titled `TODO: edit release notes for v{version}` (a deliberate placeholder reminding you to do
      step 8), whose body lists every commit since the previous release, with the commit message (but
      not description) included;
    - publishes the package to PyPI;
    - comments on every issue referenced with a closing keyword in those commits, with links to the
      release and the PyPI package, and removes their `closed in dev` labels. (The issues are
      *closed* by GitHub itself, since those commits reach `main`, the default branch.)

8. Edit the release: replace the placeholder title with the version number with a `v` prepended, e.g.,
   `v0.9.3`, and write a human summary of the release above the auto-generated `## Commits` list.
   Breaking changes belong at the top here. Saving this edit is also what regenerates `CHANGELOG.md`
   (see below).

9. Back-merge `main` into `dev`, so that `dev` picks up the release merge commit and the
   `CHANGELOG.md` commits.

### CHANGELOG.md

[CHANGELOG.md](CHANGELOG.md) is generated automatically from the GitHub releases; **never edit it by
hand**, since the next run overwrites it. Edit the release notes instead — saving an edit to a release
regenerates the whole file and commits it to `main`.

Because the file is regenerated from scratch every time, a missed update is never permanent: go to the
Actions tab, select the `changelog` workflow, and press "Run workflow". This is also the remedy for
editing a *very old* release (one tagged before this automation existed), which does not trigger the
workflow on its own.









## Styleguide

Follow the [Python style guide](https://www.python.org/dev/peps/pep-0008/), which should come along in most IDEs in the form of plugins and extensions. 

The line length should be configured to 120, as the style guide limit of 79 is a bit too restrictive.

