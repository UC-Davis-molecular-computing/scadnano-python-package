# Contributing to the scadnano Python package 
First off, thanks for taking the time to contribute!

The following is a set of guidelines for contributing to scadnano.
Feel free to propose changes to this document in a pull request, 
or post questions as issues on the [issues page](https://github.com/UC-Davis-molecular-computing/scadnano/issues).







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














## Making Contributions


### Cloning

The first step is cloning the repository so you have it available locally.

```
git clone https://github.com/UC-Davis-molecular-computing/scadnano-python-package.git
```

Changes to the scadnano package should be pushed to the
[`dev`](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/tree/dev) branch. So switch to the `dev` branch:

```
git checkout dev
```











### Pushing to the repository dev branch and documenting changes (done on all updates)

Minor changes, such as updating README, adding example files, etc., can be committed directly to the `dev` branch. (Note: currently this option is only available to administrators; other contributors should follow the instructions below.)

For any more significant change that is made (e.g., closing an issue, adding a new feature), follow these steps:

1. If there is not already a GitHub issue describing the desired change, make one. Make sure that its title is a self-contained description, and that it describes the change we would like to make to the software. For example, *"problem with importing gridless design"* is a bad title. A better title is *"fix problem where importing gridless design with negative x coordinates throws exception"*.

2. Make a new branch specifically for the issue. Base this branch off of `dev` (**WARNING**: in GitHub desktop, the default is to base it off of `master`, so switch that). The title of the issue (with appropriate hyphenation) is a good name for the branch. (In GitHub Desktop, if you paste the title of the issue, it automatically adds the hyphens.)

3. If it is about fixing a bug, *first* add tests to reproduce the bug before working on fixing it. (This is so-called [test-driven development](https://www.google.com/search?q=test-driven+development))

4. If it is about implementing a feature, first add tests to test the feature. For instance, if you are adding a new method, this involves writing code that calls the method and tests various combinations of example inputs and expected output.

5. Work entirely in that branch to fix the issue.

6. Run unit tests and ensure they pass.

7. Commit the changes. In the commit message, reference the issue using the phrase "fixes #123" or "closes #123" (see [here](https://docs.github.com/en/enterprise/2.16/user/github/managing-your-work-on-github/closing-issues-using-keywords)). Also, in the commit message, describe the issue that was fixed (one easy way is to copy the title of the issue); this message will show up in automatically generated release notes, so this is part of the official documentation of what changed.

8. Create a pull request (PR). **WARNING:** by default, it will want to merge into the `master` branch. Change the destination branch to `dev`.

9. Wait for all checks to complete (see next section), and then merge the changes from the new branch into `dev`. This will typically require someone else to review the code first and possibly request changes.

10. After merging, it will say that the branch you just merged from can be safely deleted. Delete the branch.

11. Locally, remember to switch back to the `dev` branch and pull it. (Although you added those changes locally, they revert back once you switch to your local `dev` branch, which needs to be synced with the remote repo for you to see the changes that were just merged from the now-deleted temporary branch.)









### Pushing to the repository master branch and documenting changes (done less frequently)

Less frequently, pull requests (abbreviated PR) can be made from `dev` to `master`, but make sure that `dev` is working before merging to `master` as all changes to `master` are automatically built and deployed to [PyPI](https://pypi.org/project/scadnano/), which is the site hosting the pip installation package, and [readthedocs](https://scadnano-python-package.readthedocs.io/en/latest/), which is the site hosting the API documentation. That is, changes to master immediately affect users installing via pip or reading online documention, so it is critical that these work.

**WARNING:** Always wait for the checks to complete. This is important to ensure that unit tests pass. 

We have an automated release system (through a GitHub action) that automatically creates release notes when changes are merged into the master branch.

Although the GitHub web interface abbreviates long commit messages, the full commit message is included for each commit in a PR.

However, commit descriptions are not shown in the release notes. In GitHub desktop these are two separate fields; on the command line they appear to be indicated by two separate usages of the `-m` flag: https://stackoverflow.com/questions/16122234/how-to-commit-a-change-with-both-message-and-description-from-the-command-li.

So make sure that everything people should see in the automatically generated release notes is included in the commit message. (If not, then more manual editing of the release notes is required.) GitHub lets you [automatically close](https://docs.github.com/en/enterprise/2.16/user/github/managing-your-work-on-github/closing-issues-using-keywords) an issue by putting a phrase such as "closes #14". Although the release notes will link to the issue that was closed, they [will not describe it in any other way](https://github.com/marvinpinto/actions/issues/34). So it is important, for the sake of having readable release notes, to describe briefly the issue that was closed in the commit message.

One simple way to do this is to copy/paste the title of the issue into the commit message. For this reason, issue titles should be stated in terms of what change should happen to handle an issue. For example, instead of the title being *"3D position is improperly calculated from grid position"*, a better issue title is *"calculate 3D position correctly from grid position"*. That way, when the issue is fixed in a commit, that title can simply be copied and pasted as the description of what was done for the commit message. (But you should still add "fixes #<issue_number>" in the commit message, e.g., the full commit message could be *"fixes #101; calculate 3D position correctly from grid position"* .)

Users can read the description by clicking on the link to the commit or the pull request, but anything is put there, then the commit message should say something like "click on commit/PR for more details".

Breaking changes should be announced explicitly, perhaps in the commit message, but ideally also manually added at the top of the release notes, indicating what users need to do to deal with the change.

See here for an example: https://github.com/UC-Davis-molecular-computing/scadnano-python-package/releases/tag/v0.10.0

So the steps for committing to the master branch are:

1. If necessary, follow the instructions above to merge changes from a temporary branch to the `dev` branch. There will typically be several of these. Despite GitHub's suggestions to keep commit messages short and put longer text in descriptions, because only the commit message is included in the release notes, it's okay to put more detail in the message (but very long stuff should go in the description, or possibly documentation such as the README.md file).

    One of the changes committed should change the version number. We follow [semantic versioning](https://semver.org/). This is a string of the form `"MAJOR.MINOR.PATCH"`, e.g., `"0.9.3"`
    - For the web interface repo scadnano, this is located at the top of the file https://github.com/UC-Davis-molecular-computing/scadnano/blob/master/lib/src/constants.dart
    - For the Python library repo scadnano-python-package, this is located in two places: the bottom of the file https://github.com/UC-Davis-molecular-computing/scadnano-python-package/blob/master/scadnano/_version.py (as `__version__ = "0.9.3"` or something similar) and the near the top of the file https://github.com/UC-Davis-molecular-computing/scadnano-python-package/blob/master/scadnano/scadnano.py (as `__version__ = "0.9.3"` or something similar). This latter one is only there for users who do not install from PyPI, and who simply download the file scadnano.py to put it in a directory with their script).

    The PATCH version numbers are not always synced between the two repos, but, they should stay synced on MAJOR and MINOR versions. **Note:** right now this isn't quite true since MINOR versions deal with backwards-compatible feature additions, and some features are supported on one but not the other; e.g., modifications can be made in the Python package but not the web interface, and calculating helix rolls/positions from crossovers can be done in the web interface but not the Python package. But post-version-1.0.0, the major and minor versions of the  should be enforced.

3. Ensure all unit tests pass.

4. In the Python repo, ensure that the documentation is generated without errors. First, run `pip install sphinx sphinx_rtd_theme`. Then, from within the subfolder `doc`, run the command `make html` (or `make.bat html` on Windows), ensure there are no errors, and inspect the documentation it generates in the folder `build`.

5. Create a PR to merge changes from dev into master. 

6. One the PR is reviewed and approved, do the merge.

7. Once the PR changes are merged, a release will be automatically created here: https://github.com/UC-Davis-molecular-computing/scadnano/releases or https://github.com/UC-Davis-molecular-computing/scadnano-python-package/releases. It will have a title that is a placerholder, which is a reminder to change its title and tag. Each commit will be documented, with the commit message (but not description) included in the release notes.

8. Change *both* the title *and* tag to the version number with a `v` prepended, e.g., `v0.9.3`. It is imperative to change the tag before the next merge into master, or else the release (which defaults to the tag `latest`) will be overwritten.









## Styleguide

Follow the [Python style guide](https://www.python.org/dev/peps/pep-0008/), which should come along in most IDEs in the form of plugins and extensions. 

The line length should be configured to 110, as the style guide limit of 79 is a bit too restrictive.

