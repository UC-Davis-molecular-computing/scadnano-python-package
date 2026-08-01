# Changelog

This file is generated automatically from the project's
[GitHub releases](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/releases)
by [changelog-from-release](https://github.com/rhysd/changelog-from-release).
Do not edit it manually; edit the release notes instead.

<a id="v0.21.1"></a>
## [v0.21.1](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/releases/tag/v0.21.1) - 2026-08-01

No changes, just redoing some release-related actions and switching back to `main` being the default branch.

## Commits

- b2d09d5: Merge pull request [#340](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/issues/340) from UC-Davis-molecular-computing/dev (David Doty)
- 98ae57e: Merge pull request [#339](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/issues/339) from UC-Davis-molecular-computing/bump-version-0.21.1 (David Doty)
- a616426: bump version to 0.21.1 (David Doty)
- 10cba20: Merge pull request [#338](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/issues/338) from UC-Davis-molecular-computing/automate-dev-branch-workflow (David Doty)
- 029287f: automate the dev-branch workflow now that main is the default branch again (David Doty)
- 0a88677: Merge pull request [#337](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/issues/337) from UC-Davis-molecular-computing/dependabot-target-dev (David Doty)
- b3e1408: point dependabot PRs at dev instead of the default branch main (David Doty)
- 6c3b88a: Merge pull request [#336](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/issues/336) from UC-Davis-molecular-computing/main (David Doty)
- fb50af6: update CHANGELOG.md for release v0.21.0 (dave-doty)


[Changes][v0.21.1]


<a id="v0.21.0"></a>
## [v0.21.0](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/releases/tag/v0.21.0) - 2026-08-01

Main new feature is autostaple has been implemented directly in the package instead of farming it out to cadnano. Should work the same as before but be more reliable.

## Commits

- 1dcee7e: Merge pull request [#335](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/issues/335) from UC-Davis-molecular-computing/dev (David Doty)
- c730717: Merge pull request [#334](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/issues/334) from UC-Davis-molecular-computing/document-why-not-generate-notes (David Doty)
- 06e1e40: document why release notes are built with git log, not --generate-notes (David Doty)
- 2a85cd7: Merge pull request [#332](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/issues/332) from UC-Davis-molecular-computing/remove-test-automation-scratch (David Doty)
- 39a044c: remove the scratch file used to test the closed-in-dev automation (David Doty)
- 4fe113e: Merge pull request [#331](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/issues/331) from UC-Davis-molecular-computing/330-test-closed-in-dev-automation (David Doty)
- 0e8c512: fixes [#330](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/issues/330): test that merging to dev reopens the issue and labels it closed in dev (David Doty)
- b811192: Merge pull request [#329](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/issues/329) from UC-Davis-molecular-computing/make-release-workflow-rerunnable (David Doty)
- d27f007: make the release workflow re-runnable: skip already-published PyPI files (David Doty)
- d8c13a0: Merge pull request [#327](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/issues/327) from UC-Davis-molecular-computing/323-replace-automated-release-action (David Doty)
- 9ee4912: fixes [#323](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/issues/323): replace archived automated release action with gh CLI, add CHANGELOG.md automation and dev-branch issue lifecycle (David Doty)
- b2f6241: Merge pull request [#325](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/issues/325) from UC-Davis-molecular-computing/324-update-3rd-party-action-versions (David Doty)
- 22f1b20: Update run_unit_tests.yml (David Doty)
- 954fc75: revert computed matrix; use an explicit list ending in 3.x (David Doty)
- b0ee674: run unit tests on every currently supported Python, computed rather than hardcoded (David Doty)
- 4e418cb: update third-party action versions, and add dependabot to keep them fresh (David Doty)
- 160b0f2: Merge pull request [#322](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/issues/322) from UC-Davis-molecular-computing/321-add_half_crossover-removes-modifications (David Doty)
- a77d90d: CI unit tests: upgrade checkout and setup-python actions (David Doty)
- b1be729: CI unit tests: install dependencies with pip instead of conda (David Doty)
- ad7f9ec: fixed unit tests to not mix Helix's that have an explicit idx with those that don't (David Doty)
- bf3f9bb: added Python 3.14 to CI unit test action (David Doty)
- ed5acac: fixes [#321](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/issues/321) add_half_crossover removes modifications and fixes [#320](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/issues/320) ensure names preserved with add_half_crossover (David Doty)
- ce6ec2a: Update scadnano.py (David Doty)
- 708b2d4: fixed annoyance with design error messages where newline escape sequences were printed instead of newlines (David Doty)
- 4eccf6b: fixed unit tests that were not passing (David Doty)
- 8fd799e: fixed some type errors (David Doty)
- 0ecf8fb: updated logic of helix idx's to be -1 if not specified (David Doty)
- 7bcac21: formatted tests with Ruff (David Doty)
- ee11574: Update scadnano.py (David Doty)
- ec52afe: updated to modern type hint lowercase words (David Doty)
- acd8ded: ruff formatting with line-length 120 (David Doty)
- f3ea4a5: changed line length to 120 and will use Ruff in the future (David Doty)
- f199579: Fixes [#317](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/issues/317) - Add Autostaple ([#318](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/issues/318)) (beanbeanjuice)
- 00107bd: added Python 3.9 back to CI (David Doty)


[Changes][v0.21.0]


<a id="v0.20.1"></a>
## [v0.20.1](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/releases/tag/v0.20.1) - 2025-07-15

Added back `wc` (deprecated) so as not to break existing scripts.

## Commits
- 9cf18fe: changed some refs from `wc` to `rc` (David Doty) [#314](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/314)
- 5ec6b80: Update scadnano.py (David Doty) [#314](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/314)
- 85f153c: Update scadnano.py (David Doty) [#314](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/314)
- 5327c48: Update scadnano.py (David Doty) [#314](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/314)
- bd8df31: Update scadnano.py (David Doty) [#314](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/314)
- c2f87e6: minor (Dave Doty) [#314](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/314)
- 314e848: Update scadnano.py (David Doty) [#314](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/314)
- 05acbbf: added back `wc` (deprecated) and bumped version (David Doty) [#314](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/314)
- 4eba568: Update run_unit_tests.yml (David Doty) [#314](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/314)

[Changes][v0.20.1]


<a id="v0.20.0"></a>
## [v0.20.0](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/releases/tag/v0.20.0) - 2025-03-25

# Release notes

## **BREAKING CHANGE**: Minimum Python version is 3.9
The minimum Python version has been changed to 3.9. Technically there is not yet any part of the scadnano package that actually depends on 3.9 features; mainly this is done because support is being dropped for older Python versions in other tools, for example the unit tests that are run whenever we update the github repo no longer support older Python versions. 

As of version 0.20.0, you can still install scadnano with Python version 3.7 or 3.8 by using the option `--ignore-requires-python` with pip: `pip install scadnano --ignore-requires-python` https://stackoverflow.com/questions/75726452/can-i-force-the-install-of-a-package-that-requires-a-newer-python-version-than-t However, future versions of scadnano may use newer features of Python, so it is advisable to upgrade to the latest Python version now for maximum future compatibility.

## **BREAKING CHANGE**: parameter `only_strands_with_idt` changed to `only_strands_with_vendor_fields`
Version [0.19.0](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/releases/tag/v0.19.0) changed some names of fields and classes from IDT-specific names to the more generic "vendor fields". However,  the parameter `only_strands_with_idt` in a few methods was not at that time changed to `only_strands_with_vendor_fields`; now that parameter name has been changed as well in version 0.20.0.

## HelixGroups custom geometry
See also [UC-Davis-molecular-computing/scadnano#993 (comment)](https://github.com/UC-Davis-molecular-computing/scadnano/issues/993#issuecomment-2365320437).

There is a new optional field `HelixGroup.geometry`, which overrides the `Geometry` parameters of the `Design.geometry` field. For instance, the following code:

```python
import scadnano as sc


def create_design() -> sc.Design:
    group0 = sc.HelixGroup(grid=sc.square)
    group1 = sc.HelixGroup(grid=sc.square, 
                           # NOTE: here's the custom geometry for helix 1
                           geometry=sc.Geometry(bases_per_turn=18),
                           position=sc.Position3D(0, 3, 0))
    groups = {"group 0": group0, "group 1": group1}
    helices = [sc.Helix(idx=idx, max_offset=40, group=group) for idx, group in
               [(0, "group 0"), (1, "group 1")]]
    design = sc.Design(helices=helices, groups=groups, strands=[])
    design.draw_strand(0, 0).move(40)
    design.draw_strand(0, 40).move(-40)
    design.draw_strand(1, 0).move(40)
    design.draw_strand(1, 40).move(-40)

    return design


if __name__ == '__main__':
    d = create_design()
    d.write_scadnano_file(directory='output_designs')
```

produces this design with helix 0 having 10.5 base pairs per turn, and helix 1 having 18 base pairs per turn:

<img width="889" alt="image" src="https://github.com/user-attachments/assets/f50044b4-ddc2-4c05-af37-97d499c33efa">



## Commits
- 3600d1a: changed idt references to vendor_fields (David Doty) [#313](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/313)
- e8905fa: fixes [#304](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/issues/304): Typo in `set_dna_sequence` comment (David Doty) [#313](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/313)
- d9cd906: closes [#306](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/issues/306): allow per-helix-group geometry (David Doty) [#307](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/307)
- c597888: removed xlwt from environment.yml (David Doty) [#307](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/307)
- d6a52bd: Update run_unit_tests.yml (David Doty) [#307](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/307)
- efce159: fixed group-specific geometry in oxdna/oxview export (David Doty) [#313](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/313)
- 0c435bf: Update scadnano.py (David Doty) [#313](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/313)
- 2795055: fixed issue 309 by rotating helix position by group's orientation before adding to helix group's offset (DanielHader) [#312](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/312)
- cf75e42: Update scadnano.py (David Doty) [#312](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/312)
- b7afcd2: removed assignment of `None` to `Design.helices` and `Design.groups` since constructor sets them anyway (David Doty) [#313](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/313)
- abc8435: made `StrandBuilder` a dataclass and added documentation for fields (David Doty) [#313](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/313)
- c4f64ba: bumped minimum Python version to 3.9 and bumped version to 0.20.0 (David Doty) [#313](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/313)

[Changes][v0.20.0]


<a id="v0.19.4"></a>
## [v0.19.4](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/releases/tag/v0.19.4) - 2024-03-28

# Release notes

Changed type of `Design.to_oxview_format` from `dict` to `str` to keep convention with similar other methods such as `Design.to_oxdna_format`.

This is a breaking change, although the minor version number did not change. The method `Design.to_oxview_json` has the same functionality that `Design.to_oxview_format`, namely that it returns a JSON-serializable dict. Calling `json.dumps` on this dict produces the same output as `Design.to_oxview_format`.

[Changes][v0.19.4]


<a id="v0.19.3"></a>
## [v0.19.3](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/releases/tag/v0.19.3) - 2024-03-28

# Release notes

Mostly fixes this bug with oxView export: [#297](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/issues/297).

We jumped straight from v0.19.0 to v0.19.3 due to some testing that happened in the dev branch.


## Commits
- a78b395: bumped version (David Doty) [#299](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/299)
- 5cf25e3: changed parameter type of `Strand.set_domains` to allow `Extension`'s (David Doty) [#299](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/299)
- 1128beb: added the required build section (RayBipse) [#292](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/292)
- dcede52: Removed python.version (RayBipse) [#293](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/293)
- 3d20590: Removed python.version and added back tools.python (RayBipse) [#294](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/294)
- 4380042: Changed to use mamba and added python version requirement to environment.yml (RayBipse) [#295](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/295)
- 21aa34a: Added pip packages used to environment.yml (RayBipse) [#296](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/296)
- e5572ad: fixes [#297](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/issues/297): fix oxView export exception on Loopouts (David Doty) [#298](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/298)
- a220c32: Update scadnano.py (David Doty) [#299](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/299)

[Changes][v0.19.3]


<a id="v0.19.0"></a>
## [v0.19.0](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/releases/tag/v0.19.0) - 2023-09-07

# Release notes

## BREAKING CHANGE: Renamed IDT-specific fields
Some names related to the DNA synthesis company [IDT](https://www.idtdna.com/) have been renamed to be more general. You will have to rename these in your own code for it to run:

- class `IDTFields` --> `VendorFields`
- field `Modification.idt_text` --> `Modification.vendor_code`
- field `Strand.idt` --> `Strand.vendor_fields`

Additionally, some keys in the JSON format for `.sc` files have changed as well. Scadnano (both the web interface and the Python package) should be able to read files with the old keys and convert them to the new keys upon saving. However, if you are manually editing the `.sc` file then use the new keys.

Some IDT-specific methods remain, such as [`Design.write_idt_plate_excel_file`](https://scadnano-python-package.readthedocs.io/en/latest/#scadnano.Design.write_idt_plate_excel_file). These use the values in [`Strand.vendor_fields`](https://scadnano-python-package.readthedocs.io/en/latest/#scadnano.Strand.vendor_fields) and [`Modification.vendor_code`](https://scadnano-python-package.readthedocs.io/en/latest/#scadnano.Modification.vendor_code), but since the file format actually is specific to IDT, the method name is unchanged.

Although currently there are no methods for exporting to file formats recognized by other synthesis companies, in the meantime it should be straightforward to use the values in `VendorFields` to write custom code for such exports.

## Field `Modification.id` removed
Previously, scadnano used `Modification.id` as a unique identifier for modifications.

The field `id` has been removed. Now, the field `vendor_code` is used as a unique identifier for the modification, where `id` was previously used.

Previously, if no `id` was specified, but `vendor_code`/`idt_text` was, then `id` was set to the latter. Such code should continue to work unmodified. But code referencing `id` should now refer to `vendor_code` instead. Additionally, if a script used the same `vendor_code` for different Modifications, then this will break. Each Modification should now have a unique `vendor_code` field. 

Note that some vendors such as Eurofins use the same code for 5'/3' modifications (see issue [#283](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/issues/283)). For this reason, the modifications internally are stored in a way that encodes their location (5', 3', internal). But for any two modifications of the same "location" (both 5' modifications, both 3' modifications, or both internal modifications) the `vendor_code` must be unique to the type of modification.

## Commits
- 99e16d2: bumped version (David Doty) [#285](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/285)
- bd12d4d: Update scadnano.py (David Doty) [#285](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/285)
- dbc1058: Update scadnano.py (David Doty) [#285](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/285)
- ddc8eee: Update README.md (David Doty) [#285](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/285)
- 34211de: Update scadnano.py (David Doty) [#285](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/285)
- 3ed8158: Update scadnano.py (David Doty) [#285](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/285)
- 0bd7995: Update scadnano.py (David Doty) [#285](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/285)
- d012f17: Update scadnano.py (David Doty) [#285](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/285)
- 1c39fbc: closes [#279](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/issues/279): rename IDTFields to VendorFields (David Doty) [#280](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/280)
- 5fefd41: re-ran examples to alter changed name of IDT fields (David Doty) [#280](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/280)
- f679027: fixed docstrings (David Doty) [#280](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/280)
- 65630f4: fixed docstrings (David Doty) [#280](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/280)
- 3bb1373: Update scadnano.py (David Doty) [#285](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/285)
- 6356689: closes [#281](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/issues/281): remove field `Modification.id` (David Doty) [#282](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/282)
- 4c52f52: Update scadnano.py (David Doty) [#285](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/285)
- 032ebde: Update scadnano.py (David Doty) [#285](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/285)
- 4795224: Update scadnano.py (David Doty) [#285](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/285)
- b72a5c7: Update scadnano.py (David Doty) [#285](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/285)
- db9ec53: Update scadnano.py (David Doty) [#285](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/285)
- 41d93e3: removed `Position3D.clone()`, which is unnecessary since `Position3D` is frozen (David Doty) [#285](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/285)
- 6465434: Update scadnano.py (David Doty) [#285](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/285)
- fc30862: Update scadnano.py (David Doty) [#285](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/285)
- 37fa350: fixes [#283](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/issues/283): deal with non-unique Modification vendor codes (David Doty) [#284](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/284)
- 09ad91b: Update scadnano.py (David Doty) [#285](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/285)
- 82661bb: Update scadnano.py (David Doty) [#285](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/285)

[Changes][v0.19.0]


<a id="v0.18.3"></a>
## [v0.18.3](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/releases/tag/v0.18.3) - 2023-08-26

# Release notes

## Custom delimiter between DNA sequences of domains

A custom delimiter string to go in between DNA sequences of different domains on a strand can be specified, e.g., a space (which is ignored by IDT), so that the exported IDT sequences might look like this for a 3-domain strand:

ST0[8]0[15];/5Biosg/ ACGTCGT ACGTC ACGTAC;25nm;STD

See 

https://scadnano-python-package.readthedocs.io/en/latest/#scadnano.Strand.vendor_dna_sequence
https://scadnano-python-package.readthedocs.io/en/latest/#scadnano.Design.to_idt_bulk_input_format
https://scadnano-python-package.readthedocs.io/en/latest/#scadnano.Design.write_idt_bulk_input_file


## Commits
- 301edc4: bumped version (David Doty) [#278](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/278)
- 315c0f3: fixes [#276](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/issues/276): customize delimiter between domains in exported DNA sequences (David Doty) [#277](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/277)
- 5d8b61d: updated unit test to test for non-default delimiter (David Doty) [#277](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/277)
- 96b5e0d: removed delimiters between internal modifications and rest of sequence (David Doty) [#277](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/277)
- 93431b4: Update scadnano.py (David Doty) [#277](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/277)
- 19975dd: added unit test for internal modification that goes between bases (David Doty) [#278](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/278)

[Changes][v0.18.3]


<a id="v0.18.2"></a>
## [v0.18.2](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/releases/tag/v0.18.2) - 2023-08-24

Bug fixes related to relaxing helix rolls.

## Commits
- 8165951: bumped version (David Doty) [#275](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/275)
- 3f5bb6e: fixed one unit test (David Doty) [#275](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/275)
- a37b35a: Update scadnano.py (David Doty) [#275](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/275)
- 0837366: Update scadnano.py (David Doty) [#275](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/275)
- caf20e1: fixed bug in calculating minimum strain angle for helices with no crossovers (David Doty) [#275](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/275)
- e832e53: Update scadnano_tests.py (David Doty) [#275](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/275)
- 25fa92f: fixes [#268](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/issues/268): ignore loopouts when relaxing Helix rolls (David Doty) [#269](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/269)
- 9ea350b: fixed bug with relaxing roll starting from non-0 roll (David Doty) [#275](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/275)
- 7d938b5: cleaned up unit test (David Doty) [#275](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/275)
- 15c1514: Update scadnano.py (David Doty) [#275](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/275)
- 770366a: fixes [#270](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/issues/270): ignore "same helix" crossovers when relaxing helix rolls; also ensures domains on helix are stored in increasing order of start offsets (David Doty) [#272](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/272)
- 888565c: added new example (David Doty) [#272](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/272)
- 832c503: Update scadnano_tests.py (David Doty) [#275](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/275)
- d84419f: added unit test for relaxing helix rolls with intrahelix crossovers and fixed bug in relax code that didn't account for that case (David Doty) [#275](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/275)
- d3a6d7d: change parameter name (David Doty) [#275](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/275)
- a2d4e52: fixes [#273](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/issues/273): deal with inter-group crossovers when relaxing helix rolls (David Doty) [#274](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/274)
- 27568e4: Update scadnano_tests.py (David Doty) [#275](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/275)
- 1dedd8b: fixed relax helix rolls to ignore intergroup crossovers (David Doty) [#275](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/275)
- 51fbc9c: Update tutorial.md (David Doty) [#275](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/275)
- 5189e8d: Update scadnano.py (David Doty) [#275](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/275)

[Changes][v0.18.2]


<a id="v0.18.1"></a>
## [v0.18.1](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/releases/tag/v0.18.1) - 2023-08-16

# Release notes

## automatically set Helix rolls based on crossover locations ("relax" the rolls)
In this design, the crossovers are well placed relative to each other, but the helix rolls need to be rotated. As the slice bar shows, with the default roll of 0, the crossovers do not point at the helices they connect to:

![image](https://github.com/UC-Davis-molecular-computing/scadnano/assets/19274365/325b2b42-d7bd-4253-b949-37a7bfac23e9)

![image](https://github.com/UC-Davis-molecular-computing/scadnano/assets/19274365/45ae7174-0c4f-4055-95dc-e5ded60c84b2)

By calling [`Design.relax_helix_rolls`](https://scadnano-python-package.readthedocs.io/en/latest/#scadnano.Design.relax_helix_rolls), the helix rolls are set to minimize "strain" in the crossovers:

![image](https://github.com/UC-Davis-molecular-computing/scadnano/assets/19274365/59b9cdf8-5e68-4ca0-bd8a-0c007011a745)

![image](https://github.com/UC-Davis-molecular-computing/scadnano/assets/19274365/9d1c46de-0819-4d55-913a-fb3528ebde54)

Formally, the strain is defined as the square of the angular distance between the "ideal" angle for the crossover (i.e., the relative angle of the other helix to which it is connecting) and the crossover's current angle. The roll is set so as to minimize this total strain across all crossovers in each helix.

The reason to minimize the sum of the squared angular distances is that, if we model each crossover as an angular spring exerting a rotational force on the helix proportional to its displacement from the "ideal" angle (pointing directly at the other helix), this minimizes the total energy stored in the springs.

## Commits
- 2b11342: bumped version (David Doty) [#267](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/267)
- 49867e3: Merge branch 'dev' into 257-automatically-set-helix-rolls-based-on-crossover-locations-relax-the-rolls (David Doty) [#264](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/264)
- 123e001: closes [#257](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/issues/257): automatically set Helix rolls based on crossover locations (relax the rolls) (David Doty) [#264](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/264)
- 5208be3: fixed some unit tests (David Doty) [#264](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/264)
- c370a1f: Update scadnano.py (David Doty) [#264](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/264)
- c26dbe6: fixed dependencies (David Doty) [#266](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/266)
- 62c5dd5: changed unit test versions for openpyxl and tabulate to versions matching what is in Anaconda package repository (David Doty) [#266](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/266)
- 7b963c5: added instructions for installing openpyxl and tabulate to README (David Doty) [#266](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/266)
- f45d3f0: fixes [#234](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/issues/234): export new Excel format (David Doty) [#266](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/266)
- 3da4763: added openpyxl to tests_require (David Doty) [#266](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/266)
- c177ef0: added openpyxl installation to docs-check GitHub action (David Doty) [#266](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/266)
- 50899e3: Merge branch 'dev' into 257-automatically-set-helix-rolls-based-on-crossover-locations-relax-the-rolls (David Doty) [#264](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/264)

[Changes][v0.18.1]


<a id="v0.18.0"></a>
## [v0.18.0](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/releases/tag/v0.18.0) - 2023-07-22

# Release notes

## BREAKING CHANGE: label type is now `str`
Strands, domains, loopouts, and extensions have a field `label`. Previously the declared type was arbitrary, though at runtime it was required to be JSON-serializable.

Now we have changed the type of `label` field in `Strand`, `Domain`, `Loopout`, and `Extension` to `str` instead of an arbitrary object.

This is a breaking change because existing code using non-string labels will have to be altered to change the data to a string before storing and change it back to structured data when reading.

If you would like to store "structured data" (e.g., lists or dicts) in the label, you can serialize to a string and deserialize back to structured data manually using the `json` package.

Before, this was possible:

```python
from typing import List

# previously was possible, now is not supported
nums = [1, 2, 3]
strand.label = nums   # stores strand.label as the list [1, 2, 3]; would be a mypy type error now

# and to get the structured data back out:
nums: List[int] = strand.label  # would be a mypy type error now
```

Now this is necessary to store a list of `int`'s in the label:

```python
import json
from typing import List

nums = [1, 2, 3]
strand.label = json.dumps(nums)  # stores strand.label as the string '[1, 2, 3]'

# and to get the structured data back out:
nums: List[int] = json.loads(strand.label)  # nums is now the list [1, 2, 3]
```

## added p8634 variant of M13
There is a variant of M13 mentioned in a few papers (e.g., https://doi.org/10.1038/s41565-022-01283-1) called "p8634". It can be obtained from Tilibit (though not listed on their website). This sequence is now available as a predefined sequence. See https://scadnano-python-package.readthedocs.io/en/latest/#scadnano.M13Variant.p8634



## Commits
- d955eda: bumped version (David Doty) [#258](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/258)
- 2f58b2b: updated docstrings (David Doty) [#258](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/258)
- 0a482a7: Update scadnano.py (David Doty) [#258](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/258)
- c55411e: Update scadnano.py (David Doty) [#258](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/258)
- fca673a: added test generated by GitHub Copilot (David Doty) [#258](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/258)
- 3e84a38: changed docstring for `Design.base_pairs` (David Doty) [#258](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/258)
- fc2176b: added examples (David Doty) [#258](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/258)
- 7b0603b: changed to LR newlines (David Doty) [#258](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/258)
- 01c43e3: added p8634 variant of M13 (David Doty) [#258](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/258)
- 0678d55: formatting (David Doty) [#258](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/258)
- e8211db: fixed PyCharm warnings (David Doty) [#258](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/258)
- 602575d: closes [#261](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/issues/261): change `label` type to `str` (David Doty) [#262](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/262)
- a292757: removed all string type hints and replaced with forward references (not supported in Python 3.6) (David Doty) [#262](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/262)

[Changes][v0.18.0]


<a id="v0.17.7"></a>
## [v0.17.7](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/releases/tag/v0.17.7) - 2022-12-04

# Release notes

## new method [`Design.base_pairs()`](https://scadnano-python-package.readthedocs.io/en/latest/#scadnano.Design.base_pairs)
This method returns a dict mapping helix idx's where base pairs exist to a list of offsets where those base pairs are.

## Domain/Extension/Loopout colors
Individual `Domain`'s, `Extension`'s, and `Loopout`'s now have a `color` field. If specified, then this overrides the field `Strand.color` when displayed in the web interface. There is also a new method `StrandBuilder.with_domain_color()`:

```python
helices = [sc.Helix(max_offset=100) for _ in range(4)]
design = sc.Design(helices=helices, grid=sc.square)

red = sc.Color(255, 0, 0)
dark_red = sc.Color(150, 0, 0)
green = sc.Color(0, 255, 0)
dark_green = sc.Color(0, 150, 0)
blue = sc.Color(0, 0, 255)
dark_blue = sc.Color(0, 0, 150)
black = sc.Color(0, 0, 0)

design.draw_strand(0, 0) \
    .extension_5p(num_bases=5).with_domain_color(red) \
    .move(8).with_domain_color(green) \
    .loopout(1, 5).with_domain_color(dark_blue) \
    .move(-8).with_domain_color(dark_red) \
    .cross(2) \
    .move(8).with_domain_color(dark_green) \
    .cross(3) \
    .move(-8) \
    .extension_3p(num_bases=5).with_domain_color(black) \
    .with_color(blue)
```
![image](https://user-images.githubusercontent.com/19274365/205503932-8c2c33d0-a113-4d29-8918-2d6d6e804e8e.png)


## Commits
- 5586c9f: bumped version (David Doty) [#256](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/256)
- 6463865: closes [#252](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/issues/252): add method `Design.base_pairs()` (David Doty) [#254](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/254)
- fec6f43: changed return type of `Design.base_pairs()` to dict and added `allow_mismatches` parameter (David Doty) [#254](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/254)
- e795627: fixed unit tests after updating `with_sequence` and `with_domain_sequence` to default to not assigning complement (David Doty) [#254](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/254)
- 6b6001a: removed unused _popleft() function (David Doty) [#254](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/254)
- e55609d: Update scadnano.py (David Doty) [#256](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/256)
- afc44b8: don't put Helix idx in base_pairs dict if Helix has no base pairs (David Doty) [#256](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/256)
- 1049cdc: fixed base pair calculation to allow for deletions/insertions/wildcards/None if DNA sequence is not assigned (David Doty) [#256](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/256)
- e6c0737: added unit test for no base pairs (David Doty) [#256](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/256)
- f81b8ff: Update scadnano_tests.py (David Doty) [#256](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/256)
- 7b83a43: fixed bug in `Design.base_pairs()` when no reverse domains on helix (David Doty) [#256](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/256)
- 5312a59: Update scadnano_tests.py (David Doty) [#256](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/256)
- 1b7e605: closes [#238](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/issues/238): Domain colors (David Doty) [#255](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/255)
- fea3b32: added example with domain colors (David Doty) [#256](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/256)

[Changes][v0.17.7]


<a id="v0.17.6"></a>
## [v0.17.6](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/releases/tag/v0.17.6) - 2022-11-30

# Release notes

Two big updates:

1. support Extensions in oxDNA export ([#240](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/issues/240)); note it can have problems with large extensions between two close domains, see [#253](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/issues/253)
2. preliminary oxView export ([#173](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/issues/173))

## Commits
- 5cfa122: fixed docstrings for `Extension` (David Doty) [#250](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/250)
- 9b0d3e3: Update scadnano.py (David Doty) [#250](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/250)
- f84594e: Update scadnano.py (David Doty) [#250](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/250)
- 04a4e9d: Update scadnano.py (David Doty) [#250](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/250)
- 4ae1d93: formatting (David Doty) [#250](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/250)
- 357a59a: updated extensions example (David Doty) [#241](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/241)
- b477328: closes [#240](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/issues/240): support Extensions in oxDNA export (David Doty) [#241](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/241)
- 895469e: fix bug in assigning DNA sequence to strand bound to another strand with an extension (David Doty) [#250](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/250)
- 897f66d: fixed bug in oxDNA export when extension has no DNA sequence (David Doty) [#250](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/250)
- d1c38b9: Add add_extension to file creation (Constantine Evans) [#247](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/247)
- 1400ade: add test (Constantine Evans) [#247](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/247)
- d9160d6: Initial oxview export implementation (Constantine Evans) [#246](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/246)
- 176811f: Add documentation. (Constantine Evans) [#246](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/246)
- a823fee: Add tests. (Constantine Evans) [#246](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/246)
- 7086a55: probably-working, probably-bad bp implementation (Constantine Evans) [#246](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/246)
- 39b29c0: Tests for bp (Constantine Evans) [#246](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/246)
- 6f95338: add mismatch exclusion (Constantine Evans) [#246](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/246)
- ed71923: added method `Design.add_helix` (David Doty) [#250](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/250)
- d32537b: swapped order of `directory` and `filename` in `Design.write_scadnano_file` (David Doty) [#250](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/250)
- b0fee20: changed default display angle to 35 to match Web interface default (David Doty) [#250](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/250)
- 400cc3e: bumped version (David Doty) [#250](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/250)
- b9912ee: drop Python 3.6 on unit tests, add 3.10 and 3.11 (David Doty) [#250](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/250)
- 2c3b85f: adding quotes to python versions to try to fix pip install error on github actions (David Doty) [#250](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/250)
- e133772: fixed unit test to ensure helix max offset is sufficiently large (David Doty) [#250](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/250)
- 03d54ae: fixed unit test to ensure helix max offset is sufficiently large (David Doty) [#250](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/250)
- 931f22d: delete temp file manually in unit test for oxDNA export to avoid error where it cannot be read while open for writing (David Doty) [#250](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/250)
- 8e78237: added example of oxDNA export of extension and loopouts (David Doty) [#250](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/250)

[Changes][v0.17.6]


<a id="v0.17.5"></a>
## [v0.17.5](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/releases/tag/v0.17.5) - 2022-08-17

Several small bug fixes, plus one big new feature:

## Introducing DNA `Extension` API
This feature is supported on the web interface dev branch (https://scadnano.org/dev) at the time of this writing, and in the next release will be on the stable server (https://scadnano.org). See [UC-Davis-molecular-computing/scadnano#34 (comment)](https://github.com/UC-Davis-molecular-computing/scadnano/issues/34#issuecomment-1241023201)

DNA Extensions can now be represented. This gives a nice way to specify toeholds and other single-stranded extensions common in DNA strand displacement designs, which are single-stranded (unlike most Domains on a Helix) but are at the end of a strand (unlike a Loopout that is flanked by two Domains on a Strand).

They can be created using the `Extension` class. See full API here: https://scadnano-python-package.readthedocs.io/en/latest/#scadnano.Extension

```python
import scadnano as sc

domain = sc.Domain(helix=0, forward=True, start=0, end=10)
left_toehold = sc.Extension(num_bases=6)
right_toehold = sc.Extension(num_bases=5)
strand = sc.Strand([left_toehold, domain, right_toehold])
```

They can also be created with chained method calls (see https://scadnano-python-package.readthedocs.io/en/latest/#scadnano.StrandBuilder.extension_3p and https://scadnano-python-package.readthedocs.io/en/latest/#scadnano.StrandBuilder.extension_5p)

```python
import scadnano as sc

design = sc.Design(helices=[sc.Helix(max_offset=10)])
design.draw_strand(0,0).extension_5p(6).move(10).extension_3p(5)
```

See also documentation for this feature in the web interface: https://github.com/UC-Davis-molecular-computing/scadnano/releases/tag/v0.17.6

## Commits
- c8a8970: Update scadnano.py (David Doty) [#237](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/237)
- d482066: fixed docstring with string literal `'\n\n'` (David Doty) [#237](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/237)
- 03f1142: closes [#228](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/issues/228): option to export cadnano with no whitespace (David Doty) [#229](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/229)
- f25a28f: added `with_deletions` and `with_insertions` methods to `StrandBuilder` and used them in unit tests (David Doty) [#237](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/237)
- 00451a8: added check for deletion/insertion in range in `StrandBuilder.with_deletions` and `StrandBuilder.with_insertions` (David Doty) [#237](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/237)
- 4e60f9b: Add UT for chained methods for DNA extensions (Benjamin Lee) [#230](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/230)
- 78c716d: Add TestExportCadnanoV2.test_extension (Benjamin Lee) [#230](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/230)
- f16b629: Use self.design_6helix for extension chain method tests (Benjamin Lee) [#230](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/230)
- 2cdf023: Add move extension move test (Benjamin Lee) [#230](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/230)
- 28b0ded: Add extension ligate test (Benjamin Lee) [#230](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/230)
- a0170e9: Add crossover test (Benjamin Lee) [#230](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/230)
- e8725cc: test_add_full_crossover_on_extension_error (Benjamin Lee) [#230](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/230)
- 738979b: test_add_half_crossover_on_extension_ok (Benjamin Lee) [#230](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/230)
- 3615afb: test_add_half_crossover_on_extension_error (Benjamin Lee) [#230](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/230)
- 2feed8a: test_nick_on_extension (Benjamin Lee) [#230](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/230)
- 1433022: test_from_json_extension_design (Benjamin Lee) [#230](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/230)
- de44f47: test_to_json_extension_design (Benjamin Lee) [#230](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/230)
- e06ac96: Add "relative_offset" field (Benjamin Lee) [#230](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/230)
- d2f1afb: test_strand__with_relative_offset (Benjamin Lee) [#230](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/230)
- faa0a33: test_strand__with_relative_offset_on_domain_error and test_strand__with_relative_offset_on_empty_error (Benjamin Lee) [#230](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/230)
- a1601b7: test_strand__3p_extension_forward_default_relative_offset (Benjamin Lee) [#230](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/230)
- daad677: test_strand__5p_extension_forward_default_relative_offset (Benjamin Lee) [#230](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/230)
- 4c459ea: test_strand__3p_extension_reverse_default_relative_offset (Benjamin Lee) [#230](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/230)
- 0941030: test_strand__5p_extension_reverse_default_relative_offset (Benjamin Lee) [#230](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/230)
- 8f48617: Implement 3' extension (Benjamin Lee) [#230](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/230)
- 5574331: Implement 5p extension (Benjamin Lee) [#230](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/230)
- 8b6c060: Raise StrandError in as_circular if strand contains Extension (Benjamin Lee) [#230](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/230)
- 04b25df: Add 5p_extension test case for circular strand (Benjamin Lee) [#230](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/230)
- 6a2c101: Start docstring for extension (Benjamin Lee) [#230](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/230)
- fdbdf1c: Begin replacing extension -> extension_3p and add extension_5p_length argument to draw_strand (Benjamin Lee) [#230](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/230)
- 3048af8: Redo test_strand__3p_extension (Benjamin Lee) [#230](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/230)
- 0e4c694: Redo test_strand__5p_extension (Benjamin Lee) [#230](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/230)
- b8077c6: delete default relative_offset tests (Benjamin Lee) [#230](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/230)
- 122445c: In to, check for 3' extension; fixed some tests (Benjamin Lee) [#230](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/230)
- c45b916: Check for extension in cross (Benjamin Lee) [#230](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/230)
- cd8edaf: Fix test_strand__cross_after_3p_extension_should_raise_error (Benjamin Lee) [#230](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/230)
- 68466c8: Pass test_strand__extension_3p_after_loopout_should_raise_error (Benjamin Lee) [#230](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/230)
- de3175c: Pass test_strand__extension_3p_after_extension_should_raise_error (Benjamin Lee) [#230](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/230)
- 173fc34: Pass test_strand__update_to_after_3p_extension_should_raise_error (Benjamin Lee) [#230](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/230)
- 5610cf8: Pass test_strand__extension_3p_on_circular_strand_should_raise_error (Benjamin Lee) [#230](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/230)
- 185d097: Extension name, label, and sequence (Benjamin Lee) [#230](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/230)
- 43e06bf: Remove relative_offset and ExtensionBuilder (Benjamin Lee) [#230](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/230)
- 9b9e401: Fix some json and cadnano tests (Benjamin Lee) [#230](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/230)
- 819d0d8: Implement Extension from_json, and fix bug with default color from json (Benjamin Lee) [#230](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/230)
- 860be06: Implement Extension.to_json_serializable (Benjamin Lee) [#230](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/230)
- 06fe222: Fix nick,ligate,and crossover tests (Benjamin Lee) [#230](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/230)
- 2b485be: Add ligate and half crossover error cases (Benjamin Lee) [#230](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/230)
- 63fbd53: Fix middle domain bug in ligate (Benjamin Lee) [#230](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/230)
- 4a2023a: Handle ligate and crossover error case (Benjamin Lee) [#230](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/230)
- a27bcc7: Merge branch 'dev' into 2-support-loopouts-on-the-end-of-a-strand (Benjamin Lee) [#230](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/230)
- 75e4740: fixed errors in plate map documentation (David Doty) [#237](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/237)
- 996e59b: Update tutorial.md (David Doty) [#237](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/237)
- bb152cd: Add docstring for Extension (Benjamin Lee) [#230](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/230)
- d52f6ce: Add docstrings for rest of Extension fields and functions (Benjamin Lee) [#230](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/230)
- 25b3c56: fixes [#231](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/issues/231): Design.add_full_crossover and Design.add_half_crossover should check for existing crossovers (David Doty) [#232](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/232)
- 33b027a: added unit tests for add_half_crossover (David Doty) [#232](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/232)
- 290a8db: added example of consecutive crossovers to examples/ directory and added Address information to error message when calling `Design.add_full_crossover` (David Doty) [#237](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/237)
- 0f61dd4: added example of how to assign `StrandBuilder` instance to a variable to use loops to specify a long strand with its methods (David Doty) [#237](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/237)
- c24fbb7: Update README.md (David Doty) [#237](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/237)
- da0b087: added link to method `Design.draw_strand` in README when discussing `StrandBuilder` (David Doty) [#237](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/237)
- 8b52060: Set Position3D as frozen (Constantine Evans) [#233](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/233)
- 63da326: Update README.md (Dave Doty) [#237](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/237)
- 3d1b025: Update scadnano.py (Dave Doty) [#237](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/237)
- a158621: fixed some documentation and added constants for `display_angle` and `display_length` for `Extension` (David Doty) [#237](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/237)
- 73cd06c: added extensions example (David Doty) [#237](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/237)
- 380321e: updated extensions example (David Doty) [#237](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/237)
- f23b5f8: fixed oxDNA export bug that made length-0 normal vectors (David Doty) [#237](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/237)
- 9c93c1a: added no_git examples subfolder for examples that I don't want on the git repo (David Doty) [#237](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/237)
- 36e1725: fixed issue where get_normal_vector_to was being calculated incorrectly and loopout normal vectors weren't being normalized (Daniel Hader) [#236](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/236)
- 62f626d: bumped version (David Doty) [#237](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/237)
- bd63724: bumped version (David Doty) [#237](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/237)
- b1ff5b8: bumped version (David Doty) [#237](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/237)
- 20f50d5: fixed unit tests (David Doty) [#237](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/237)
- badd6b4: fixed last unit test (David Doty) [#237](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/237)

[Changes][v0.17.5]


<a id="v0.17.3"></a>
## [v0.17.3](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/releases/tag/v0.17.3) - 2022-03-29

# Release notes


## allow other table formats besides Markdown in `Design.plate_maps`

In a jupyter notebook cell, put this code:

```python
import scadnano as sc

helices = [sc.Helix(max_offset=100)]
design = sc.Design(helices=helices, strands=[], grid=sc.square)
design.draw_strand(0, 0).move(10).with_name('strand 0').with_idt(plate='plate 1', well='A1')
design.draw_strand(0, 10).move(10).with_name('strand 1').with_idt(plate='plate 1', well='A2')
design.draw_strand(0, 20).move(10).with_name('strand 2').with_idt(plate='plate 1', well='B2')
design.draw_strand(0, 30).move(10).with_name('strand 3').with_idt(plate='plate 1', well='B3')
design.draw_strand(0, 40).move(10).with_name('strand 4').with_idt(plate='plate 1', well='D7')

from IPython.display import display, Markdown
def dm(o):
    display(Markdown(o))

plate_map = design.plate_maps()[0]
dm(plate_map.to_table(tablefmt='html', vertical_borders=True))
```

It should render

<img width=500 src=https://user-images.githubusercontent.com/19274365/160621269-905ef7a0-051f-4593-b421-ddb55bf5c535.png></img>

The returned HTML uses inline styles to ensure there are vertical borders between columns of the table.             The vertical borders make it easier to see which column a well is in. This is useful when rendering in a Jupyter notebook, since the inline styles will be preserved when saving the Jupyter notebook using the nbconvert tool: https://nbconvert.readthedocs.io/en/latest/ 

Any format supported by the `tabular` package is supported as `tablefmt` for the method `PlateMap.to_table()`. See API for more details: https://scadnano-python-package.readthedocs.io/en/latest/#scadnano.PlateMap.to_table

## allow `Design.plate_maps` parameter `well_marker` to be function of well position

This code (note that strands do not require a name if using `well_marker`)

```python
import scadnano as sc

helices = [sc.Helix(max_offset=100)]
design = sc.Design(helices=helices, strands=[], grid=sc.square)
design.draw_strand(0, 0).move(10).with_idt(plate='plate 1', well='A1')
design.draw_strand(0, 10).move(10).with_idt(plate='plate 1', well='A2')
design.draw_strand(0, 20).move(10).with_idt(plate='plate 1', well='B2')
design.draw_strand(0, 30).move(10).with_idt(plate='plate 1', well='B3')
design.draw_strand(0, 40).move(10).with_idt(plate='plate 1', well='D7')

plate_map = design.plate_maps()[0]
print(plate_map.to_table(well_marker=lambda x:x))
```

prints

```

### plate "plate 1"
|     | 1   | 2   | 3   | 4   | 5   | 6   | 7   | 8   | 9   | 10   | 11   | 12   |
|:----|:----|:----|:----|:----|:----|:----|:----|:----|:----|:-----|:-----|:-----|
| A   | A1  | A2  |     |     |     |     |     |     |     |      |      |      |
| B   |     | B2  | B3  |     |     |     |     |     |     |      |      |      |
| C   |     |     |     |     |     |     |     |     |     |      |      |      |
| D   |     |     |     |     |     |     | D7  |     |     |      |      |      |
| E   |     |     |     |     |     |     |     |     |     |      |      |      |
| F   |     |     |     |     |     |     |     |     |     |      |      |      |
| G   |     |     |     |     |     |     |     |     |     |      |      |      |
| H   |     |     |     |     |     |     |     |     |     |      |      |      |
```



## Commits
- 1efd9ad: Closes [#149](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/issues/149); store DNA sequence in domains and loopouts, not in strand ([#223](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/issues/223)) (Benjamin Lee) [#223](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/223)
- 9e1a287: corrected docstring for `Design.draw_strand` (David Doty) [#226](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/226)
- 243562d: Update scadnano.py (David Doty) [#226](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/226)
- 7ed1610: Update tutorial.md (David Doty) [#226](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/226)
- f9e10b6: Update tutorial.md (David Doty) [#226](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/226)
- 45bf3ed: Update tutorial.md (David Doty) [#226](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/226)
- 48f6b54: Update tutorial.md (David Doty) [#226](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/226)
- 6610ef3: Update tutorial.md (David Doty) [#226](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/226)
- d6aa2ba: Update tutorial.md (David Doty) [#226](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/226)
- 0addd6d: Update tutorial.md (David Doty) [#226](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/226)
- 9ab30af: Update tutorial.md (David Doty) [#226](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/226)
- b09632f: Update scadnano.py (David Doty) [#226](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/226)
- 3fd2a86: Update scadnano.py (David Doty) [#226](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/226)
- e9ab4bc: Update scadnano.py (David Doty) [#226](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/226)
- 7d9764b: added documentation about idt_dna_sequence in docstring for `Strand.dna_sequence` property. (David Doty) [#225](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/225)
- 1150ea0: closes [#224](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/issues/224): allow other table formats besides Markdown in `Design.plate_maps` and closes [#222](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/issues/222): allow `Design.plate_maps` parameter `well_marker` to be function of well position (David Doty) [#225](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/225)
- 9ec66bb: added tabulate as dependency to tests (David Doty) [#225](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/225)
- 08bbf35: Update scadnano.py (David Doty) [#226](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/226)
- 2fb3d35: bumped version (David Doty) [#226](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/226)

[Changes][v0.17.3]


<a id="v0.17.2"></a>
## [v0.17.2](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/releases/tag/v0.17.2) - 2022-03-06

# Release notes

## Change modification connector length

This mainly affects how modifications are viewed in the web interface.

You can now change the length of the modification connector. The connector is drawn like a hydrocarbon chain with a default of 4 links:

![image](https://user-images.githubusercontent.com/19274365/151944535-37739dbc-7a61-4d59-90b3-7e713bf0a325.png)

The base class `Modification` now has a `connector_length` field: https://scadnano-python-package.readthedocs.io/en/latest/#scadnano.Modification.connector_length. Setting it will modify the displayed connector length in the web interface.


This is especially helpful to display modifications close to each other on a helix without their text overlapping:

![image](https://user-images.githubusercontent.com/19274365/151944658-93bf0187-a583-4ea6-a380-e9cbb734855d.png)




## Plate maps

You can now print a plate map of all or some of the strands in the design that have `Strand.idt.plate` and `Strand.idt.well` specified.

For example, the following code:

```python
import scadnano as sc

helices = [sc.Helix(max_offset=100)]
design = sc.Design(helices=helices, strands=[], grid=sc.square)
design.draw_strand(0, 0).move(10).with_name('strand 0').with_idt(plate='plate 1', well='A1')
design.draw_strand(0, 10).move(10).with_name('strand 1').with_idt(plate='plate 1', well='A2')
design.draw_strand(0, 20).move(10).with_name('strand 2').with_idt(plate='plate 1', well='B2')
design.draw_strand(0, 30).move(10).with_name('strand 3').with_idt(plate='plate 1', well='B3')
design.draw_strand(0, 40).move(10).with_name('strand 4').with_idt(plate='plate 1', well='D7')

# plate_maps = design.plate_maps_markdown(well_marker='X', strands=[design.strands[0], design.strands[3]])
plate_maps = design.plate_maps_markdown()
plate_map = plate_maps['plate 1']
print(plate_map)
```

prints the following Markdown representation of the plate:

```
## plate 1
|     | 1        | 2        | 3        | 4   | 5   | 6   | 7        | 8   | 9   | 10   | 11   | 12   |
|-----|----------|----------|----------|-----|-----|-----|----------|-----|-----|------|------|------|
| A   | strand 0 | strand 1 |          |     |     |     |          |     |     |      |      |      |
| B   |          | strand 2 | strand 3 |     |     |     |          |     |     |      |      |      |
| C   |          |          |          |     |     |     |          |     |     |      |      |      |
| D   |          |          |          |     |     |     | strand 4 |     |     |      |      |      |
| E   |          |          |          |     |     |     |          |     |     |      |      |      |
| F   |          |          |          |     |     |     |          |     |     |      |      |      |
| G   |          |          |          |     |     |     |          |     |     |      |      |      |
| H   |          |          |          |     |     |     |          |     |     |      |      |      |
```

One can specify only a subset of strands, and use a different entry than the strand's name:

```python
plate_maps = design.plate_maps_markdown(well_marker='X', strands=[design.strands[0], design.strands[3]])
```

which prints
```
## plate 1
|     | 1   | 2   | 3   | 4   | 5   | 6   | 7   | 8   | 9   | 10   | 11   | 12   |
|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|------|------|------|
| A   | X   |     |     |     |     |     |     |     |     |      |      |      |
| B   |     |     | X   |     |     |     |     |     |     |      |      |      |
| C   |     |     |     |     |     |     |     |     |     |      |      |      |
| D   |     |     |     |     |     |     |     |     |     |      |      |      |
| E   |     |     |     |     |     |     |     |     |     |      |      |      |
| F   |     |     |     |     |     |     |     |     |     |      |      |      |
| G   |     |     |     |     |     |     |     |     |     |      |      |      |
| H   |     |     |     |     |     |     |     |     |     |      |      |      |
```



## method `Design.strand` changed to `Design.draw_strand`
The method `Design.strand` is not well-named. It has been changed to `Design.draw_strand`. The method `Design.strand` still exists for now but is deprecated. It prints a warning that it is deprecated and then simply calls `Design.draw_strand`.

## Commits
- 5051e87: added some detail about unique IDs to Modification docstring (David Doty) [#221](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/221)
- dde8d5c: Update scadnano.py (David Doty) [#221](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/221)
- 020631d: fixes [#212](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/issues/212): change length of modification "connector" (David Doty) [#213](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/213)
- ca2522d: added parameter warn_duplicate_strand_names to scadnano and oxdna export (David Doty) [#221](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/221)
- 35d0e56: Update scadnano.py (David Doty) [#221](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/221)
- 13d11f5: Solving the bug pointed by test_crossover_to_same_helix.json (Tristan Stérin) [#214](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/214)
- cce904c: stap_loop -> stapLoop (Tristan Stérin) [#215](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/215)
- 6fe6447: Update scadnano.py (David Doty) [#221](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/221)
- 2c88d95: fixed bug in JSON serialization of IDT fields and made PlateCoordinate public (though still undocumented) (David Doty) [#221](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/221)
- ddd841d: closes [#210](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/issues/210): generate plate map (David Doty) [#216](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/216)
- ea9357a: added strands and well_markers parameters to plate_maps_markdown to enable a subset of strands to be put in the plate map, and to enable the table to contain another marker besides the strand name such as an X (David Doty) [#216](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/216)
- 1a8aa53: updated docstrings for plate_maps_markdown (David Doty) [#221](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/221)
- b330d3c: Bug is fixed (Tristan Stérin) [#217](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/217)
- 6170f98: Bug is fixed (Tristan Stérin) [#217](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/217)
- 9d9538d: closes [#219](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/issues/219): change `Design.strand` to `Design.draw_strand` (David Doty) [#220](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/220)
- 11702b2: added back method `Design.strand` (David Doty) [#220](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/220)
- 4f8c28c: removed several docstring references to deprecated method `Design.strand` (David Doty) [#220](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/220)

[Changes][v0.17.2]


<a id="v0.16.3"></a>
## [v0.16.3](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/releases/tag/v0.16.3) - 2021-08-03

## fixed bug in oxDNA export

We fixed a bug ([#188](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/issues/188)) in oxDNA export that made the helix roll go the opposite rotational direction from where it should have been.

## handling paranemic crossovers in cadnano export

The previous cadnano export did not correctly handle [paranemic crossovers](https://doi.org/10.1021/acs.chemrev.8b00207) ([#185](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/issues/185)): crossovers that join domains going in the same direction:

![image](https://user-images.githubusercontent.com/19274365/128058365-d52857b7-7379-44f2-9642-f061f4b5d439.png)


## Commits
- ce5163a: bumped version (David Doty) [#191](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/191)
- c83501f: changed .conf to .dat file extension in oxDNA export (David Doty) [#191](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/191)
- a5b3ad1: Update scadnano.py (David Doty) [#191](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/191)
- 309a8c1: Update scadnano.py (David Doty) [#191](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/191)
- 686abfc: fixed bug with file writing methods that allow use to specify optional filename with no extension (David Doty) [#191](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/191)
- 337fd1b: Bug fix helix number was badly set in cadnano export (Tristan Stérin) [#189](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/189)
- a474796: Adding cases in _cadnano_v2_place_crossover for paranemic crossovers (Tristan Stérin) [#189](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/189)
- 865b947: Un-ignoring test file (Tristan Stérin) [#189](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/189)
- 18aee03: Un-ignoring test file (Tristan Stérin) [#189](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/189)
- e128a3a: Merge branch 'dev' into bug_paranemic_crossovers (Cosmo Stérin) [#189](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/189)
- bf1a5e7: fixed sign error causing roll to be applied in the wrong direction (Daniel Hader) [#190](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/190)
- 447c999: fixes [#188](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/issues/188) (Daniel Hader) [#190](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/190)

[Changes][v0.16.3]


<a id="v0.16.2"></a>
## [v0.16.2](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/releases/tag/v0.16.2) - 2021-07-30

Fixes bugs in oxDNA export that would cause exported files not to be simulatable on oxDNA (due to non-cubic bounding box) and that, even with a cubic bounding box, would cause stable structures to melt.

## Commits
- f01be0f: corrected bad forward vector in oxDNA export that caused DNA origami to melt upon relaxation in oxDNA (Daniel Hader) [#183](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/183)
- 7bee4b1: updated oxDNA export to use cubic bounding box (taking all three sides to be max of computed bounding box) (David Doty) [#187](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/187)
- a4af148: bumped version (David Doty) [#187](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/187)
- cf18744: added Python 3.9 to list of versions unit tested on GitHub CI action (David Doty) [#187](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/187)

[Changes][v0.16.2]


<a id="v0.16.1"></a>
## [v0.16.1](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/releases/tag/v0.16.1) - 2021-07-18

New features:

## export to oxDNA
See [#139](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/issues/139). There is now a method [Design.write_oxdna_files()](https://scadnano-python-package.readthedocs.io/en/latest/#scadnano.Design.write_oxdna_files) that can write files representing the scadnano Design readable by the tool [oxDNA](https://oxdna.org/) and the web app [oxView](https://sulcgroup.github.io/oxdna-viewer/).

[oxView](https://sulcgroup.github.io/oxdna-viewer/) can be used to visualize the intended 3D structure of the design, with base positions inferred based on Helix (x,y,z) position and (pitch, yaw, roll) orientation angles. [oxDNA](https://oxdna.org/) can be used to simulate physical stress and motion to predict an expected 3D structure.

You can also call [Design.to_oxdna_format()](https://scadnano-python-package.readthedocs.io/en/latest/#scadnano.Design.to_oxdna_format) to get the two strings that are written to the two oxDNA files when calling [Design.write_oxdna_files()](https://scadnano-python-package.readthedocs.io/en/latest/#scadnano.Design.write_oxdna_files).

## export to IDT plates rebalances strands among last two plates if the last plate has too few strands 
See [#177](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/issues/177). Calling [Design.write_idt_plate_excel_file()](https://scadnano-python-package.readthedocs.io/en/latest/#scadnano.Design.write_idt_plate_excel_file) will now rebalance strands in the last two plates if necessary to ensure that each plate has the minimum required not to be charged extra by IDT. (24 strands for a 96-well plate, and 96 strands for a 384-well plate).

For instance if there are 202 strands, previously the three 96-well plates would have had 96, 96, 10 strands respectively. Now they will have 96, 82, 24, respectively, moving 14 strands from the second plate to the third to ensure all have at least 24 strands.


## Commits
- 312e761: Implemented oxdna export, but still need to add unit tests (Benjamin Lee) [#179](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/179)
- 37d4795: fixed bug where headers of conf file where not written (untested) (David Doty) [#179](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/179)
- a77e695: Bug fix and cleanup of oxdna conversion functions and very basic unit test case (Daniel Hader) [#179](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/179)
- d59c22d: added note to top of README that relative links won't work on PyPI site (David Doty) [#182](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/182)
- 72137fd: Added Oxdna unit tests and additional comments (Anelise Cho) [#179](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/179)
- 4cb98a2: Fixed OxDNA export issue where insertions would throw error, finished OxDNA export loopout unit test (Anelise Cho) [#179](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/179)
- 3d62c92: Update scadnano_tests.py (Anelise Cho) [#179](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/179)
- 49a259e: fixed Sphinx docstring error due to *'s being interpreted as emphasis (David Doty) [#179](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/179)
- e6a9c8d: updated variable names in oxDNA export code (David Doty) [#182](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/182)
- 1937acd: removed some PEP and mypy warnings and bumped version (David Doty) [#182](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/182)
- 52a2a5e: fix purification/scale typo in to_idt_bulk_input_format (Constantine Evans) [#180](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/180)
- c9401f7: closes [#177](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/issues/177); export to IDT plates rebalances strands among last two plates if the last plate has too few strands (David Doty) [#181](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/181)
- 92699a7: updated IDT Excel export API documentation to discuss rebalancing among last two plates to ensure minimum strand count (David Doty) [#181](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/181)
- 1f0d2c5: Update CONTRIBUTING.md (David Doty) [#182](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/182)
- 7f8f336: Update CONTRIBUTING.md (David Doty) [#182](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/182)
- c6dcd00: Update CONTRIBUTING.md (David Doty) [#182](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/182)
- 4c8b0e0: Update CONTRIBUTING.md (David Doty) [#182](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/182)
- f996834: closes [#139](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/issues/139) (not this commit actually, but I want this to show up in the release notes) (David Doty) [#182](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/182)

[Changes][v0.16.1]


<a id="v0.16.0"></a>
## [v0.16.0](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/releases/tag/v0.16.0) - 2021-05-23

## Breaking Changes

### `Helix.pitch` and `Helix.yaw` have been removed

The fields `Helix.pitch` and `Helix.yaw` have been removed to ensure that all helices within a `HelixGroup` are parallel. 

If you want to set `pitch` and `yaw` for a helix, put them into a `HelixGroup` with the desired pitch/yaw.

```python
# old
helix = sc.Helix(pitch=90)


# new
pitch90_helix_group = sc.HelixGroup(pitch=90)
helix = sc.Helix(group='pitch90')
```

Note that the name of a `HelixGroup`, such as `'pitch90'` above, is not stored in the `HelixGroup` itself, but is the key in the map [Design.groups](https://scadnano-python-package.readthedocs.io/en/latest/#scadnano.Design.groups).

If you want to look up `pitch` and `yaw` of a helix, use `design.pitch_of_helix` or `design.yaw_of_helix` respectively:

```python
pitch = design.pitch_of_helix(helix)
yaw = design.yaw_of_helix(helix)
```

which returns the `pitch` and `yaw` of `helix`'s `HelixGroup`.

`Helix.roll` is still a valid field. Use `design.roll_of_helix` to compute the roll of the helix *added to* the roll specified in its `HelixGroup`:

```python
roll = design.roll_of_helix(helix)
```


### `grid` fields are of type `Grid`, not `str` or `None`

Library now checks that `grid` fields are of type `Grid`, not `str` or `None`. This will break existing code that was manually passing in strings for the `grid` parameter.

## Commits
- [[a1ee18b](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/a1ee18ba47739f886fbd1fdedaf4066c6eb9aa16)]: bumped version to 0.15.2 (David Doty) [#175](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/175)
- [[60a07da](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/60a07da0f9fcbe01da918bf0251e9fce4babe3cd)]: remove debug print hello (Constantine Evans) [#165](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/165)
- [[e322fbe](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/e322fbe701b23735ea9d4b5b032cced19bbf12fd)]: Fixes [#163](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/issues/163); Remove Helix.pitch and Helix.yaw ([#166](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/issues/166)) (Benjamin Lee) [#166](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/166)
- [[bccc3ad](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/bccc3ad4cd03ee75856f0fb7b9934882c9464af3)]: added docstring to Strand.idt_dna_sequence and stub for to_oxdna_format (David Doty) [#175](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/175)
- [[860ae9a](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/860ae9ad407a65750571a99e84bb75a533be99b9)]: Encoder indent suppression: str.replace → re.sub (fixes [#170](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/issues/170)) (Constantine Evans) [#171](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/171)
- [[93a8135](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/93a8135fd112c0d3160ea2b073a18280e9d69162)]: add optional parameter suppress_indent to method Design.write_scadnano_file; closes [#170](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/issues/170) (David Doty) [#175](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/175)
- [[1f0a233](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/1f0a2335ee8ec23a33543b59c7cff8e10129be2b)]: Fix Sphinx warning on idt_dna_sequence (Benjamin Lee) [#172](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/172)
- [[8ca07f9](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/8ca07f9245ba02fbd59f4753efff896392595ef8)]: fixed indenting on idt_dna_sequence docstring to suppress Sphinx warning (David Doty) [#175](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/175)
- [[da1eadd](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/da1eadd8697c38350a7eeb82a076f91458072e8e)]: fixed write_scadnano_file and to_json docstrings (David Doty) [#175](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/175)
- [[418d4b5](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/418d4b5742ceac7c044efa48a0c98905d81cbb1b)]: Merge branch 'dev' into dev-fix-sphinx-warning (David Doty) [#172](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/172)
- [[f0fbe53](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/f0fbe53179f7783dd7fca849c70099991cfaaee9)]: updated Sphinx to add view code option, and put type hints in parameter and return description (David Doty) [#175](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/175)
- [[cad6674](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/cad667416f5a1e369d126b0124bd3dbf670c4d7d)]: fixed some docstrings and autoformatted code (David Doty) [#175](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/175)
- [[73bb35d](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/73bb35d019dd26a9146f8bd4e0251b61c71c7976)]: fixed PEP warnings (David Doty) [#175](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/175)
- [[35b33aa](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/35b33aaaddcb173a2ede8f659f080b81aa8afa0b)]: fixes [#168](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/issues/168): ensure grid fields are of type Grid, not str or None (David Doty) [#175](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/175)
- [[7f1f591](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/7f1f5915cc3e5515ba22e2b13c227f518174cadb)]: Update tutorial.md (David Doty) [#175](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/175)
- [[737570e](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/737570e8d63c2328b23f3d900be023d8570e19cb)]: Redirect workflow from master to main (Benjamin Lee) [#175](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/175)
- [[08ecf96](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/08ecf96562df7b139b13cb8290e9640f264aa2ae)]: Fixes [#176](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/issues/176); Rewrite doc test (Benjamin Lee) [#178](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/178)

[Changes][v0.16.0]


<a id="v0.15.1"></a>
## [v0.15.1](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/releases/tag/v0.15.1) - 2021-01-15

Not much new, mainly this release is to alter the documentation to clarify that scadnano is a separate project from cadnano.

## Commits
- [[fb249d8](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/fb249d8181dff1b9fe81b6d52953f39a7540e2d5)]: bumped version to 0.15.1 (David Doty) [#162](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/162)
- [[6ce55c2](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/6ce55c277c956efbf280ba770f67ff7670afaef1)]: added to README note clarifying difference between scadnano and cadnano, and table of contents (David Doty) [#162](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/162)
- [[cb09bfa](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/cb09bfadc2352be0418ac57c925a6f6728b6008e)]: added to API documentation note clarifying difference between scadnano and cadnano (David Doty) [#162](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/162)
- [[bed7b3c](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/bed7b3c40807d786e6bc7c2c7fad3f9ebee85e7a)]: Update tutorial.md (David Doty) [#162](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/162)
- [[05e76b9](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/05e76b95c845f59a2ab23960ec525002ced5b354)]: updated API docs to link from assign_m13_to_scaffold method to M13Variant enumerated type (David Doty) [#162](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/162)

[Changes][v0.15.1]


<a id="v0.15.0"></a>
## [v0.15.0](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/releases/tag/v0.15.0) - 2021-01-14

**BREAKING CHANGE:** The field `Strand.idt.name` has been removed. Since each `Strand` now has a `name` field, that is used instead. So replace this:

```python
# old
idt = IDTFields(name='staple1', plate='plate1', well='A1')
strand = Strand(domains=[domain], idt=idt)
```

with this:

```python
# new
idt = sc.IDTFields(plate='plate1', well='A1')
strand = sc.Strand(domains=[domain], idt=idt, name='staple1')
```

and replace this:

```python
# old
design.strand(0,0).move(8).with_idt(name='staple1', scale='25nm')
```

with this:

```python
# new
design.strand(0,0).move(8).with_name('staple1').with_idt(scale='25nm')
```

## Commits
- [[449f516](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/449f516ddceb18d5530def85cbe0565fdf653df6)]: updated docstrings (David Doty) [#161](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/161)
- [[658e9e4](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/658e9e4df129b0a2b0dbefe2d9997e05cf83703b)]: closes [#154](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/issues/154); re-make tutorial; The tutorial is now up-to-date with scadnano features as of Dec. 2020. (David Doty) [#161](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/161)
- [[30cbfcc](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/30cbfcc90a5e53d9091feeeb09f6d578b0b9fdba)]: Update tutorial.md (David Doty) [#161](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/161)
- [[6a3dea8](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/6a3dea8d4139b1efbec1383718e7f319661950de)]: Update tutorial.md (David Doty) [#161](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/161)
- [[4a122e9](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/4a122e9d3849986f4a469dcdc473076700c3bfb4)]: added note about offsets being inclusive for half crossovers (David Doty) [#161](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/161)
- [[1de2aa8](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/1de2aa876594204eb091b074f3a2f6043ffab63c)]: Update tutorial.md (David Doty) [#161](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/161)
- [[2c2cb1e](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/2c2cb1eb863a4789e453769baf3822d195da4604)]: Update tutorial.md (David Doty) [#161](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/161)
- [[a5969f0](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/a5969f0a8f0e7782e5e44a81197a0d53e92a069f)]: Update tutorial.md (David Doty) [#161](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/161)
- [[7c03a5d](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/7c03a5df9a6f8d3e3820c19481ada60140ea4887)]: Update tutorial.md (David Doty) [#161](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/161)
- [[1c09d6f](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/1c09d6f5a0b0cfedb29162aa3ef8d5753353961f)]: Update tutorial.md (David Doty) [#161](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/161)
- [[c4ea9c0](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/c4ea9c0555a170024e9f8a27f305b3b4fd7dbcc9)]: added note explaining type hints in tutorial (David Doty) [#161](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/161)
- [[77f3a05](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/77f3a059c08f6871299d44b48296176b67b254a9)]: Update tutorial.md (David Doty) [#161](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/161)
- [[0fe9bf6](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/0fe9bf6881ea1d672b73fd7503253a8c379bf504)]: added screenshot of Excel file to tutorial (David Doty) [#161](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/161)
- [[a563263](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/a563263a9c2c457d32e41bb62dddd93b3a3f70c6)]: Update tutorial.md (David Doty) [#161](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/161)
- [[3254605](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/3254605c3b7c7e9bbfb972dcce0a52466fa473c7)]: closes [#159](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/issues/159); remove IDT name field; now Strand.name is the official name of the strand; bumped minor version since this is a breaking change (David Doty) [#160](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/160)
- [[b096e5d](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/b096e5da5ab88322ee802565c9d52e080d19be7f)]: updated installation instructions to indicate how to test if scadnano was installed successfully (David Doty) [#161](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/161)
- [[cfb2ffd](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/cfb2ffd2b14d2becba3f58754d95fe5307a6d93b)]: Update README.md (David Doty) [#161](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/161)

[Changes][v0.15.0]


<a id="v0.14.0"></a>
## [v0.14.0](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/releases/tag/v0.14.0) - 2020-12-25

# Exporting IDT files with default values for idt fields no longer requires first modifying the Strand.idt field
It is no longer necessary to modify the design to set the field `Strand.idt` in each Strand before calling the methods that export DNA sequences in IDT-formatted files. For staple strands with no idt field, a reasonable default for each value will be chosen.

So it is now possible to do this:

```python
import scadnano as sc

design = sc.Design(helices=[sc.Helix(max_offset=100) for _ in range(2)], grid=sc.square)
design.strand(0, 0).move(8).cross(1).move(-8)
design.strand(0, 16).move(-8).cross(1).move(8)
design.assign_dna(strands[0], 'A'*16)
design.assign_dna(strands[1], 'C'*16)

# before the change, the next line would have skipped writing the two strands since they have no idt field set,
# now, reasonable defaults are used, without requiring the side-effect of writing the field Strand.idt
design.write_idt_plate_excel_file()

# to skip exporting strands that lack an idt field, specify the parameter only_strands_with_idt
# below, only the newly added strand with T's will be exported; the previous two will not
design.strand(0, 24).move(8).cross(1).move(-8).with_idt('only_strand_to_export')
design.assign_dna(strands[2], 'T'*16)
design.write_idt_plate_excel_file(only_strands_with_idt=True)
```

This implies several changes in the API

- **BREAKING CHANGE:** Changed the export methods so that, by default (with no parameters specified), they behave differently. In particular, now by default they will export DNA sequences for *all staple strands* (i.e., non-scaffold), using the `idt` field of the Strand if it is present, and otherwise using reasonable defaults, the same defaults that previously were stored in the Strand by calling `Strand.set_default_idt()`.

- **BREAKING CHANGE:** Removed the following:
  - field `Strand.use_default_idt`
  - method `Strand.set_default_idt()` 
  - method `Design.set_default_idt()`
  - parameter `use_idt_defaults` in function `origami_rectangle.create()`
  Now, if you want to set a Strand to have an `idt` field, it must be explicit, although the `IDTFields` constructor only requires a `name` parameter, so it's as easy as `strand.idt = IDTFields('name_of_strand')` if you are happy with the defaults for other `idt` fields such as `idt.purification`.

- **BREAKING CHANGE:** Removed parameter `warning_on_non_idt_strands` from the IDT export methods on `Design`. Now, you can either ask those methods to skip exporting Strands lacking an `idt` field by setting the parameter `only_strands_with_idt` to True, or let all (non-scaffold) strands be exported by setting `only_strands_with_idt` to False (the default).

- Added parameter `export_scaffold` to DNA sequence export methods to allow the scaffold(s) to be exported (False by default).


# `Crossover` class and bulk `Design.add_crossovers` method removed

- **BREAKING CHANGE:** (This one is unrelated to exporting IDT files; it is related to the circular strands implemented in [v0.13.4](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/releases/tag/v0.13.4).) Since circular strands make it easier to use the `Design.add_half_crossover` and `Design.add_full_crossover` methods, we have removed the method `Design.add_crossovers` and the type `Crossover`. Previously, that method helped avoid creating circular strands by allowing one to build up a list of Crossovers and add them in bulk, where adding them one at a time would have resulted in an intermediate circular strand, even if the final result had all linear strands. Now that circular strands are supported, this is no longer needed. The recommended method of adding many crossovers at once is simply to call `Design.add_half_crossover` and/or `Design.add_full_crossover` repeatedly, i.e., replace

  ```python
  crossovers = [
      Crossover(helix=0, helix2=1, offset=16, forward=True, half=True),
      Crossover(helix=0, helix2=1, offset=24, forward=False, half=True),
      Crossover(helix=2, helix2=3, offset=32, forward=True),
      Crossover(helix=2, helix2=3, offset=40, forward=False),
  ]
  design.add_crossovers(crossovers)
  ```

  with this instead:

  ```python
  design.add_half_crossover(helix=0, helix2=1, offset=16, forward=True)
  design.add_half_crossover(helix=0, helix2=1, offset=24, forward=False)
  design.add_full_crossover(helix=2, helix2=3, offset=32, forward=True)
  design.add_full_crossover(helix=2, helix2=3, offset=40, forward=False)
  ```



## Commits
- [[d2e924c](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/d2e924cb79a75bea75e05d9e4fa226ea460e3eb8)]: closes [#111](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/issues/111): write_idt_plate_excel_file uses reasonable defaults even when some strands have no IDT field set; bumped version (David Doty) [#157](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/157)
- [[b52a046](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/b52a0463164aefcd7e1a9eb6e9f796a91e18bbf4)]: updated examples to fully type annotate all functions and avoid name shadowing in if __name__=="__main__" blocks (David Doty) [#157](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/157)
- [[87a47d5](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/87a47d5da066a45dbc915b304bf92a494257ed24)]: fixed bug where scaffold property being lost when joining two strands (at least one of which was scaffold)by crossover; also reworked write_idt_plate_excel_file to work properly with default idt name if idt fields are not present in some strands (David Doty) [#157](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/157)
- [[b17a6a1](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/b17a6a1a3555001844ff13c87d4c51716ace76f8)]: minor documentation and identifier name updates (David Doty) [#157](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/157)
- [[c6620df](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/c6620dfbb87097892bac48a447e97cd25b2218e9)]: cleaned up old links in package docstring (David Doty) [#157](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/157)
- [[e23f361](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/e23f361d95ab14e20f92c09248b1c5ce41e40026)]: re-ran examples after last commit (David Doty) [#157](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/157)
- [[0667935](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/0667935c6e1fbe7d7b6436aa2cb6324dc7918ba1)]: inlined creation of empty design in tutorial script (David Doty) [#157](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/157)
- [[424f478](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/424f478cc49744c2397a62a4255aef6248e45c4d)]: Update scadnano.py (David Doty) [#157](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/157)
- [[dd63418](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/dd63418bd08a54331bde68c5fee444c2989d1ec7)]: made Loopout generic type parameterized by DomainLabel (David Doty) [#157](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/157)
- [[7574058](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/7574058a5cdb7b717a5da8d80b53f1295742e046)]: Update scadnano.py (David Doty) [#157](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/157)
- [[9e32cd4](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/9e32cd404cce18c1ac63c33517e2da2dedbe1527)]: fixed documentation of Strand.rotate_domains (David Doty) [#157](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/157)
- [[3da17ec](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/3da17ec4d3508b4b0d074ae5d7d902cf92f4c3a9)]: Merge branch 'dev' into write_idt_plate_excel_file-uses-reasonable-defaults (David Doty) [#157](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/157)

[Changes][v0.14.0]


<a id="v0.13.4"></a>
## [v0.13.4](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/releases/tag/v0.13.4) - 2020-12-24

# Circular strands
Circular strands are now supported. (See [#14](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/issues/14).) In general, it is not recommended that a final design have circular strands. In particular, there are aspects of scadnano, such as naming conventions for strands and conventions for assigning DNA sequences, that assume the strand has a 5' and 3' end. Under the hood, the domains of a circular strand are still listed in some order, with the same constraint as before that a strand cannot begin or end with a loopout; see [#2](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/issues/2). In particular, this means a circular strand must have at least one crossover; it cannot be all loopouts linking the domains.

However, circular strands are convenient for the intermediate steps of a design, allowing one to add crossovers and join strands by ligation without worrying whether it will create a circular strand. But it is recommended, particularly before assigning DNA sequences, to linearize all circular strands (i.e., add a nick to break some domain into two, or remove a crossover somewhere). This includes even strands such as those representing M13 that are naturally circular. Otherwise the effect of assigning a DNA sequence is undefined. Operations that circularize and linearize strands with DNA sequences already assigned are similarly undefined and may change the DNA sequence in unexpected ways.

# Strand DNA sequence export order
Previously, when exporting DNA sequences, they would be exported in whatever order the strands appear in the Design. Now, there are a few reasonable options to choose from. See [#147](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/issues/147).

## Commits
- [[f524a56](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/f524a56e1b78e9db83a42510edb78fa19853693c)]: closes [#147](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/issues/147): allow row-major and column/major order of export of strands (David Doty) [#156](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/156)
- [[481d966](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/481d9669211506cb4f966d6995689cbc2b63f35a)]: updated unit tests (David Doty) [#156](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/156)
- [[76b4b64](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/76b4b64fb9ad4c345501c732c7cd14247671882c)]: fixed documentation of hex coordinate system (David Doty) [#156](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/156)
- [[a1e7467](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/a1e74674b7bdcc61208a06ac4e420938cbf9868d)]: fixed documentation of hex coordinate system (David Doty) [#156](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/156)
- [[9d1ac96](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/9d1ac962405a73bb314fcb9d3ee161d04a0c0fce)]: fixed documentation of Loopout to say they are not allowed on ends of Strand, and added example of creating one with chained method calls. (David Doty) [#156](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/156)
- [[4b13ef2](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/4b13ef26ccf1fee51e0f63e4571b7861fe608f3e)]: closes [#14](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/issues/14); add support for circular Strands (David Doty) [#152](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/152)
- [[893075f](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/893075f37846b2c153a725f975afb3c76da68765)]: fixed docstring errors and bumped version (David Doty) [#152](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/152)
- [[7f5d3c9](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/7f5d3c96556721d825dd4ae8f94042369a627214)]: updated add_half_crossover and add_nick to handle circular strands (David Doty) [#153](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/153)
- [[84db6a8](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/84db6a8bc93cbdde38daa02d2e640193c721d937)]: added ligate method to design and associated unit tests (David Doty) [#156](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/156)
- [[4a7dd1d](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/4a7dd1d6235519ed26044ef7585ec34eb881672d)]: added unit tests for nicking on circular strand with 3 domains and a loopout, nicking on all three domains (David Doty) [#156](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/156)
- [[8eda390](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/8eda390b180f2f27ba924f74d17ef07dab5cf96e)]: Update scadnano_tests.py (David Doty) [#156](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/156)
- [[c58b832](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/c58b8320da8a0607a198c1980e41d846def165c9)]: modified circular strand example script to have length 8 domains (David Doty) [#156](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/156)
- [[5b1d084](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/5b1d08421fc6aa1c66f209dd393012e023dff9cc)]: bumped version to match web interface (David Doty) [#156](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/156)
- [[22d446e](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/22d446eced1fccae06974ce84342af9421842e9c)]: updated tutorial script (still not completely in line with [#154](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/issues/154), however, since issue [#111](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/issues/111) still need to be implemented) (David Doty) [#156](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/156)
- [[5f91962](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/5f91962e9aab8e31281b4fea5edc38200d1ab9cf)]: Modifications not needed in circular example (Cosmo) [#155](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/155)
- [[4ee5e28](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/4ee5e288d4c94fb57791845b2c5df79f7e7470ed)]: Importing circular strands from cadnano to scadnano is now functional. Two tests added:   which come from cadnano Autostaple's output (Cosmo) [#155](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/155)
- [[6e27f5f](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/6e27f5f97ac3defafe3a754c127f7f51211ffdfa)]: Exporting designs with circular strands is working, un ignoring tests_inputs/cadnano_v2_export as some .sc files are not generated but directly read by the tests (Cosmo) [#155](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/155)
- [[675b542](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/675b5428385bb8885cb25ba7ad45a225223608a7)]: fixed mypy and PEP errors in cadnano import/export code (David Doty) [#155](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/155)

[Changes][v0.13.4]


<a id="v0.13.0"></a>
## [v0.13.0](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/releases/tag/v0.13.0) - 2020-11-22

## Breaking Change: Swap x and z coordinates

Python scripts using the "none" grid and specifying 3D coordinates will need to be rewritten so that x and z coordinates are swapped.

Previously, positive x moved right in the main view and into the screen in the side view, and positive z moved right in the side view and out of the screen in the main view. Now these are swapped.

If you are only using the web interface, you shouldn't see any change. But Python scripts that specify (*x*,*y*,*z*) coordinates will need to be updated to swap the roles of *x* and *z*.

If you are curious why this was done, read here: [UC-Davis-molecular-computing/scadnano#488](https://github.com/UC-Davis-molecular-computing/scadnano/issues/488)

## Commits
- [[c6385f0](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/c6385f0b49d982a3edf321d6c627e170beac0b2e)]: Fixes [#138](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/issues/138); add Python 3.6 CI tests (Benjamin Lee) [#140](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/140)
- [[0546749](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/0546749c85001edc77ed26f22a69851a1af0ff3c)]: Fixes [#125](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/issues/125); CI checks for docs and PyPI packaging (Benjamin Lee) [#141](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/141)
- [[76a85be](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/76a85be1c4922099a12a1f9fcfd0c3d80b0c92d2)]: Remove misplaced name field (Benjamin Lee) [#141](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/141)
- [[8699cb7](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/8699cb7d9350b5f87e308058d45820c6305b6898)]: Remove publish from task name for clarity (Benjamin Lee) [#141](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/141)
- [[b137c80](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/b137c8080cc663ff82406fe792998383a869500e)]: added FAM, ROX, and Fluorescein modifications; added code to automatically population field Modification.id with Modification.idt_text if latter is specified and former is not (David Doty) [#146](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/146)
- [[7486223](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/74862239bbf6fce21e5935b6c2722eb6fe47da3c)]: fixed error in setting Modification.id from Modification.idt_text (David Doty) [#146](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/146)
- [[b6f75d4](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/b6f75d488d93d94abb14ecadb61510d735ab34a1)]: included everything in modifications and origami_rectangle when importing scadnano (David Doty) [#146](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/146)
- [[61ba537](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/61ba5377bb9bed178431553a776749e5469e56cb)]: updated defaults for Design.set_default_idt (David Doty) [#146](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/146)
- [[11d6606](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/11d6606b9131a2de46273373dc945f361539887b)]: fixed some mypy errors, and added unique_names argument to set_default_idt to break with cadnano's naming convention and ensure strand names are unique (David Doty) [#146](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/146)
- [[11d9d67](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/11d9d675028d82c2dd8e0c7be88c909653547be6)]: changed default IDT purification with modifications from PAGE to HPLC (David Doty) [#146](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/146)
- [[c24784f](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/c24784f79ac53fab16344ba24521d13325bc24c5)]: added code to import IDT fields from JSON, along with unit tests (David Doty) [#146](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/146)
- [[ecf5ced](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/ecf5ced2ebd2d194f07630d6145c4583c5148596)]: Adds PR to CI workflow (Benjamin Lee) [#143](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/143)
- [[0bd0082](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/0bd0082f824c7088f5b69a81a48b5482ec594d35)]: Remove push events from CI workflow (Benjamin Lee) [#143](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/143)
- [[8787c6d](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/8787c6da83a01aed7897b7abbb18a26bbfadc0d2)]: Update README.md (David Doty) [#146](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/146)
- [[0e5b323](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/0e5b323c3788a981961b63e68ce3ca45dbdbbca7)]: corrected relative link to .sc file in tutorial (David Doty) [#146](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/146)
- [[2ad4868](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/2ad4868e3f26ed3c8cf0967fb50cac6d6703186a)]: Fixes [#144](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/issues/144); swap position x z coordinate interpretation (Benjamin Lee) [#145](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/145)
- [[34d24fe](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/34d24fe4116fd1dcda04e610c89b7738e0b825b5)]: Rewrite none-grid example scripts for [#144](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/issues/144) (Benjamin Lee) [#145](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/145)
- [[a2a31f6](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/a2a31f60baf5001936cd8f613cb6a7bfe2350fa9)]: Update proposal example for [#144](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/issues/144) (Benjamin Lee) [#145](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/145)
- [[3382247](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/33822476167006e60bcc896abf328a60e68b017c)]: updated docstrings for Position3D (David Doty) [#146](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/146)
- [[10eca78](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/10eca786ccfc87611230b047c97f0af0ab571a0a)]: updated StrandBuilder docstrings (David Doty) [#146](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/146)
- [[1c36eeb](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/1c36eeb3595e74bf982a57fd07b1c41511336985)]: Update scadnano.py (David Doty) [#146](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/146)
- [[b91ef99](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/b91ef990e6296f9b439a954423f316317b54a8b7)]: reverted previous edit of docstrings (David Doty) [#146](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/146)
- [[a94102c](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/a94102cf0d72548f3d09163b606d451d82f596ff)]: fixed erroneous reference to method to() to be move() (David Doty) [#146](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/146)

[Changes][v0.13.0]


<a id="v0.12.2"></a>
## [v0.12.2](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/releases/tag/v0.12.2) - 2020-09-11

For Python 3.6, dataclasses backport library is now automatically installed when installing via pip.

The README still contains instructions for installing it manually, in case a user wants to simply download the scadnano.py file to use with Python 3.6. But they should be updated to emphasize that it is unnecessary to install it manually when installing scadnano via pip.

## Commits
- [[1102804](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/11028041680e1f8e5595fed67df9f8ba31e95d9e)]: Adding installation instruction for python dataclasses for Python 3.6 (Cosmo) [#135](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/135)
- [[1999445](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/1999445cc1e09db2e15fb43a3b7763ef242109db)]: Require python 3.6 min (Cosmo) [#135](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/135)
- [[2061a9f](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/2061a9fb7e67e99e8c98f5e601a65d4744111375)]: added instructions about commit messages to CONTRIBUTING; bumped version (David Doty) [#137](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/137)

[Changes][v0.12.2]


<a id="v0.12.1"></a>
## [v0.12.1](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/releases/tag/v0.12.1) - 2020-09-10

### cadnano export upgraded
Can now export to cadnano from design that uses multiple helix groups.

### rotate domains
Added method `Strand.rotate_domains` allowing domains to be "rotated" on a strand. Think of it like adding a crossover between the 5' and 3' ends, and removing another crossover.

### type hints
Beefed up type hints, in particular making type variables `StrandLabel` and `DomainLabel` for the types of strand labels and domain labels, allowing `Design`, `Strand`, and `Domain` to be parameterized by them, so that a static type checker such as [mypy](http://mypy-lang.org/) will catch errors such as

```python
domain: Domain[str] = Domain(label=123)    # error
domain2: Domain[str] = Domain(label='123') # fine
domain2.label = 123                        # error
```

NOTE: this particular type hint no longer makes sense since labels are now assumed to be strings: https://github.com/UC-Davis-molecular-computing/scadnano-python-package/releases/tag/v0.18.0

Also got rid of all previous [mypy](http://mypy-lang.org/) type errors.

## Commits
- [[41a7f70](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/41a7f701aabf56fba4e48dc04c2859b3e9cf470d)]: updated paper URL now that DNA 2020 paper is published (David Doty) [#136](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/136)
- [[304e120](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/304e120c1a0de989ce375cacbe3ab83c3c918cd5)]: added rotate_domains method to Strand to "rotate" domains of strand (i.e., like adding a crossover between the 5' and 3' ends, and removing another crossover) (David Doty) [#136](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/136)
- [[b7e41d0](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/b7e41d0de73d5b94cbf398f8e2ec0a9d6efcef89)]: Design.write_scadnano_file now warns if a Loopout is the first or last substrand on a Strand (still allowed in intermediate designs) (David Doty) [#136](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/136)
- [[903510d](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/903510d7f27b19b6275ef6c1d9475cc7e72b9ba7)]: fixed all mypy warnings; closes [#109](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/issues/109) (David Doty) [#136](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/136)
- [[b2012e9](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/b2012e9964437532cc5ee81bd9e0707961f3c50f)]: fixed mypy warnings in origami_rectangle (David Doty) [#136](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/136)
- [[a3248bc](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/a3248bca145ddde9725cc0a16b191ed692248c4a)]: made Design, Strand and Domain generic parameterized by StrandLabel and DomainLabel. Made strand label not indented in serialized JSON (David Doty) [#136](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/136)
- [[080bf6d](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/080bf6de7f005cbdde6b6613953a536f03b9e5a4)]: added example with domains names (some mismatching) (David Doty) [#136](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/136)
- [[a2a6b07](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/a2a6b072a21f83578c9aea5802eb25cc7e5eb829)]: updated names example to have more kinds of mismatches (David Doty) [#136](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/136)
- [[27baee2](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/27baee2d4e27ee7d8d5995418105acb6b7610ea9)]: annotated variable to quiet mypy (David Doty) [#136](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/136)
- [[c5bdb1c](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/c5bdb1c39744612df82f16d8a5cd3b32509a605b)]: Update names_domains_strands.py (David Doty) [#136](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/136)
- [[9b54b1e](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/9b54b1e5fc8d18960aef068fae487862a5d8eac9)]: Export code supports helix groups and associated unittest. test_6_helix_bundle_honeycomb restored. (Cosmo) [#133](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/133)
- [[4d62212](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/4d622123c3bbc11bbb99c2c826a9d3504fd9ac49)]: Export code supports helix groups and associated unittest. test_6_helix_bundle_honeycomb restored. (Cosmo) [#133](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/133)
- [[399f905](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/399f9055f16a80e01f22bcdd1fb767c07711ff73)]: Export code supports helix groups and associated unittest. test_6_helix_bundle_honeycomb restored. (Cosmo) [#133](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/133)
- [[6233075](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/62330752fbe511c95c6a1d1e8688f03467c99639)]: Export code supports helix groups and associated unittest. test_6_helix_bundle_honeycomb restored. (Cosmo) [#133](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/133)
- [[a274a98](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/a274a98999c94f7c5bf11b9c64a1b4af1ef50b0c)]: Forcing the add of test_6_helix_bundle_honeycomb.sc which was ignored (Cosmo) [#133](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/133)
- [[3753e36](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/3753e36022364d91fc9f4d6b55e739ee00f4234c)]: Correct type annotation (Cosmo) [#133](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/133)
- [[627d45b](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/627d45b993c9ebbbecd9aeaa493dafbadf2102ca)]: Correcting syntax error (Cosmo) [#133](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/133)
- [[db82d47](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/db82d47fd3fdb7f208afe78d978a30a11217a083)]: bumped version (David Doty) [#136](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/136)

[Changes][v0.12.1]


<a id="v0.12.0"></a>
## [v0.12.0](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/releases/tag/v0.12.0) - 2020-09-01

Optional field `name` now supported in Strand, Domain, and Loopout.

The domain/loopout name is displayed (optionally) in the scadnano web interface main view. The strand name is displayed on mouseover, in the tooltip that pops up, and if "backbone mode" is selected, in the footer at the bottom of the page.

All three are also used with the dsd DNA sequence designer (not public yet).

## Commits
- [[3ebcb86](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/3ebcb86aaabaa970ee8610ad86b79619cd340f20)]: added move (relative offset) chained method to README (David Doty) [#131](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/131)
- [[eb7563f](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/eb7563f549e4cadc4db58a4a2cea8884a55f89e7)]: added link to json.dumps documentation in docstrings for Loopout.label and Domain.label (David Doty) [#131](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/131)
- [[ca482ba](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/ca482ba8708cb8a775dd2ac5bb9058a906cd0daa)]: added parameter `check_length` to Design.assign_dna that enforces the sequence is exactly the length of the Strand/Domain being assigned (David Doty) [#131](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/131)
- [[2afc921](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/2afc921db30bdedac712d1b4932c78dab420b3a4)]: added optional name fields to Strand, Domain, and Loopout; these are used now instead of label to assign names to Strands, Domains, and Loopouts in the dsd DNA strand designer (though Strand labels are still used to assign Strand groups in the dsd sequence designer) (David Doty) [#131](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/131)
- [[db5d31f](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/db5d31ff46824607c95667390a2b85a45ce53a43)]: bumped version to 0.12.0 (David Doty) [#131](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/131)
- [[5bec304](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/5bec30477e0a5b42e452294eb8d87df459b68acd)]: removed _version.py file (David Doty) [#131](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/131)

[Changes][v0.12.0]


<a id="v0.11.2"></a>
## [v0.11.2](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/releases/tag/v0.11.2) - 2020-08-25

Added loopout labels. Now you can do this:

```python
(design.strand(0, 0).to(8)
 .loopout(1, 4).with_domain_label('loopout label')
 .to(0))
```



## Commits
- [[ed0ea3b](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/ed0ea3bebf7ea1ae3d6af62e0872977c69c0fc5d)]: removed example; made StrandBuilder._strand private and added getter that raiases exception if Strand has not been created yet (David Doty) [#130](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/130)
- [[23d5a1c](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/23d5a1c3e302f2e81e62461b3e3493d3a8e72585)]: fixed unit test after changing indenting JSON behavior if Helix.position is specified (David Doty) [#130](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/130)
- [[9a969ce](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/9a969ce3a538f342339d824d5198c5c0ddd06a58)]: added loopout labels; bumped version (David Doty) [#130](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/130)

[Changes][v0.11.2]


<a id="v0.11.1"></a>
## [v0.11.1](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/releases/tag/v0.11.1) - 2020-08-25

The main new feature is a method `StrandBuilder.move` and parameter `move` in the existing methods `StrandBuilder.cross` and `StrandBuilder.loopout`. `move` is like `to`, but is a relative offset rather than an absolute offset.

For example, to make a strand that, starting at helix 0, offset 123, goes 8 bases forward, then crosses over to the same offset on helix 1 and goes 16 bases in reverse, then crosses over to helix 2 (but the crossover jumps 2 offsets back), then goes 10 bases forward:

```python
(design.strand(0, 123)
 .move(8)
 .cross(1)
 .move(-16)
 .cross(2, move=-2)
 .move(10))
```

![image](https://user-images.githubusercontent.com/19274365/91471756-55ece080-e84b-11ea-82c3-c345c42e1485.png)

where the parameter to the method `move`, as well as the parameter named `move`, are relative to the current offset. This is equivalent to the more cumbersome absolute offsets with `to`:

```python
(design.strand(0, 123)
 .to(131)
 .cross(1)
 .to(115)
 .cross(2, offset=113)
 .to(123))
```

## Commits
- [[37db033](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/37db033c9827c1ea9390210c8bc26ec4dc0976cd)]: updated tutorial with correct parameter names for Crossover constructor (David Doty) [#129](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/129)
- [[0808d88](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/0808d88b6b298286224cd2b63350cb4479cfb436)]: fixed code display in docstrings (David Doty) [#129](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/129)
- [[2499733](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/2499733a64fe328befe948e1fd3577ff2f2b9851)]: added move method to StrandBuilder and move parameter to StrandBuilder.cross and StrandBuilder.loopout to enable relative (instead of absolute) specification of offsets when creating Domains through chained methodso (David Doty) [#129](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/129)
- [[0b43e14](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/0b43e14f36464d5672857f5ace455d9886ffa7ac)]: added examples; bumped version (David Doty) [#129](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/129)
- [[4760354](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/4760354cb189915fd379f33bcf5fd2e0d15d8b3f)]: fixed unit test now that positions are automatically created in Helix when grid is None (David Doty) [#129](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/129)

[Changes][v0.11.1]


<a id="v0.11.0"></a>
## [v0.11.0](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/releases/tag/v0.11.0) - 2020-08-06

Note there are two breaking changes:

-  removed helix_template and num_helices from Design constructor. This feature wasn't carrying its weight, given the complexity it introduced to parse optional parameters in the `Design` constructor. It's simple enough to simply type something like 
    ```python
    helices=[Helix(idx=idx, <other properties you want>) for idx in range(num_helices)]
    ```

- removed `major_tick_distance` from `Design`. Similarly, this feature wasn't carrying its weight and REALLY complicated parsing optional parameters. To assign the same `major_tick_distance` to every `Helix`, simply give the same value in every `Helix` constructor, e.g., 
  ```python
  helices=[Helix(idx=idx, major_tick_distance=10) for idx in range(num_helices)]
  ```

## Commits
- [[b9e7eb6](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/b9e7eb6f9ea9386cba65001b7b56b68c6afab9c3)]: ensures default helices view order is assigned properly in each helix group (David Doty) [#128](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/128)
- [[6874906](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/6874906a7f2fa0bb1714b66f929a86a7831ff919)]: BREAKING CHANGE: removed helix_template and num_helices from Design constructor (David Doty) [#128](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/128)
- [[419caeb](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/419caeb637de5587157915b2af73576b20ae504d)]: bumped version for breaking change (removed helix_template and num_helices from Design constructor parameters) (David Doty) [#128](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/128)
- [[e7e223f](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/e7e223fc408052557119be6d3371ad7377541ba3)]: ignoring .sc files in tests_inputs/cadnano_v2_export directory (David Doty) [#128](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/128)
- [[fb78aa1](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/fb78aa16a1e476bf9a026fa6d2348bd76b4122ad)]: BREAKING CHANGE: removed major_tick_distance from Design (David Doty) [#128](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/128)
- [[d36d457](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/d36d457d53b00c0278e97e1be14209a1ac53742f)]: added example with helix groups (David Doty) [#128](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/128)
- [[507c601](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/507c60166fb254696229b23acf403e631bccd7a2)]: closes [#104](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/issues/104); add support for Helix.major_tick_start and Helix.major_tick_periodic_distances (David Doty) [#128](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/128)
- [[3f59e6d](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/3f59e6dc25e4fecbc71a9910d12d4fef9b816aa5)]: re-ran examples to produce latest version (David Doty) [#128](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/128)
- [[107b677](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/107b677381e0e98c60c275e92bb1af154f8d934a)]: removed DomainLabel and StrandLabel as types in typing hints (David Doty) [#128](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/128)
- [[dceefc1](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/dceefc169b6c25d5788162342d00899795ff24b1)]: fixes [#126](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/issues/126); fix bug where strands are not always assigned a color (David Doty) [#128](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/128)

[Changes][v0.11.0]


<a id="v0.10.3"></a>
## [v0.10.3](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/releases/tag/v0.10.3) - 2020-07-24

This was mainly to test some functionality with auto-generating docs and PyPI releases.

## Commits
- [[42e6fbd](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/42e6fbdfe5e7cffbcdcbab6a60d3bad125638801)]: fixed version-finding function in setup.py (David Doty) [#124](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/124)
- [[2e81a57](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/2e81a57cc78c3929d631608139ced75055a0e478)]: moved __version__ up closer to top of scadnano.py (David Doty) [#124](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/124)
- [[ab79777](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/ab7977746a34e8c83ae360bc333da585a37c6190)]: removed print statement from setup.py (David Doty) [#124](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/124)
- [[a2c4763](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/a2c476383e7ba403e00c37cc962c436bee26d220)]: Update conf.py (David Doty) [#124](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/124)

[Changes][v0.10.3]


<a id="v0.10.2"></a>
## [v0.10.2](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/releases/tag/v0.10.2) - 2020-07-24

Introduces Helix groups in preparation for implementing them in the web interface: [UC-Davis-molecular-computing/scadnano#249](https://github.com/UC-Davis-molecular-computing/scadnano/issues/249)

Helix groups allow groups of helices to be grouped and given their own position, orientation, and grid, to help with designs where not all helices are parallel.

## Commits
- [[acd4e81](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/acd4e819712e43a9909e03fe63ba7b214c340f1f)]: Update README.md (David Doty) [#123](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/123)
- [[f5bbd38](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/f5bbd38c5baa6a3dc474b292655e5eea0997ba56)]: updated version in test files (David Doty) [#123](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/123)
- [[05db22d](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/05db22d3a7fc782b98b5db54777f110ffc509c40)]: removed in_browser test from examples and re-ran with new version (David Doty) [#122](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/122)
- [[24bd68e](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/24bd68e89a0ef068bfc30fdea78f77a9ada08f23)]: closes [#121](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/issues/121); introduces Helix groups to allow groups of helices to be grouped and given their own position, orientation, and grid, to help with designs where not all helices are parallel (David Doty) [#122](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/122)
- [[299e9ea](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/299e9eababcf91f7e3ae4c0ec04f995e40828362)]: removed _version.py so it only needs to be specified in scadnano.py (David Doty) [#123](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/123)
- [[631c1f3](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/631c1f3a8bf3f73907fa65f99ef2f8c7a3ca0464)]: minor docstring changes (David Doty) [#123](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/123)

[Changes][v0.10.2]


<a id="v0.10.1a"></a>
## [v0.10.1a](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/releases/tag/v0.10.1a) - 2020-07-23

Calling this release 0.10.1a, but really it's just to correct for the fact that I messed up the version number 0.10.1 in the source code on the last PR.

## Commits
- [[4d82e70](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/4d82e702d815f9522d6fc6fddcb27f7d977a5a11)]: incorrectly bumped version last time; bumping it in _version.py now (David Doty) [#120](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/120)

[Changes][v0.10.1a]


<a id="v0.10.1"></a>
## [v0.10.1](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/releases/tag/v0.10.1) - 2020-07-23

Mostly documentation was updated.

## Commits
- [[018fdf8](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/018fdf8a4e641171abe668a0e77dd78420a07ddc)]: added CONTRIBUTING doc (David Doty) [#119](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/119)
- [[ac3de4e](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/ac3de4ef8241847bf8682ecdef04d2464a7aa26e)]: removed bad link to images of checks pending/passing (from web repo) (David Doty) [#119](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/119)
- [[3871de6](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/3871de61fb18b5255def4cf790c6bc76cb813672)]: made contributing a subheader (David Doty) [#119](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/119)
- [[b908561](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/b9085618e39b048d5ce85dee62394c871bf3e295)]: Update CONTRIBUTING.md (David Doty) [#119](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/119)
- [[a601ea5](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/a601ea5a03fd6383375788dd3bc1b3eba1076e05)]: added screenshot showing how to download raw Python file from GitHub (David Doty) [#119](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/119)
- [[b5849b9](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/b5849b95c73b439509cfc28e1e15b5da3188645e)]: changed .dna extension to .sc in README and tutorial (David Doty) [#119](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/119)
- [[fbe4c20](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/fbe4c20facd4c29d244bc28e15772350b38f280d)]: lined up tutorial example source code (David Doty) [#119](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/119)
- [[56a5564](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/56a556445ee3ac6ae025ea71634906cb5478d05c)]: replaced --> with &rarr; in tutorial (David Doty) [#119](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/119)
- [[66db277](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/66db2774a82eeaadd618d8675c701d2089333f49)]: in tutorial, changed mistaken references from .py files to .sc files (David Doty) [#119](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/119)
- [[e287a88](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/e287a887adb1a096d56a45ab247e6d29ff443b23)]: added link to API in tutorial (David Doty) [#119](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/119)
- [[edadadc](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/edadadc6b235c3bde79bd97a50920a6b2c1ba7e4)]: in tutorial explained why we use add_strand, and changed default rotation for assigning M13 to 5587 from 5588 (David Doty) [#119](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/119)
- [[8c9daf6](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/8c9daf677a75aaac3191c4dc186e3faec33db4bf)]: updated relative links in README and tutorial (David Doty) [#119](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/119)
- [[4b27c51](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/4b27c510263e615962dd300a395b0e8f3436d4e8)]: added explanation in tutorial for why you need to manually reload .sc file after making changes with scripting library locally (David Doty) [#119](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/119)
- [[6b78b81](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/6b78b816b442aadd7249572036c6b2c10f44f3d5)]: removed local variable scaffold in tutorial code (David Doty) [#119](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/119)
- [[dfe38fa](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/dfe38fac9a429ce05e2bebcaacbe493f38efd079)]: added detail to tutorial about automatically adding nicks when calling add_full_crossover (David Doty) [#119](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/119)
- [[068162a](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/068162ad1be49add7eea11b15100f20b257f878c)]: changed title in README (David Doty) [#119](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/119)
- [[ae2c223](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/ae2c223c29c7401c4a3a00719da17b72e11b0c0c)]: Update README.md (David Doty) [#119](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/119)
- [[96fadf9](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/96fadf994c714bf29052a524fd47326cb98dcc19)]: removed incorrect word "literal" from docstrings (David Doty) [#119](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/119)
- [[7266c62](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/7266c62808f9c4c4700d51bf38cd4090c36827cc)]: corrected reference to StrandBuilder methods in docstrings (David Doty) [#119](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/119)
- [[97e133c](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/97e133cdb24c493b1146ed7fd4a296ff61a9ed36)]: updated README to be explicit about what steps are needed for pip install and what is troubleshooting (David Doty) [#119](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/119)
- [[cab05a4](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/cab05a483685642392099e6a4d918d2a77d8c6d2)]: added link to PYTHONPATH documentation in README (David Doty) [#119](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/119)
- [[3e17fb0](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/3e17fb041100d94b914597f7ead9433bd0ea446e)]: separated out installation instructions for Python from pip in README (David Doty) [#119](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/119)
- [[b15a33d](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/b15a33def5271e9dc204fadabb6dbdb613e9978c)]: Update README.md (David Doty) [#119](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/119)
- [[256bcb0](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/256bcb092fd0e2678329307cc391cdc3268efd02)]: bumped version (David Doty) [#119](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/119)

[Changes][v0.10.1]


<a id="v0.10.0"></a>
## [v0.10.0](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/releases/tag/v0.10.0) - 2020-07-19

**BREAKING CHANGE**: The name of the class `DNADesign` has been changed to `Design`. You will have to update scripts to reflect this new name. From now on, it would be good practice to name local variables `design` instead of `dna_design`, for instance. I found it was common for novices to abbreviate these as only `dna`, which is confusing since that could be DNA sequence, for example.

Also, the default file extension has been changed from `.dna` to `.sc`. Most short file extensions are taken, but in this case, the most popular for `.dna` is [SnapGene](https://www.snapgene.com/), which is software for biological DNA, and for `.sc`, it's [SuperCollider](https://supercollider.github.io/), a music synthesis tool. There's much less chance that we'll be confused with SuperCollider than SnapGene. However, the web interface still recognizes `.sc`, `.dna`, and `.json` as file extensions for scadnano files that can be opened using the Menu->Load option, or dragged from a file browser onto the web browser with scadnano open.

## Commits
- [[bb8d62d](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/bb8d62de30c7730868db4be253f8ae3311638cc1)]: updated name of function from main() to create_design() in example script in README (David Doty) [#118](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/118)
- [[6db03e8](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/6db03e85314df3e549d345dd5589fb3fd75b9e40)]: changed name of function main() to create_design() in tutorial (David Doty) [#118](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/118)
- [[b4548da](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/b4548daa8d3fefa28ccc9d952751c719e46ca413)]: modified intro to README (David Doty) [#118](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/118)
- [[d0b58d2](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/d0b58d28f06b64be6c28a5e0fc77e7bf8526333e)]: updated tutorial to link directly to raw design file (David Doty) [#118](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/118)
- [[f87b219](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/f87b219bc3e8a93c9560c0092d36e46cdb5c1729)]: fixes [#110](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/issues/110); updated default file extension from .dna to .sc (David Doty) [#118](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/118)
- [[5f9bf45](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/5f9bf4512470652ee4cfc9556a943ee289f279d4)]: closes [#91](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/issues/91); BREAKING CHANGE: change DNADesign class name to Design (David Doty) [#118](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/118)

[Changes][v0.10.0]


<a id="v0.9.10"></a>
## [v0.9.10](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/releases/tag/v0.9.10) - 2020-07-13

Not much changed; this release was to test the automated publishing of the package to PyPI through a new GitHub Action.

## Commits
- [[a0e9f7f](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/a0e9f7fe09d227ad51ec7f81342bd7cf8726f7fe)]: updated installation instructions with detailed pip cases (David Doty) [#115](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/115)
- [[358e389](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/358e3892a8d517445fdeade33515d8e7630657b5)]: corrected grammar in README (David Doty) [#115](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/115)
- [[d8249ae](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/d8249ae549a091bbf2f4de90fcd969d581b8024a)]: updated README instructions on using Python 3.6 (David Doty) [#115](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/115)
- [[4992423](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/4992423ee51b46f10b9be534689a978931450fd2)]: Setup workflow to publish to PyPI (Benjamin Lee) [#114](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/114)
- [[a3e6b1c](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/a3e6b1cb4562d9bd2214ed59fff369ae6fce2284)]: bumped version to test automated package updating on PyPI (David Doty) [#115](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/115)
- [[b3a880e](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/b3a880e2238b2d622cb175bef61f4eef23fefec4)]: Append python publish to release workflow (Benjamin Lee) [#116](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/116)

[Changes][v0.9.10]


<a id="v0.9.9"></a>
## [v0.9.9](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/releases/tag/v0.9.9) - 2020-07-11

This release supports Python version 3.6, if the [dataclasses backport package](https://pypi.org/project/dataclasses/) is also installed. See the [installation instructions](https://github.com/UC-Davis-molecular-computing/scadnano-python-package#installation) for details.

## Commits
- [[c2b3b1d](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/c2b3b1d4a76a54151cf38f1a1a399896438e7945)]: removed `<ins>` tag from README (for underlining name of paper) since PyPI doesn't recognize it (David Doty) [#112](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/112)
- [[54f0e12](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/54f0e12b782d687f58341718bc9df4a3445e488c)]: added note to README that relative links do not work when viewing on PyPI (David Doty) [#112](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/112)
- [[385fb29](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/385fb29d114ee81847fff81fd7f842ff13e5a29a)]: added mypy directory to .gitignore, added squarenut design to examples, updated README on using Python 3.6, and fixed a few mypy warnings in scadnano.py (David Doty) [#112](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/112)
- [[8be2817](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/8be2817d4a54a7c6efac99b3fe439815f002eed4)]: to support Python 3.6, replaced forward annotations with pre-3.7 string-based hack; bumped version (David Doty) [#112](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/112)
- [[054a8bf](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/054a8bf0f96bbd66e13050cce5c262da3482d0ec)]: updated examples to v0.9.9 (David Doty) [#112](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/112)

[Changes][v0.9.9]


<a id="v0.9.8"></a>
## [v0.9.8](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/releases/tag/v0.9.8) - 2020-07-09

Forgot to bump version number in source code in last release, so doing it now.

## Commits
- [[4a3abae](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/4a3abae1681173a133c51b7aaca01920ddc65c82)]: bumped version (David Doty) [#106](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/106)

[Changes][v0.9.8]


<a id="v0.9.7"></a>
## [v0.9.7](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/releases/tag/v0.9.7) - 2020-07-09

Just trying to test out some new repo settings. Note the version number is not updated in the source code because I forgot. The next version will fix that.

## Commits
- [[57a1eae](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/57a1eae5da6bdd7298de5b8912f1cc8a1674951c)]: changed comment on rise_per_base_pair (David Doty) [#105](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/105)
- [[a081e59](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/a081e5956f9691b47a81799f3fc9c5fbf3ebae89)]: fixed newline in README (David Doty) [#105](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/105)
- [[12ff50d](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/12ff50d585ea338560d33b974428b8002b6d0d14)]: removed requirements.txt (David Doty) [#105](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/105)
- [[a4d2bcd](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/a4d2bcd520bb0b922b707e531bfb1b31530a0452)]: added CODEOWNERS for protected branches PR review requests (David Doty) [#105](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/105)
- [[dff2ec8](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/dff2ec85530152ef24de6017fea5dd0429b94ef1)]: added very large origami example (David Doty) [#105](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/105)
- [[58bb6bc](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/58bb6bcfd7ca0274afa54a65aed63c5304b76b8f)]: fixed example code in README to reference design.scaffold (David Doty) [#105](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/105)
- [[85d0f9d](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/85d0f9d716b712b1cc5af4eff4898c5414cdeeeb)]: changed to relative links to Python files in installation by download method in README, so viewing the README in another branch will link to the version in that branch (David Doty) [#105](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/105)

[Changes][v0.9.7]


<a id="v0.9.6"></a>
## [v0.9.6](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/releases/tag/v0.9.6) - 2020-07-03

## Commits
- [[24ff64f](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/24ff64f555d558d00f48664adcacb3170c2d435e)]: bumped version and updated doc/requirements.txt in response to dependabot alerts (David Doty) [#103](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/103)

[Changes][v0.9.6]


<a id="v0.9.5"></a>
## [v0.9.5](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/releases/tag/v0.9.5) - 2020-07-03

Based on this: readthedocs/readthedocs.org#4868, seems like "latest" tag is confusing readthedocs. This release fixes that.

[Changes][v0.9.5]


<a id="v0.9.4"></a>
## [v0.9.4](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/releases/tag/v0.9.4) - 2020-07-03

## Commits
- [[128baa0](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/128baa0302122640d1835bd50ce4a393c8878e6f)]: ignoring dist/ directory from now on (where PyPI .tar.gz files are created) (David Doty) [#98](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/98)
- [[0b830a2](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/0b830a2d78666ea9c30f5798533006c3d286403c)]: removed dist/ directory from tracking (David Doty) [#98](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/98)
- [[c3e2c12](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/c3e2c12406b105c42f557e0cb530312250ed4d5a)]: add with_domain_label for method chaining; closes [#93](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/issues/93) (David Doty) [#98](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/98)
- [[18ba037](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/18ba0375480a9cc41dc1f1b517a4a604063f4bdf)]: corrected README link to write_idt_plate_excel_file (David Doty) [#98](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/98)
- [[8d86dff](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/8d86dff1188c408ad11707955e4f366b618ce2cc)]: updated citation for paper from arXiv to DNA 26 (David Doty) [#98](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/98)
- [[fed5149](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/fed5149271e0f6183eb544588dcbf61c312544ba)]: bumped version (David Doty) [#98](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/98)

[Changes][v0.9.4]


<a id="v0.9.3"></a>
## [v0.9.3](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/releases/tag/v0.9.3) - 2020-06-28

## Commits
- [[c4db84a](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/commit/c4db84a647e5cd3528347996b97a4ebae2a29a0f)]: bumped version from 0.9.2 to 0.9.3 (David Doty) [#97](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/pull/97)

[Changes][v0.9.3]


<a id="v0.9.2"></a>
## [v0.9.2](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/releases/tag/v0.9.2) - 2020-06-28

This is the first release that has automated release notes, but since the previous one didn't they weren't automatically generated. Future releases should be automatically documented here. I'm manually pasting in changes from PR [#96](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/issues/96) 

[@dave-doty](https://github.com/dave-doty)
changed name of parameters in Crossover, add_half_crossover, and add_… …
6a4d537
…full_crossover
[@dave-doty](https://github.com/dave-doty)
fixed name shadowing bug
0a93dbd
[@dave-doty](https://github.com/dave-doty)
cleaned up docstrings
ffdd62f
[@dave-doty](https://github.com/dave-doty)
cleaned up docstrings
d33dbae
[@dave-doty](https://github.com/dave-doty)
replaced parameter names in calls to add_*_crossover
9dc504f
[@dave-doty](https://github.com/dave-doty)
added code to check for overlapping domains in post_init of Strand
0426e10
[@dave-doty](https://github.com/dave-doty)
fixes [#85](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/issues/85): allow chained commands for less verbose way to add Strands to DNADesign
f36d223
[@dave-doty](https://github.com/dave-doty)
added a few chained methods for modifying strands
534fc8a
[@dave-doty](https://github.com/dave-doty)
updated docstrings
aa0e16d
[@dave-doty](https://github.com/dave-doty)
cleaned up docstrings
bfe6150
[@dave-doty](https://github.com/dave-doty)
cleaned up some PEP warnings
7ef61cf
[@dave-doty](https://github.com/dave-doty)
closes [#87](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/issues/87): add chained command for adding new domain on current helix
b10e6ee
[@dave-doty](https://github.com/dave-doty)
fixes [#92](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/issues/92): make ColorCycler part of DNADesign
5b001fa
[@dave-doty](https://github.com/dave-doty)
added Kelly colors, but not using them yet
700c95a
[@dave-doty](https://github.com/dave-doty)
closes [#84](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/issues/84): add Geometry parameters field to DNADesign
033529f
[@dave-doty](https://github.com/dave-doty)
closes [#90](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/issues/90): custom extension in write_scadnano_file
2849caa
[@dave-doty](https://github.com/dave-doty)
closes [#86](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/issues/86): support domain labels, strand labels
6c67b5a
[@dave-doty](https://github.com/dave-doty)
closes [#88](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/issues/88): add StrandBuilder.with_domain_sequence
0462097
[@dave-doty](https://github.com/dave-doty)
updated docstrings
fb38cb9
[@dave-doty](https://github.com/dave-doty)
cleaned up docstrings
b6bc810
[@dave-doty](https://github.com/dave-doty)
removed note about web interface in tutorial since none of the menus … …
8c66ac0
…are shown
[@dave-doty](https://github.com/dave-doty)
moved Python code from web interface README to here
30a60db
[@dave-doty](https://github.com/dave-doty)
minor
c8b3a21
[@dave-doty](https://github.com/dave-doty)
reordered sections
af09199
[@dave-doty](https://github.com/dave-doty)
changed section title
056ad56
[@dave-doty](https://github.com/dave-doty)
changed name of z_step to rise_per_base_pair.
9f9a013
[@dave-doty](https://github.com/dave-doty)
modified example with many DNA modifications to have every combinatio… …
1dc9479
…n of modifications on forward and reverse strands, on both helices, to help with testing "invert y axis" view option in web interface.
[@UnHumbleBen](https://github.com/UnHumbleBen)
Create release.yml
3d4f93b
[@dave-doty](https://github.com/dave-doty)
changed name of z_step to rise_per_base_pair in unit tests
5742687

[Changes][v0.9.2]


<a id="v0.9.1"></a>
## [v0.9.1](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/releases/tag/v0.9.1) - 2020-06-20

removed `forward` and `reverse` helper variables from scadnano. So you'll have to replace this:

```python
import scadnano as sc

domain1 = sc.Domain(0, sc.forward, 0, 10)
domain2 = sc.Domain(0, sc.reverse, 0, 10)
```

with this

```python
import scadnano as sc

domain1 = sc.Domain(0, True, 0, 10)
domain2 = sc.Domain(0, False, 0, 10)
```


[Changes][v0.9.1]


<a id="v0.9.0"></a>
## [v0.9.0](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/releases/tag/v0.9.0) - 2020-06-13

***Breaking change:***
Changed name of Helix.position3d to Helix.position

Changed documentation explaining meaning of Helix.position to match new interpretation in scadnano of how it translates to side view and main view display, for issue [#78](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/issues/78), which swaps the role of Helix.position.x and Helix.position.z.

The Python library itself did not do any such translation, but scripts will have to be re-written to change x and z to display properly in the scadnano web interface once issue [UC-Davis-molecular-computing/scadnano#307](https://github.com/UC-Davis-molecular-computing/scadnano/issues/307) is closed.

[Changes][v0.9.0]


<a id="v0.8.3"></a>
## [v0.8.3](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/releases/tag/v0.8.3) - 2020-06-09

Removed need to have xlwt installed if downloading the scadnano.py file and importing locally, unless an Excel file needs to be written.

[Changes][v0.8.3]


<a id="v0.8.0"></a>
## [v0.8.0](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/releases/tag/v0.8.0) - 2020-06-04

Removed support for Helix.rotation and Helix.rotation_anchor.

Removed pitch, roll, and yaw from Helix.position, and made them top-level fields in Helix.

See https://github.com/UC-Davis-molecular-computing/scadnano/blob/master/README.md for description of current data model.

[Changes][v0.8.0]


<a id="v0.7.3"></a>
## [v0.7.3](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/releases/tag/v0.7.3) - 2020-05-31



[Changes][v0.7.3]


<a id="v0.7.2"></a>
## [v0.7.2](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/releases/tag/v0.7.2) - 2020-05-31



[Changes][v0.7.2]


<a id="v0.7.1"></a>
## [v0.7.1](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/releases/tag/v0.7.1) - 2020-05-30



[Changes][v0.7.1]


<a id="v0.7.0"></a>
## [v0.7.0](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/releases/tag/v0.7.0) - 2020-05-30



[Changes][v0.7.0]


<a id="v0.6.8"></a>
## [v0.6.8](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/releases/tag/v0.6.8) - 2020-05-27



[Changes][v0.6.8]


<a id="v0.6.5"></a>
## [v0.6.5](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/releases/tag/v0.6.5) - 2020-05-22



[Changes][v0.6.5]


<a id="v0.6.0"></a>
## [v0.6.0](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/releases/tag/v0.6.0) - 2020-05-21



[Changes][v0.6.0]


<a id="v0.5.0"></a>
## [v0.5.0](https://github.com/UC-Davis-molecular-computing/scadnano-python-package/releases/tag/v0.5.0) - 2020-05-20



[Changes][v0.5.0]


[v0.21.1]: https://github.com/UC-Davis-molecular-computing/scadnano-python-package/compare/v0.21.0...v0.21.1
[v0.21.0]: https://github.com/UC-Davis-molecular-computing/scadnano-python-package/compare/v0.20.1...v0.21.0
[v0.20.1]: https://github.com/UC-Davis-molecular-computing/scadnano-python-package/compare/v0.20.0...v0.20.1
[v0.20.0]: https://github.com/UC-Davis-molecular-computing/scadnano-python-package/compare/v0.19.4...v0.20.0
[v0.19.4]: https://github.com/UC-Davis-molecular-computing/scadnano-python-package/compare/v0.19.3...v0.19.4
[v0.19.3]: https://github.com/UC-Davis-molecular-computing/scadnano-python-package/compare/v0.19.0...v0.19.3
[v0.19.0]: https://github.com/UC-Davis-molecular-computing/scadnano-python-package/compare/v0.18.3...v0.19.0
[v0.18.3]: https://github.com/UC-Davis-molecular-computing/scadnano-python-package/compare/v0.18.2...v0.18.3
[v0.18.2]: https://github.com/UC-Davis-molecular-computing/scadnano-python-package/compare/v0.18.1...v0.18.2
[v0.18.1]: https://github.com/UC-Davis-molecular-computing/scadnano-python-package/compare/v0.18.0...v0.18.1
[v0.18.0]: https://github.com/UC-Davis-molecular-computing/scadnano-python-package/compare/v0.17.7...v0.18.0
[v0.17.7]: https://github.com/UC-Davis-molecular-computing/scadnano-python-package/compare/v0.17.6...v0.17.7
[v0.17.6]: https://github.com/UC-Davis-molecular-computing/scadnano-python-package/compare/v0.17.5...v0.17.6
[v0.17.5]: https://github.com/UC-Davis-molecular-computing/scadnano-python-package/compare/v0.17.3...v0.17.5
[v0.17.3]: https://github.com/UC-Davis-molecular-computing/scadnano-python-package/compare/v0.17.2...v0.17.3
[v0.17.2]: https://github.com/UC-Davis-molecular-computing/scadnano-python-package/compare/v0.16.3...v0.17.2
[v0.16.3]: https://github.com/UC-Davis-molecular-computing/scadnano-python-package/compare/v0.16.2...v0.16.3
[v0.16.2]: https://github.com/UC-Davis-molecular-computing/scadnano-python-package/compare/v0.16.1...v0.16.2
[v0.16.1]: https://github.com/UC-Davis-molecular-computing/scadnano-python-package/compare/v0.16.0...v0.16.1
[v0.16.0]: https://github.com/UC-Davis-molecular-computing/scadnano-python-package/compare/v0.15.1...v0.16.0
[v0.15.1]: https://github.com/UC-Davis-molecular-computing/scadnano-python-package/compare/v0.15.0...v0.15.1
[v0.15.0]: https://github.com/UC-Davis-molecular-computing/scadnano-python-package/compare/v0.14.0...v0.15.0
[v0.14.0]: https://github.com/UC-Davis-molecular-computing/scadnano-python-package/compare/v0.13.4...v0.14.0
[v0.13.4]: https://github.com/UC-Davis-molecular-computing/scadnano-python-package/compare/v0.13.0...v0.13.4
[v0.13.0]: https://github.com/UC-Davis-molecular-computing/scadnano-python-package/compare/v0.12.2...v0.13.0
[v0.12.2]: https://github.com/UC-Davis-molecular-computing/scadnano-python-package/compare/v0.12.1...v0.12.2
[v0.12.1]: https://github.com/UC-Davis-molecular-computing/scadnano-python-package/compare/v0.12.0...v0.12.1
[v0.12.0]: https://github.com/UC-Davis-molecular-computing/scadnano-python-package/compare/v0.11.2...v0.12.0
[v0.11.2]: https://github.com/UC-Davis-molecular-computing/scadnano-python-package/compare/v0.11.1...v0.11.2
[v0.11.1]: https://github.com/UC-Davis-molecular-computing/scadnano-python-package/compare/v0.11.0...v0.11.1
[v0.11.0]: https://github.com/UC-Davis-molecular-computing/scadnano-python-package/compare/v0.10.3...v0.11.0
[v0.10.3]: https://github.com/UC-Davis-molecular-computing/scadnano-python-package/compare/v0.10.2...v0.10.3
[v0.10.2]: https://github.com/UC-Davis-molecular-computing/scadnano-python-package/compare/v0.10.1a...v0.10.2
[v0.10.1a]: https://github.com/UC-Davis-molecular-computing/scadnano-python-package/compare/v0.10.1...v0.10.1a
[v0.10.1]: https://github.com/UC-Davis-molecular-computing/scadnano-python-package/compare/v0.10.0...v0.10.1
[v0.10.0]: https://github.com/UC-Davis-molecular-computing/scadnano-python-package/compare/v0.9.10...v0.10.0
[v0.9.10]: https://github.com/UC-Davis-molecular-computing/scadnano-python-package/compare/v0.9.9...v0.9.10
[v0.9.9]: https://github.com/UC-Davis-molecular-computing/scadnano-python-package/compare/v0.9.8...v0.9.9
[v0.9.8]: https://github.com/UC-Davis-molecular-computing/scadnano-python-package/compare/v0.9.7...v0.9.8
[v0.9.7]: https://github.com/UC-Davis-molecular-computing/scadnano-python-package/compare/v0.9.6...v0.9.7
[v0.9.6]: https://github.com/UC-Davis-molecular-computing/scadnano-python-package/compare/v0.9.5...v0.9.6
[v0.9.5]: https://github.com/UC-Davis-molecular-computing/scadnano-python-package/compare/v0.9.4...v0.9.5
[v0.9.4]: https://github.com/UC-Davis-molecular-computing/scadnano-python-package/compare/v0.9.3...v0.9.4
[v0.9.3]: https://github.com/UC-Davis-molecular-computing/scadnano-python-package/compare/v0.9.2...v0.9.3
[v0.9.2]: https://github.com/UC-Davis-molecular-computing/scadnano-python-package/compare/v0.9.1...v0.9.2
[v0.9.1]: https://github.com/UC-Davis-molecular-computing/scadnano-python-package/compare/v0.9.0...v0.9.1
[v0.9.0]: https://github.com/UC-Davis-molecular-computing/scadnano-python-package/compare/v0.8.3...v0.9.0
[v0.8.3]: https://github.com/UC-Davis-molecular-computing/scadnano-python-package/compare/v0.8.0...v0.8.3
[v0.8.0]: https://github.com/UC-Davis-molecular-computing/scadnano-python-package/compare/v0.7.3...v0.8.0
[v0.7.3]: https://github.com/UC-Davis-molecular-computing/scadnano-python-package/compare/v0.7.2...v0.7.3
[v0.7.2]: https://github.com/UC-Davis-molecular-computing/scadnano-python-package/compare/v0.7.1...v0.7.2
[v0.7.1]: https://github.com/UC-Davis-molecular-computing/scadnano-python-package/compare/v0.7.0...v0.7.1
[v0.7.0]: https://github.com/UC-Davis-molecular-computing/scadnano-python-package/compare/v0.6.8...v0.7.0
[v0.6.8]: https://github.com/UC-Davis-molecular-computing/scadnano-python-package/compare/v0.6.5...v0.6.8
[v0.6.5]: https://github.com/UC-Davis-molecular-computing/scadnano-python-package/compare/v0.6.0...v0.6.5
[v0.6.0]: https://github.com/UC-Davis-molecular-computing/scadnano-python-package/compare/v0.5.0...v0.6.0
[v0.5.0]: https://github.com/UC-Davis-molecular-computing/scadnano-python-package/tree/v0.5.0

<!-- Generated by https://github.com/rhysd/changelog-from-release v3.9.1 -->
