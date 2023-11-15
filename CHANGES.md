# Changelog

## v1.2.0

- Molecules as "big supercells" is now implemented using `Database.from_files_molecule`.
- Fixed a bug with the Materials Project API where `mpr.materials.summary.search` no longer works. This was changed to `mpr.materials.search`, which seems to provide the same functionality.

## v1.1.0

- Added compatibility with the [new Materials Project API (v2)](https://next-gen.materialsproject.org/api). Old users of Lightshow will have to update their API key and might notice that some materials that were previously available no longer are.
- Fixed a bug ([#199](https://github.com/AI-multimodal/Lightshow/issues/199)) that in rare ocassions caused mismatches between the supercell and cluster input files.
