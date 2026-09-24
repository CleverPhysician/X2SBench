# X2SBench

Evaluation resources for **X2SBench: an open benchmark for evaluating crystal structure determination from powder diffraction**.

This initial repository contains the element-standardized candidate scorer, its regression tests and the fixed 558-target experimental evaluation list. The accompanying experimental data archive contains 591 exact-profile-deduplicated measurements and their archived reference structures. The simulated dataset and model-specific inference runners are not included in this initial package.

**[Download experimental data (591 records, approximately 42 MB)](https://github.com/CleverPhysician/X2SBench/releases/download/v0.1.0-experimental/experimental_591.zip)** · [Release notes](https://github.com/CleverPhysician/X2SBench/releases/tag/v0.1.0-experimental)

## Included

- `scorer/candidate_score.py`: scores packed candidate structures against CIF references and retains failed candidates in their sampling positions.
- `scorer/vendor/element_reference.py`: element-only species normalization and reference crystal-system assignment.
- `scorer/archive_layout.py`: candidate-budget aggregation and archived packing helpers.
- `tests/test_candidate_score.py`: ten existing regression checks copied with the scoring implementation.
- `splits/evaluation_558_ids.json`: ordered identifiers of the common experimental cohort.

The scorer is copied unchanged from the manuscript's archived common-rescoring implementation. It uses pymatgen StructureMatcher with lattice tolerance 0.3, site tolerance 0.5, angle tolerance 10 degrees, primitive-cell reduction and volume scaling. RMSD is dimensionless and is averaged over successfully matched targets. Reference strata use crystal system, reference-cell atom count and distinct element count.

## Smoke test

With Python 3.12 in a virtual environment:

```sh
python -m pip install -r requirements.txt
PYTHONPATH=scorer python -m unittest discover -s tests -p 'test_*.py'
```

The pinned versions match the local regression-test environment: pymatgen 2024.5.1, spglib 2.7.0 and NumPy 1.26.4. Tests cover species equivalence, invalid candidates, candidate order, periodic cell equivalence and reference-based stratification. A fresh installation and complete reproduction of manuscript results remain to be validated.

## Scoring API

Add `scorer` to `PYTHONPATH`, then call `candidate_score.score_target(item)`. An item contains:

- `material_id`: target identifier; `index`: target index.
- `cif`: reference CIF text, used for scoring only.
- `candidates`: an ordered list of structures, each with `lengths` (three values in angstrom), `angles` (three values in degrees), `atom_types` (atomic numbers) and `frac_coords` (N-by-3 fractional coordinates).
- `budgets`: candidate budgets, for example `[1, 5, 10, 20]` for 20 candidates.

Keep failed generations in their original sampling positions. An empty candidate uses empty `atom_types` and `frac_coords` with three placeholder lengths and angles. Malformed packing raises an error rather than silently removing the target. See the synthetic structures in `tests/test_candidate_score.py` for executable examples.

The result contains candidate-level validity and RMSD, best-of-k matches, and reference crystal system, atom count and element count. Match rate is the fraction of targets with at least one match within the budget. Average RMSD is the mean of each matched target's best RMSD; unmatched targets remain in the match-rate denominator and are excluded from conditional RMSD. Report matched counts alongside RMSD.

This is a scoring library, not an inference service or an end-to-end leaderboard client. Input conditions, reference-cell conventions, training history and reference-lattice access must be specified for comparisons.

## Data and attribution

See [experimental data](docs/experimental_data.md) for loading instructions and field definitions. The measured collection draws on RRUFF and opXRD. Preserve source attribution and applicable source terms. No new license is assigned to third-party data, and a software license has not yet been selected by the authors.

The benchmark platform is [X2SBench at CMPDC](https://cmpdc.iphy.ac.cn/benchmarks/detail/2f). The experimental archive is versioned in [GitHub Releases](https://github.com/CleverPhysician/X2SBench/releases/tag/v0.1.0-experimental). A verified manuscript identifier will be added when available.
