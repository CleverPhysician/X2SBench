# Experimental data

## Package and cohort

[Download the versioned experimental archive](https://github.com/CleverPhysician/X2SBench/releases/download/v0.1.0-experimental/experimental_591.zip). SHA-256 of the ZIP archive: `2e0e81f8413dbee9e18e35dddf6383b09e7f3d8e6547e5fa276b568843fbd604`.

`experimental_591_exact_deduplicated.pkl` contains 591 records. Every profile has 7,501 pairs of two-theta angle and intensity. Archived values are preserved without new smoothing, normalization or resampling. SHA-256 of the uncompressed file:

```text
3ec2b1b358c00d20e38f14b6b0aaa92d21837c677fc4f4f5922fc2dffd08083c
```

The common evaluation uses the 558 identifiers in `splits/evaluation_558_ids.json`; 33 larger structures are outside this cohort. Use the ordered identifier list rather than row positions from earlier dataset versions.

## Loading

Unpack the experimental archive beside this repository. Pickle files can execute code: load only a trusted release after checking its checksum.

```python
import json
import pickle
import numpy as np

with open("experimental_591_exact_deduplicated.pkl", "rb") as f:
    records = pickle.load(f)
with open("splits/evaluation_558_ids.json") as f:
    ids = json.load(f)
lookup = {record["material_id"]: record for record in records}
evaluation = [lookup[identifier] for identifier in ids]
profile = np.asarray(evaluation[0]["xrd"]["xrd_plot_data"], dtype=float)
two_theta, intensity = profile[:, 0], profile[:, 1]
```

The measured full profile is **`xrd.xrd_plot_data`**. Separate `xrd.two_theta`, `xrd.intensity`, `hkls` and `d_hkls` fields are archived peak-list metadata, not the full measured profile.

Records also include `material_id`, composition metadata, `conventional.cif`, `primitive.cif` and lattice metadata. Explicitly choose the reference representation required by the evaluation protocol. Reference CIFs are scoring targets, not model inputs in the no-reference-lattice task.

Legacy `spacegroup` fields are retained for provenance. Recompute evaluation symmetry using the scorer's element-only representation rather than grouping by legacy species labels. CIF text is unchanged by packaging; species normalization is performed by the scorer.

## Measurement conventions and sources

The nominal grid is 5–80 degrees in two-theta with 0.01-degree spacing and a Cu K-alpha convention. Some wavelength assignments are inherited or assumed in the curated sources. Packaging preserves those assignments and does not independently verify every pattern–structure pairing. Deduplication keeps one retained record per exact-profile group; alternative phase assignments require source-level checking.

The measured collection draws on RRUFF and opXRD:

- Lafuente et al., *The power of databases: the RRUFF project* (2015), [RRUFF](https://rruff.info/).
- Hollarek et al., *opXRD: Open Experimental Powder X-ray Diffraction Database*, [dataset record](https://zenodo.org/records/15298026), DOI: 10.5281/zenodo.15298026.

These are source-resource citations, not an assertion that every archived record came from the cited opXRD version. Preserve source attribution and consult applicable source terms for redistribution and reuse. No new blanket data license is assigned here. Source-specific provenance is not reconstructed where archived records do not supply it.
