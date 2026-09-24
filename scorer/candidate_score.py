"""Common element-only scoring of variable-size archived crystal candidates."""
import warnings
import numpy as np
from pymatgen.core import Lattice, Structure
from pymatgen.analysis.structure_matcher import StructureMatcher
from vendor.element_reference import element_reference, crystal_system
from archive_layout import summarize_candidates

def score_target(item):
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        reference = element_reference(Structure.from_str(item['cif'], fmt='cif'))
    matcher = StructureMatcher(stol=.5, angle_tol=10, ltol=.3, scale=True)
    results, distances = [], []
    for i, c in enumerate(item['candidates']):
        types = c['atom_types']
        xyz = np.asarray(c['frac_coords'], dtype=float)
        if len(xyz) != len(types) or (len(types) and xyz.shape != (len(types), 3)):
            raise ValueError('Malformed candidate coordinate/species packing')
        lengths, angles = np.asarray(c['lengths']), np.asarray(c['angles'])
        if lengths.shape != (3,) or angles.shape != (3,):
            raise ValueError('Malformed candidate lattice packing')
        distance, valid, error = None, False, ''
        try:
            if not types:
                raise ValueError('empty_candidate')
            if not (np.isfinite(lengths).all() and np.isfinite(angles).all() and np.isfinite(xyz).all()):
                raise ValueError('nonfinite')
            if lengths.min() < .5 or angles.min() < 5 or angles.max() > 175:
                raise ValueError('invalid_lattice')
            candidate = element_reference(Structure(Lattice.from_parameters(*lengths, *angles), types, xyz))
            if candidate.volume < .1:
                raise ValueError('invalid_volume')
            valid = True
            match = matcher.get_rms_dist(candidate, reference)
            if match is not None:
                distance = float(match[0])
                if not np.isfinite(distance) or distance < 0:
                    raise ValueError('invalid_match_distance')
        except Exception as exc:
            distance = None
            error = type(exc).__name__+': '+str(exc)[:160]
        distances.append(distance)
        results.append(dict(candidate=i+1, candidate_num_atoms=len(types), valid=valid,
                            rmsd=distance, error=error))
    return dict(material_id=item['material_id'], index=item['index'],
                num_atoms=len(reference), num_elements=len(reference.composition.elements),
                crystal_system=crystal_system(reference), formula=reference.composition.formula,
                parser_warning_count=len(caught), candidates=results,
                topk=summarize_candidates(distances, item['budgets']))
