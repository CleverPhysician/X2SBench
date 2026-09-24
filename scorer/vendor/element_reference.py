"""Element-only references: preserve geometry and occupancies, never relax tolerances."""
import numpy as np
from pymatgen.core import Element, Structure
from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

REFERENCE_POLICY = 'element_identity_v1_symprec0.01_angle5'


def element_reference(structure):
    if not structure.is_ordered:
        raise ValueError('Disordered/partial-occupancy reference needs explicit handling; not converted to full occupancy')
    original_z = [int(site.specie.Z) for site in structure]
    original_occupancies = [tuple(site.species.values()) for site in structure]
    result = structure.copy()
    result.remove_oxidation_states()
    if not all(isinstance(site.specie, Element) for site in result):
        raise ValueError('Non-element species after standardization')
    assert np.array_equal(result.lattice.matrix, structure.lattice.matrix)
    assert np.array_equal(result.frac_coords, structure.frac_coords)
    assert original_z == [int(site.specie.Z) for site in result]
    assert original_occupancies == [tuple(site.species.values()) for site in result]
    return result


def crystal_system(structure):
    return SpacegroupAnalyzer(structure, symprec=0.01, angle_tolerance=5).get_crystal_system()


def identity_check(structure):
    # Construct the exact geometry through the model's atomic-number representation.
    candidate = Structure(structure.lattice, [int(s.specie.Z) for s in structure],
                          structure.frac_coords)
    match = StructureMatcher(stol=0.5, angle_tol=10, ltol=0.3, scale=True).get_rms_dist(candidate, structure)
    # Primitive-cell reduction can leave ~6e-6 numerical residual for C8.
    # This near-zero diagnostic threshold does not alter StructureMatcher tolerances.
    if match is None or not np.isfinite(match[0]) or match[0] > 1e-4:
        raise AssertionError(f'Exact model-representation reference failed identity check: formula={structure.composition}, match={match}')
    return float(match[0])
