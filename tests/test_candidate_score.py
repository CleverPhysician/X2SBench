import copy
import unittest
from pymatgen.core import Lattice, Structure
from pymatgen.io.cif import CifWriter
from candidate_score import score_target

def fixture():
    s = Structure(Lattice.cubic(4), ['Cs+', 'Cl-'], [[0,0,0],[.5,.5,.5]])
    candidate = dict(lengths=[4,4,4], angles=[90,90,90], atom_types=[55,17],
                     frac_coords=[[0,0,0],[.5,.5,.5]])
    return dict(material_id='toy', index=0, cif=str(CifWriter(s)),
                candidates=[candidate], budgets=[1])

class CandidateScoreTest(unittest.TestCase):
    def test_element_equivalent_perfect_prediction_is_matched(self):
        result = score_target(fixture())
        self.assertTrue(result.get('topk', [{}])[0].get('matched'))
        self.assertAlmostEqual(result['topk'][0]['best_rmsd'], 0.0, places=8)

    def test_failure_is_preserved_before_later_success(self):
        item=fixture(); bad=copy.deepcopy(item['candidates'][0]);bad['atom_types']=[55,9]
        item['candidates'].insert(0,bad);item['budgets']=[1,2]
        result=score_target(item)
        self.assertEqual([r['matched'] for r in result.get('topk',[])],[False,True])
        self.assertEqual(len(result['candidates']),2)

    def test_empty_candidate_is_not_removed_from_budget(self):
        item=fixture();item['candidates'][0].update(atom_types=[],frac_coords=[])
        result=score_target(item)
        self.assertEqual(result.get('topk',[{}])[0].get('matched'),False)
        self.assertEqual(result['candidates'][0]['valid'],False)

    def test_nonfinite_candidate_is_retained_as_invalid(self):
        item=fixture();item['candidates'][0]['frac_coords'][0][0]=float('nan')
        self.assertFalse(score_target(item).get('candidates',[{'valid':True}])[0]['valid'])

    def test_invalid_lattice_is_rejected_before_matching(self):
        item=fixture();item['candidates'][0]['angles'][0]=0
        self.assertFalse(score_target(item).get('candidates',[{'valid':True}])[0]['valid'])

    def test_malformed_packing_raises_instead_of_changing_denominator(self):
        item=fixture();item['candidates'][0]['frac_coords'].pop()
        with self.assertRaises(ValueError):score_target(item)

    def test_disordered_reference_stops_scoring(self):
        item=fixture();s=Structure(Lattice.cubic(4),[{'Cs':.5,'Rb':.5},'Cl'],[[0,0,0],[.5,.5,.5]])
        item['cif']=str(CifWriter(s))
        with self.assertRaises(ValueError):score_target(item)

    def test_periodically_equivalent_larger_cell_is_not_rejected_by_atom_count(self):
        item=fixture();s=Structure(Lattice.cubic(4),['Cs','Cl'],[[0,0,0],[.5,.5,.5]])
        s.make_supercell([2,1,1])
        item['candidates'][0]=dict(lengths=s.lattice.abc,angles=s.lattice.angles,
                                    atom_types=[a.specie.Z for a in s],frac_coords=s.frac_coords.tolist())
        self.assertTrue(score_target(item).get('topk',[{}])[0].get('matched'))

    def test_prediction_scoring_retains_the_common_volume_scaling_policy(self):
        item=fixture();item['candidates'][0]['lengths']=[5,5,5]
        self.assertTrue(score_target(item).get('topk',[{}])[0].get('matched'))

    def test_reference_strata_do_not_follow_prediction_size(self):
        item=fixture();item['candidates'][0]['atom_types']=[55]
        item['candidates'][0]['frac_coords']=[[0,0,0]]
        result=score_target(item)
        self.assertEqual(result.get('num_atoms'),2)
        self.assertEqual(result.get('num_elements'),2)
        self.assertEqual(result.get('crystal_system'),'cubic')

if __name__=='__main__':unittest.main()
