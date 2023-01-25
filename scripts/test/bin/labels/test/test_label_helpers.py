#!/nfs/turbo/umms-petefred/jaschwa/.conda/HDPRED/bin/python
import sys
import unittest
import numpy as np
from tools.alt_calcs import calc_dihedral, dist, calc_angle_b

sys.path.append('..')
import labels as label_funcs
        
        
lines_example = [ \
    'ATOM      1  N   LEU A  12      39.822 125.897  22.200', \
    'ATOM      2  CA  LEU A  12      40.423 125.448  23.495', \
    'ATOM      3  C   LEU A  12      40.211 126.555  24.528', \
    'ATOM      4  O   LEU A  12      41.112 127.348  24.785', \
    'ATOM      5  CB  LEU A  12      41.928 125.188  23.299', \
    'ATOM      6  CG  LEU A  12      42.637 124.139  24.174', \
    'ATOM      7  CD1 LEU A  12      42.527 124.511  25.642', \
    'ATOM      8  CD2 LEU A  12      42.028 122.766  23.922', \
    'ATOM      9  N   ASP A  13      39.017 126.605  25.114', \
    'ATOM     10  CA  ASP A  13      38.688 127.621  26.119', \
    'ATOM     11  C   ASP A  13      38.823 127.066  27.536', \
    'ATOM     12  O   ASP A  13      39.293 125.948  27.723', \
    'ATOM     13  CB  ASP A  13      37.254 128.143  25.914', \
    'ATOM     14  CG  ASP A  13      37.114 129.043  24.688', \
    'ATOM     15  OD1 ASP A  13      37.483 130.236  24.768', \
    'ATOM     16  OD2 ASP A  13      36.632 128.559  23.641', \
    'ATOM     17  N   ARG A  14      38.402 127.843  28.530', \
    'ATOM     18  CA  ARG A  14      38.479 127.401  29.920', \
    'ATOM     19  CB  ARG A  14      37.094 127.002  30.406', \
    'ATOM     20  O   ARG A  14      36.815 125.830  30.648'  \
]


class TestLabelHelpers(unittest.TestCase):
    
    def test_read_atom_lines(self):
        lines = label_funcs.read_atom_lines('data/small_example.pdb')
        self.assertEqual(len(lines), 13)

    def test_collect_residues_from_lines(self):
        residues = label_funcs.collect_residues_from_lines(lines_example)
        self.assertEqual(len(residues), 3)

        res12 = residues[0]
        self.assertEqual(len(res12.keys()), 9)
        self.assertEqual(res12['res_name'], 'LEU')
        self.assertTrue(np.array_equal(res12['CA'], np.array([[40.423, 125.448, 23.495]])))

        res14 = residues[-1]
        self.assertEqual(len(res14.keys()), 5)
        self.assertEqual(res14['res_name'], 'ARG')
        self.assertTrue(np.array_equal(res14['O'], np.array([[36.815, 125.830, 30.648]])))

    def test_filter_valid_residues_removeCA(self):
        lines_copy = lines_example.copy()
        lines_copy.pop(9) # remove res13's CA
        residues = label_funcs.collect_residues_from_lines(lines_copy)

        # should only get two residues after filter
        residues = label_funcs.filter_valid_residues(residues)
        self.assertEqual(len(residues), 2)
        self.assertEqual(residues[0]['res_name'], 'LEU')
        self.assertEqual(residues[1]['res_name'], 'ARG')

    def test_filter_valid_residues_exception(self):
        # now we want to confirm an exception occurs for a missing CB
        lines_copy = lines_example.copy()
        lines_copy.pop(18) # res14's CB
        residues = label_funcs.collect_residues_from_lines(lines_copy)
        self.assertRaises(ValueError, lambda: label_funcs.filter_valid_residues(residues))
          

    def test_convert_to_resdic_format(self):
        residues = label_funcs.collect_residues_from_lines(lines_example)
        resdic, nums, seq = label_funcs.convert_to_resdic_format(residues)
        
        self.assertEqual(len(resdic), 20)
        self.assertEqual(len(nums), 3)
        self.assertEqual(len(seq), 3)
        self.assertEqual(seq, ['LEU', 'ASP', 'ARG'])

    def test_read_chain(self):
        resdic, nums, seq = label_funcs.read_chain('data/large_example.pdb')
        self.assertEqual(len(resdic), 1053)
        self.assertEqual(len(nums), 153)
        self.assertEqual(len(seq), 153)
        self.assertEqual(nums[0], 0)
        self.assertEqual(nums[-1], 152)
        self.assertEqual(seq[:3], 'LDR')
        self.assertEqual(seq[-3:], 'KTR')

    def test_compute_dis(self):
        resdic, _, _ = label_funcs.read_chain('data/small_example.pdb')
        a = resdic['0_CA']
        b = resdic['1_O']
        dis = dist(a, b)
        self.assertEqual(dis, label_funcs.compute_dis(0, 1, resdic, 'CA', 'O'))

    def test_compute_omg(self):
        resdic, _, _ = label_funcs.read_chain('data/small_example.pdb')
        a = resdic['0_CA'].flatten()
        b = resdic['0_CB'].flatten()
        c = resdic['1_CB'].flatten()
        d = resdic['1_CA'].flatten()
        omg = calc_dihedral(a, b, c, d)
        self.assertAlmostEqual(omg, label_funcs.compute_omg(0, 1, resdic))

    def test_compute_theta(self):
        resdic, _, _ = label_funcs.read_chain('data/small_example.pdb')
        a = resdic['0_N'].flatten()
        b = resdic['0_CA'].flatten()
        c = resdic['0_CB'].flatten()
        d = resdic['1_CB'].flatten()
        theta = calc_dihedral(a, b, c, d)
        self.assertAlmostEqual(theta, label_funcs.compute_theta(0, 1, resdic))

    def test_compute_phi(self):
        resdic, _, _ = label_funcs.read_chain('data/small_example.pdb')
        a = resdic['0_CA'].flatten()
        b = resdic['0_CB'].flatten()
        c = resdic['1_CB'].flatten()
        phi = calc_angle_b(a, b, c)
        self.assertAlmostEqual(phi, label_funcs.compute_phi(0, 1, resdic))

    def test_unit_vector_1(self):
         vec = np.array([10,0,0])
         unit = label_funcs.unit_vector(vec)
         self.assertTrue(np.allclose(unit, [1,0,0]))

    def test_unit_vector_2(self):
         vec = np.array([53.12, -4, 28.6])
         unit = label_funcs.unit_vector(vec)
         self.assertTrue(np.allclose(unit, [0.878563, -0.0661569, 0.473022]))


if __name__ == '__main__':
    unittest.main()
