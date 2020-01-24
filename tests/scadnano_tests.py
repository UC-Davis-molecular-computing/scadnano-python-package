import os
# import sys
# sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import math
import unittest
import re
import json
from typing import Iterable

import scadnano as sc
import scadnano.origami_rectangle as rect


# TODO: add tests for mutation methods on DNADesign

# TODO: mutator methods let me create this strand (which should be illegal); add a test for it
# {
#   "color": {"r": 51, "g": 51, "b": 51},
#   "substrands": [
#     {"helix": 2, "forward": false, "start": 40, "end": 48},
#     {"helix": 2, "forward": false, "start": 32, "end": 48, "deletions": [44]},
#     {"helix": 3, "forward": true, "start": 32, "end": 40}
#   ]
# }

def strand_matching(strands: Iterable[sc.Strand], helix: int, forward: bool, start: int, end: int):
    """
    Finds strand whose first bound substrand matches the given parameters.
    """
    return next(s for s in strands if
                s.first_bound_substrand().helix == helix and
                s.first_bound_substrand().forward == forward and
                s.first_bound_substrand().start == start and
                s.first_bound_substrand().end == end)


def remove_whitespace(sequence):
    sequence = re.sub(r'\s*', '', sequence)
    return sequence

class TestImportCadnanoV2(unittest.TestCase):
    """
    Tests the import feature to cadnano v2 (see misc/cadnano-format-specs/v2.txt).
    """
    io_folder = "cadnano_v2_import"
    def test_32_helix_rectangle(self):
        design = sc.DNADesign.from_cadnano_v2(directory=os.path.join('tests_outputs',self.io_folder), 
                           filename='test_32_helix_rectangle.json')
        design.write_scadnano_file(directory=os.path.join('tests_outputs',self.io_folder), 
                                   filename='test_32_helix_rectangle.dna')

    def test_helices_order(self):
        design = sc.DNADesign.from_cadnano_v2(directory=os.path.join('tests_outputs',self.io_folder), 
                           filename='test_helices_order.json')
        design.write_scadnano_file(directory=os.path.join('tests_outputs',self.io_folder), 
                                   filename='test_helices_order.dna')


class TestExportCadnanoV2(unittest.TestCase):
    """
    Tests the export feature to cadnano v2 (see misc/cadnano-format-specs/v2.txt).
    """
    output_folder = "cadnano_v2_export"
    def test_2_stape_2_helix_origami_extremely_simple(self):
        helices = [sc.Helix(max_offset=32), sc.Helix(max_offset=32)]
        scaf_part = sc.Substrand(helix=0, forward=True, start=0, end=32)
        scaf = sc.Strand(substrands=[scaf_part])
        design = sc.DNAOrigamiDesign(helices=helices, strands=[scaf], grid=sc.square, scaffold=scaf)
        design.write_scadnano_file(directory=os.path.join('tests_outputs',self.output_folder), 
                                   filename='test_2_stape_2_helix_origami_extremely_simple.dna')
        design.export_cadnano_v2(directory=os.path.join('tests_outputs',self.output_folder), 
                                 filename='test_2_stape_2_helix_origami_extremely_simple.json')

    def test_2_stape_2_helix_origami_extremely_simple_2(self):
        helices = [sc.Helix(max_offset=32), sc.Helix(max_offset=32)]
        scaf_part1 = sc.Substrand(helix=0, forward=True, start=0, end=32)
        scaf_part2 = sc.Substrand(helix=1, forward=False, start=0, end=32)
        scaf = sc.Strand(substrands=[scaf_part1,scaf_part2])
        design = sc.DNAOrigamiDesign(helices=helices, strands=[scaf], grid=sc.square, scaffold=scaf)
        design.write_scadnano_file(directory=os.path.join('tests_outputs',self.output_folder), 
                                   filename='test_2_stape_2_helix_origami_extremely_simple_2.dna')
        design.export_cadnano_v2(directory=os.path.join('tests_outputs',self.output_folder), 
                                 filename='test_2_stape_2_helix_origami_extremely_simple_2.json')


    def test_2_stape_2_helix_origami_deletions_insertions(self):
        # left staple
        stap_left_ss1 = sc.Substrand(helix=1, forward=True, start=0, end=16)
        stap_left_ss0 = sc.Substrand(helix=0, forward=False, start=0, end=16)
        stap_left = sc.Strand(substrands=[stap_left_ss1, stap_left_ss0])
    
        # right staple
        stap_right_ss0 = sc.Substrand(helix=0, forward=False, start=16, end=32)
        stap_right_ss1 = sc.Substrand(helix=1, forward=True, start=16, end=32)
        stap_right = sc.Strand(substrands=[stap_right_ss0, stap_right_ss1])
    
        # scaffold
        scaf_ss1_left = sc.Substrand(helix=1, forward=False, start=0, end=16)
        scaf_ss0 = sc.Substrand(helix=0, forward=True, start=0, end=32)
        #loopout = sc.Loopout(length=3) No loopout in cadnano
        scaf_ss1_right = sc.Substrand(helix=1, forward=False, start=16, end=32)
        scaf = sc.Strand(substrands=[scaf_ss1_left, scaf_ss0, scaf_ss1_right])
    
        # whole design
        design = sc.DNAOrigamiDesign(strands=[scaf, stap_left, stap_right], grid=sc.square, scaffold=scaf)
    
        # deletions and insertions added to design so they can be added to both strands on a helix
        design.add_deletion(helix=0, offset=11)
        design.add_deletion(helix=0, offset=12)
        design.add_deletion(helix=0, offset=24)
        design.add_deletion(helix=1, offset=12)
        design.add_deletion(helix=1, offset=24)
    
        design.add_insertion(helix=0, offset=6, length=1)
        design.add_insertion(helix=0, offset=18, length=2)
        design.add_insertion(helix=1, offset=6, length=3)
        design.add_insertion(helix=1, offset=18, length=4)
    
        # also assigns complement to strands other than scaf bound to it
        design.assign_dna(scaf, 'AACGT' * 18)
        design.write_scadnano_file(directory=os.path.join('tests_outputs',self.output_folder), 
                                   filename='test_2_stape_2_helix_origami_deletions_insertions.dna')
        design.export_cadnano_v2(directory=os.path.join('tests_outputs',self.output_folder), 
                                 filename='test_2_stape_2_helix_origami_deletions_insertions.json')
    
    def test_6_helix_origami_rectangle(self):
        design = rect.create(num_helices=6, num_cols=10, nick_pattern=rect.staggered, twist_correction_deletion_spacing=3)
        design.write_scadnano_file(directory=os.path.join('tests_outputs',self.output_folder), 
                                   filename='test_6_helix_origami_rectangle.dna')
        design.export_cadnano_v2(directory=os.path.join('tests_outputs',self.output_folder), 
                                 filename='test_6_helix_origami_rectangle.json')

    def test_16_helix_origami_rectangle_no_twist(self):
        design = rect.create(num_helices=16, num_cols=26, assign_seq=True, twist_correction_deletion_spacing=3)
        design.write_scadnano_file(directory=os.path.join('tests_outputs',self.output_folder), 
                                   filename='test_16_helix_origami_rectangle_no_twist.dna')
        design.export_cadnano_v2(directory=os.path.join('tests_outputs',self.output_folder), 
                                 filename='test_16_helix_origami_rectangle_no_twist.json')

    def test_bad_cases(self):
        """ We do not handle Loopouts and design where the parity of the helix
        does not correspond to the direction.
        """
        
        # Bad case one: parity issue in design (see cadnano v2 format spec, v2.txt)
        helices = [sc.Helix(max_offset=32), sc.Helix(max_offset=32)]
        scaf_part = sc.Substrand(helix=1, forward=True, start=0, end=32)
        scaf = sc.Strand(substrands=[scaf_part])
        design = sc.DNAOrigamiDesign(helices=helices, strands=[scaf], grid=sc.square, scaffold=scaf)
        
        with self.assertRaises(ValueError) as context:
            design.export_cadnano_v2(directory=os.path.join('tests_outputs',self.output_folder), 
                                 filename='test_parity_issue.json')
        self.assertTrue('forward' in context.exception.args[0])

        # Bad case two: Loopouts
        helices = [sc.Helix(max_offset=48), sc.Helix(max_offset=48)]

        # left staple
        stap_left_ss1 = sc.Substrand(helix=1, forward=True, start=8, end=24)
        stap_left_ss0 = sc.Substrand(helix=0, forward=False, start=8, end=24)
        stap_left = sc.Strand(substrands=[stap_left_ss1, stap_left_ss0])

        # right staple
        stap_right_ss0 = sc.Substrand(helix=0, forward=False, start=24, end=40)
        stap_right_ss1 = sc.Substrand(helix=1, forward=True, start=24, end=40)
        stap_right = sc.Strand(substrands=[stap_right_ss0, stap_right_ss1])

        # scaffold
        scaf_ss1_left = sc.Substrand(helix=1, forward=False, start=8, end=24)
        scaf_ss0 = sc.Substrand(helix=0, forward=True, start=8, end=40)
        loopout = sc.Loopout(length=3)
        scaf_ss1_right = sc.Substrand(helix=1, forward=False, start=24, end=40)
        scaf = sc.Strand(substrands=[scaf_ss1_left, scaf_ss0, loopout, scaf_ss1_right])

        # whole design
        design = sc.DNAOrigamiDesign(helices=helices, strands=[scaf, stap_left, stap_right], grid=sc.square,
                                    scaffold=scaf)

        # deletions and insertions added to design are added to both strands on a helix
        design.add_deletion(helix=1, offset=20)
        design.add_insertion(helix=0, offset=14, length=1)
        design.add_insertion(helix=0, offset=26, length=2)

        with self.assertRaises(ValueError) as context:
            design.export_cadnano_v2(directory=os.path.join('tests_outputs',self.output_folder), 
                                 filename='test_loopout_issue.json')
        self.assertTrue('Loopouts' in context.exception.args[0])



class TestDesignFromJson(unittest.TestCase):
    """
    Tests reading a design from a dict derived from JSON.
    """

    def test_from_json__three_strands(self):
        """
        0       8       16
        |       |       |
    0   +--X-----------+
       /<--X---++------]\
       |       ||       loopout(3)
       \[---2--++------>/
    1   +---2--]<------+
        """
        st_l = sc.Strand([
            sc.Substrand(1, True, 0, 8, insertions=[(4, 2)]),
            sc.Substrand(0, False, 0, 8, deletions=[3]),
        ])
        st_r = sc.Strand([
            sc.Substrand(0, False, 8, 16),
            sc.Substrand(1, True, 8, 16),
        ])
        scaf = sc.Strand([
            sc.Substrand(1, False, 0, 8, insertions=[(4, 2)]),
            sc.Substrand(0, True, 0, 16, deletions=[3]),
            sc.Loopout(3),
            sc.Substrand(1, False, 8, 16, deletions=[]),
        ])
        design_pre_json = sc.DNADesign(strands=[st_l, st_r, scaf], grid=sc.square)
        design_pre_json.assign_dna(scaf, 'A' * 36)

        json_str = design_pre_json.to_json()
        json_map = json.loads(json_str)
        design = sc.DNADesign.from_json(json_map)

        self.assertEquals(sc.Grid.square, design.grid)

        self.assertEquals(2, len(design.helices))
        helix0 = design.helices[0]
        helix1 = design.helices[1]
        self.assertEquals(0, helix0.idx)
        self.assertEquals(0, helix0.min_offset)
        self.assertEquals(16, helix0.max_offset)
        self.assertEquals((0, 0, 0), helix0.grid_position)
        self.assertEquals(1, helix1.idx)
        self.assertEquals(0, helix1.min_offset)
        self.assertEquals(16, helix1.max_offset)
        self.assertEquals((0, 1, 0), helix1.grid_position)

        self.assertEquals(3, len(design.strands))
        st_l = design.strands[0]
        st_r = design.strands[1]
        scaf = design.strands[2]

        self.assertEquals(2, len(st_l.substrands))
        self.assertEquals(2, len(st_r.substrands))
        self.assertEquals(4, len(scaf.substrands))

        self.assertEquals('A' * 36, scaf.dna_sequence)
        self.assertEquals('T' * 17, st_l.dna_sequence)
        self.assertEquals('T' * 16, st_r.dna_sequence)

        st_l_ss0 = st_l.substrands[0]
        st_l_ss1 = st_l.substrands[1]
        st_r_ss0 = st_r.substrands[0]
        st_r_ss1 = st_r.substrands[1]
        scaf_ss0 = scaf.substrands[0]
        scaf_ss1 = scaf.substrands[1]
        scaf_loop = scaf.substrands[2]
        scaf_ss2 = scaf.substrands[3]

        self.assertEquals(3, scaf_loop.length)

        self.assertEquals(1, st_l_ss0.helix)
        self.assertEquals(0, st_l_ss1.helix)
        self.assertEquals(0, st_r_ss0.helix)
        self.assertEquals(1, st_r_ss1.helix)
        self.assertEquals(1, scaf_ss0.helix)
        self.assertEquals(0, scaf_ss1.helix)
        self.assertEquals(1, scaf_ss2.helix)

        self.assertEquals(True, st_l_ss0.forward)
        self.assertEquals(False, st_l_ss1.forward)
        self.assertEquals(False, st_r_ss0.forward)
        self.assertEquals(True, st_r_ss1.forward)
        self.assertEquals(False, scaf_ss0.forward)
        self.assertEquals(True, scaf_ss1.forward)
        self.assertEquals(False, scaf_ss2.forward)

        self.assertEquals(0, st_l_ss0.start)
        self.assertEquals(8, st_l_ss0.end)
        self.assertEquals(0, st_l_ss1.start)
        self.assertEquals(8, st_l_ss1.end)
        self.assertEquals(8, st_r_ss0.start)
        self.assertEquals(16, st_r_ss0.end)
        self.assertEquals(8, st_r_ss1.start)
        self.assertEquals(16, st_r_ss1.end)
        self.assertEquals(0, scaf_ss0.start)
        self.assertEquals(8, scaf_ss0.end)
        self.assertEquals(0, scaf_ss1.start)
        self.assertEquals(16, scaf_ss1.end)
        self.assertEquals(8, scaf_ss2.start)
        self.assertEquals(16, scaf_ss2.end)

        self.assertListEqual([(4, 2)], st_l_ss0.insertions)
        self.assertListEqual([], st_l_ss0.deletions)
        self.assertListEqual([], st_l_ss1.insertions)
        self.assertListEqual([3], st_l_ss1.deletions)
        self.assertListEqual([], st_r_ss0.insertions)
        self.assertListEqual([], st_r_ss0.deletions)
        self.assertListEqual([], st_r_ss1.insertions)
        self.assertListEqual([], st_r_ss1.deletions)
        self.assertListEqual([(4, 2)], scaf_ss0.insertions)
        self.assertListEqual([], scaf_ss0.deletions)
        self.assertListEqual([], scaf_ss1.insertions)
        self.assertListEqual([3], scaf_ss1.deletions)
        self.assertListEqual([], scaf_ss2.insertions)
        self.assertListEqual([], scaf_ss2.deletions)


class TestStrandReversePolarity(unittest.TestCase):
    """
    Tests reversing polarity of all strands in the system.
    """

    def test_reverse_all__one_strand(self):
        """
        before
        0       8
        |       |
    0   [---X-->

        after
    0   <---X--]
        """
        design = sc.DNADesign(
            strands=[sc.Strand([sc.Substrand(0, True, 0, 8, deletions=[4])])],
            grid=sc.square)
        design.reverse_all()

        self.assertEqual(1, len(design.strands))
        strand = design.strands[0]
        self.assertEqual(1, len(strand.substrands))
        self.assertEqual(7, strand.offset_5p())
        self.assertEqual(0, strand.offset_3p())
        ss = strand.substrands[0]
        self.assertEqual(0, ss.helix)
        self.assertEqual(False, ss.forward)
        self.assertEqual(0, ss.start)
        self.assertEqual(8, ss.end)

    def test_reverse_all__three_strands(self):
        """
        before
        0       8       16
        |       |       |

    0   +--X-----------+
       /<--X---++------]\
       |       ||       |
       \[---2--++------>/
    1   +---2--]<------+

        after
    0   +--X-----------+
       /[--X---++------>\
       |       ||       |
       \<---2--++------]/
    1   +---2-->[------+
        """
        design = sc.DNADesign(
            strands=[
                sc.Strand([
                    sc.Substrand(1, True, 0, 8, insertions=[(4, 2)]),
                    sc.Substrand(0, False, 0, 8, deletions=[3]),
                ]),
                sc.Strand([
                    sc.Substrand(0, False, 8, 16),
                    sc.Substrand(1, True, 8, 16),
                ]),
                sc.Strand([
                    sc.Substrand(1, False, 0, 8, insertions=[(4, 2)]),
                    sc.Substrand(0, True, 0, 16, deletions=[3]),
                    sc.Substrand(1, False, 8, 16, deletions=[]),
                ]),
            ],
            grid=sc.square)
        design.reverse_all()

        self.assertEqual(3, len(design.strands))
        stapL = design.strands[0]
        stapR = design.strands[1]
        scaf = design.strands[2]

        self.assertEqual(0, stapL.offset_5p())
        self.assertEqual(0, stapL.offset_3p())
        self.assertEqual(15, stapR.offset_5p())
        self.assertEqual(15, stapR.offset_3p())
        self.assertEqual(8, scaf.offset_5p())
        self.assertEqual(7, scaf.offset_3p())

        self.assertEqual(2, len(stapL.substrands))
        stapL_ss0 = stapL.substrands[0]
        stapL_ss1 = stapL.substrands[1]

        self.assertEqual(0, stapL_ss0.helix)
        self.assertEqual(True, stapL_ss0.forward)
        self.assertEqual(0, stapL_ss0.start)
        self.assertEqual(8, stapL_ss0.end)
        self.assertEqual([3], stapL_ss0.deletions)
        self.assertEqual([], stapL_ss0.insertions)

        self.assertEqual(1, stapL_ss1.helix)
        self.assertEqual(False, stapL_ss1.forward)
        self.assertEqual(0, stapL_ss1.start)
        self.assertEqual(8, stapL_ss1.end)
        self.assertEqual([], stapL_ss1.deletions)
        self.assertEqual([(4, 2)], stapL_ss1.insertions)

        self.assertEqual(2, len(stapR.substrands))
        stapR_ss0 = stapR.substrands[0]
        stapR_ss1 = stapR.substrands[1]

        self.assertEqual(1, stapR_ss0.helix)
        self.assertEqual(False, stapR_ss0.forward)
        self.assertEqual(8, stapR_ss0.start)
        self.assertEqual(16, stapR_ss0.end)
        self.assertEqual([], stapR_ss0.deletions)
        self.assertEqual([], stapR_ss0.insertions)

        self.assertEqual(0, stapR_ss1.helix)
        self.assertEqual(True, stapR_ss1.forward)
        self.assertEqual(8, stapR_ss1.start)
        self.assertEqual(16, stapR_ss1.end)
        self.assertEqual([], stapR_ss1.deletions)
        self.assertEqual([], stapR_ss1.insertions)

        self.assertEqual(3, len(scaf.substrands))
        scaf_ss0 = scaf.substrands[0]
        scaf_ss1 = scaf.substrands[1]
        scaf_ss2 = scaf.substrands[2]

        self.assertEqual(1, scaf_ss0.helix)
        self.assertEqual(True, scaf_ss0.forward)
        self.assertEqual(8, scaf_ss0.start)
        self.assertEqual(16, scaf_ss0.end)
        self.assertEqual([], scaf_ss0.deletions)
        self.assertEqual([], scaf_ss0.insertions)

        self.assertEqual(0, scaf_ss1.helix)
        self.assertEqual(False, scaf_ss1.forward)
        self.assertEqual(0, scaf_ss1.start)
        self.assertEqual(16, scaf_ss1.end)
        self.assertEqual([3], scaf_ss1.deletions)
        self.assertEqual([], scaf_ss1.insertions)

        self.assertEqual(1, scaf_ss2.helix)
        self.assertEqual(True, scaf_ss2.forward)
        self.assertEqual(0, scaf_ss2.start)
        self.assertEqual(8, scaf_ss2.end)
        self.assertEqual([], scaf_ss2.deletions)
        self.assertEqual([(4, 2)], scaf_ss2.insertions)


class TestInlineInsDel(unittest.TestCase):
    """
    Tests inlining of insertions/deletions.
    """

    def test_inline_deletions_insertions__one_deletion(self):
        """
        before
        0   4   8       16      24
        |       |       |       |
    0   [---X-->

        after
        0   4  7       15      23
        |      |       |       |
    0   [----->
        """
        design = sc.DNADesign(
            helices=[sc.Helix(max_offset=24, major_tick_distance=8)],
            strands=[sc.Strand([sc.Substrand(0, True, 0, 8, deletions=[4])])],
            grid=sc.square)
        design.inline_deletions_insertions()
        self.helix0_strand0_inlined_test(design, max_offset=23, major_ticks=[0, 7, 15, 23], start=0, end=7)

    def helix0_strand0_inlined_test(self, design, max_offset, major_ticks, start, end):
        self.assertEqual(1, len(design.helices))
        self.assertEqual(1, len(design.strands))
        helix = design.helices[0]
        strand = design.strands[0]
        self.assertEqual(max_offset, helix.max_offset)
        self.assertEqual(major_ticks, helix.major_ticks)
        self.assertEqual(start, strand.substrands[0].start)
        self.assertEqual(end, strand.substrands[0].end)
        self.assertListEqual([], strand.substrands[0].deletions)
        self.assertListEqual([], strand.substrands[0].insertions)

    def test_inline_deletions_insertions__two_deletions(self):
        """
        before
        0 2 4   8       16      24
        |       |       |       |
    0   [-X-X-->

        after
        0     6       14      22
        |     |       |       |
    0   [---->
        """
        design = sc.DNADesign(
            helices=[sc.Helix(max_offset=24, major_tick_distance=8)],
            strands=[sc.Strand([sc.Substrand(0, True, 0, 8, deletions=[2, 4])])],
            grid=sc.square)
        design.inline_deletions_insertions()
        self.helix0_strand0_inlined_test(design, max_offset=22, major_ticks=[0, 6, 14, 22], start=0, end=6)

    def test_inline_deletions_insertions__one_insertion(self):
        """
        before
        0   4   8       16      24
        |       |       |       |
    0   [---1-->

        after
        0        9       17      25
        |        |       |       |
    0   [------->
        """
        design = sc.DNADesign(
            helices=[sc.Helix(max_offset=24, major_tick_distance=8)],
            strands=[sc.Strand([sc.Substrand(0, True, 0, 8, insertions=[(4, 1)])])],
            grid=sc.square)
        design.inline_deletions_insertions()
        self.helix0_strand0_inlined_test(design, max_offset=25, major_ticks=[0, 9, 17, 25], start=0, end=9)

    def test_inline_deletions_insertions__two_insertions(self):
        """
        before
        0 2 4   8       16      24
        |       |       |       |
    0   [-3-1-->

        after
        0           12      20      28
        |           |       |       |
    0   [---------->
        """
        design = sc.DNADesign(
            helices=[sc.Helix(max_offset=24, major_tick_distance=8)],
            strands=[sc.Strand([sc.Substrand(0, True, 0, 8, insertions=[(2, 3), (4, 1)])])],
            grid=sc.square)
        design.inline_deletions_insertions()
        self.helix0_strand0_inlined_test(design, max_offset=28, major_ticks=[0, 12, 20, 28], start=0, end=12)

    def test_inline_deletions_insertions__one_deletion_one_insertion(self):
        """
        before
        0 2 4   8       16      24
        |       |       |       |
    0   [-3-X-->

        after
        0         10      18      26
        |         |       |       |
    0   [-------->
        """
        design = sc.DNADesign(
            helices=[sc.Helix(max_offset=24, major_tick_distance=8)],
            strands=[sc.Strand([sc.Substrand(0, True, 0, 8, deletions=[4], insertions=[(2, 3)])])],
            grid=sc.square)
        design.inline_deletions_insertions()
        self.helix0_strand0_inlined_test(design, max_offset=26, major_ticks=[0, 10, 18, 26], start=0, end=10)

    def test_inline_deletions_insertions__one_deletion_right_of_major_tick(self):
        """
        before
        0       89      16      24
        |       |       |       |
    0   [--------X->

        after
        0       8      15      23
        |       |      |       |
    0   [--------->
        """
        design = sc.DNADesign(
            helices=[sc.Helix(max_offset=24, major_tick_distance=8)],
            strands=[sc.Strand([sc.Substrand(0, True, 0, 12, deletions=[9])])],
            grid=sc.square)
        design.inline_deletions_insertions()
        self.helix0_strand0_inlined_test(design, max_offset=23, major_ticks=[0, 8, 15, 23], start=0, end=11)

    def test_inline_deletions_insertions__one_deletion_on_major_tick(self):
        """
        | is major tick, and . is minor tick
        before
         0               8               16              24
        | . . . . . . . | . . . . . . . | . . . . . . . |
         [ - - - - - - - X - - >

        after
         0               8             15              23
        | . . . . . . . | . . . . . . | . . . . . . . |
         [ - - - - - - - - - >
        """
        design = sc.DNADesign(
            helices=[sc.Helix(max_offset=24, major_tick_distance=8)],
            strands=[sc.Strand([sc.Substrand(0, True, 0, 12, deletions=[8])])],
            grid=sc.square)
        design.inline_deletions_insertions()
        self.helix0_strand0_inlined_test(design, max_offset=23, major_ticks=[0, 8, 15, 23], start=0, end=11)

    def test_inline_deletions_insertions__one_deletion_left_of_major_tick(self):
        """
        before
        0      78       16      24
        |       |       |       |
    0   [------X--->

        after
        0       8      15      23
        |       |      |       |
    0   [--------->
        """
        design = sc.DNADesign(
            helices=[sc.Helix(max_offset=24, major_tick_distance=8)],
            strands=[sc.Strand([sc.Substrand(0, True, 0, 12, deletions=[7])])],
            grid=sc.square)
        design.inline_deletions_insertions()
        self.helix0_strand0_inlined_test(design, max_offset=23, major_ticks=[0, 7, 15, 23], start=0, end=11)

    def test_inline_deletions_insertions__one_insertion_right_of_major_tick(self):
        """
        before
        0       89      16      24
        |       |       |       |
    0   [--------1->

        after
        0       8        17      25
        |       |        |       |
    0   [----------->
        """
        design = sc.DNADesign(
            helices=[sc.Helix(max_offset=24, major_tick_distance=8)],
            strands=[sc.Strand([sc.Substrand(0, True, 0, 12, insertions=[(9, 1)])])],
            grid=sc.square)
        design.inline_deletions_insertions()
        self.helix0_strand0_inlined_test(design, max_offset=25, major_ticks=[0, 8, 17, 25], start=0, end=13)

    def test_inline_deletions_insertions__one_insertion_on_major_tick(self):
        """
        before
        0       8       16      24
        |       |       |       |
    0   [-------1-->

        after
        0       8        17      25
        |       |        |       |
    0   [----------->
        """
        design = sc.DNADesign(
            helices=[sc.Helix(max_offset=24, major_tick_distance=8)],
            strands=[sc.Strand([sc.Substrand(0, True, 0, 12, insertions=[(8, 1)])])],
            grid=sc.square)
        design.inline_deletions_insertions()
        self.helix0_strand0_inlined_test(design, max_offset=25, major_ticks=[0, 8, 17, 25], start=0, end=13)

    def test_inline_deletions_insertions__one_insertion_left_of_major_tick(self):
        """
        before
        0      78       16      24
        |       |       |       |
    0   [------1--->

        after
        0        9       17      25
        |        |       |       |
    0   [----------->
        """
        design = sc.DNADesign(
            helices=[sc.Helix(max_offset=24, major_tick_distance=8)],
            strands=[sc.Strand([sc.Substrand(0, True, 0, 12, insertions=[(7, 1)])])],
            grid=sc.square)
        design.inline_deletions_insertions()
        self.helix0_strand0_inlined_test(design, max_offset=25, major_ticks=[0, 9, 17, 25], start=0, end=13)

    def test_inline_deletions_insertions__deletions_insertions_in_multiple_domains(self):
        """
        before
        0    5  8  11   16 19   24
        |       |       |       |
    0   [----2-----1-------X---->

        after
        0         10      19      26
        |         |       |       |
    0   [------------------------->
        """
        design = sc.DNADesign(
            helices=[sc.Helix(max_offset=24, major_tick_distance=8)],
            strands=[sc.Strand([sc.Substrand(0, True, 0, 24, deletions=[19], insertions=[(5, 2), (11, 1)])])],
            grid=sc.square)
        design.inline_deletions_insertions()
        self.helix0_strand0_inlined_test(design, max_offset=26, major_ticks=[0, 10, 19, 26], start=0, end=26)

    def test_inline_deletions_insertions__deletions_insertions_in_multiple_domains_two_strands(self):
        """
        | is major tick, . is minor tick
        before
         0   2     5     8   10          16    19        24
        | . . . . . . . | . . . . . . . | . . . . . . . |
         [ - X - - 2 - - - - 1 - - > [ - - - - X - - - >

        after
         0   2     5       9             16  18            25
        | . . . . . . . . | . . . . . . . . | . . . . . . |
         [ - - - - - - - - - - - - - - > [ - - - - - - - >
        """
        design = sc.DNADesign(
            helices=[sc.Helix(max_offset=24, major_tick_distance=8)],
            strands=[
                sc.Strand([sc.Substrand(0, True, 0, 14, deletions=[2], insertions=[(5, 2), (10, 1)])]),
                sc.Strand([sc.Substrand(0, True, 14, 24, deletions=[19])]),
            ],
            grid=sc.square)
        design.inline_deletions_insertions()
        self.assertEqual(1, len(design.helices))
        self.assertEqual(2, len(design.strands))
        helix = design.helices[0]
        strand0 = design.strands[0]
        strand1 = design.strands[1]
        self.assertEqual(25, helix.max_offset)
        self.assertEqual([0, 9, 18, 25], helix.major_ticks)
        self.assertEqual(0, strand0.substrands[0].start)
        self.assertEqual(16, strand0.substrands[0].end)
        self.assertEqual(16, strand1.substrands[0].start)
        self.assertEqual(25, strand1.substrands[0].end)
        self.assertListEqual([], strand0.substrands[0].deletions)
        self.assertListEqual([], strand0.substrands[0].insertions)
        self.assertListEqual([], strand1.substrands[0].deletions)
        self.assertListEqual([], strand1.substrands[0].insertions)


class TestNickAndCrossover(unittest.TestCase):
    """
    Tests add_nick() and add_*_crossover() methods on DNADesign as an easier way of specifying an origami.
    """

    r"""
    small_design:
    0        8        16
    ACGTACGA AACCGGTA
0   [------- ------->
    <------- -------]
    TGCATGCT TTGGCCAT

    AAACCCGG TTTGGGCC
1   [------- ------->
    <------- -------]
    TTTGGGCC AAACCCGG


    design:
    0        8        16       24       32       40       48       56       64       72       80       88       96
0   [------- -------- -------- -------- -------- -------- -------- -------- -------- -------- -------- ------->
    <------- -------- -------- -------- -------- -------- -------- -------- -------- -------- -------- -------]

1   [------- -------- -------- -------- -------- -------- -------- -------- -------- -------- -------- ------->
    <------- -------- -------- -------- -------- -------- -------- -------- -------- -------- -------- -------]

2   [------- -------- -------- -------- -------- -------- -------- -------- -------- -------- -------- ------->
    <------- -------- -------- -------- -------- -------- -------- -------- -------- -------- -------- -------]

3   [------- -------- -------- -------- -------- -------- -------- -------- -------- -------- -------- ------->
    <------- -------- -------- -------- -------- -------- -------- -------- -------- -------- -------- -------]

4   [------- -------- -------- -------- -------- -------- -------- -------- -------- -------- -------- ------->
    <------- -------- -------- -------- -------- -------- -------- -------- -------- -------- -------- -------]

5   [------- -------- -------- -------- -------- -------- -------- -------- -------- -------- -------- ------->
    <------- -------- -------- -------- -------- -------- -------- -------- -------- -------- -------- -------]
    """

    def setUp(self):
        strands_small_design = [
            sc.Strand([sc.Substrand(0, True, 0, 16)]),
            sc.Strand([sc.Substrand(0, False, 0, 16)]),
            sc.Strand([sc.Substrand(1, True, 0, 16)]),
            sc.Strand([sc.Substrand(1, False, 0, 16)]),
        ]
        self.small_design = sc.DNADesign(strands=strands_small_design, grid=sc.square)
        self.small_design.assign_dna(strands_small_design[0], "ACGTACGA AACCGGTA")
        self.small_design.assign_dna(strands_small_design[2], "AAACCCGG TTTGGGCC")

        self.max_offset: int = 8 * 12
        scafs = []
        staps = []
        for helix in range(6):
            scaf_ss = sc.Substrand(helix, helix % 2 == 0, 0, self.max_offset)
            stap_ss = sc.Substrand(helix, helix % 2 == 1, 0, self.max_offset)
            scaf = sc.Strand([scaf_ss])
            stap = sc.Strand([stap_ss])
            scafs.append(scaf)
            staps.append(stap)
        self.design: sc.DNADesign = sc.DNADesign(strands=scafs + staps, grid=sc.square)

    def test_add_nick__twice_on_same_substrand(self):
        """
        before
        0        8        16       24
    0   [------- -------- ------->

        after
        0        8        16       24
    0   [------> [------> [------>
        """
        design = sc.DNADesign(strands=[
            sc.Strand([sc.Substrand(0, True, 0, 24)]),
        ], grid=sc.square)
        design.add_nick(helix=0, offset=8, forward=True)
        design.add_nick(helix=0, offset=16, forward=True)
        self.assertEqual(3, len(design.strands))
        self.assertIn(sc.Strand([sc.Substrand(0, True, 0, 8)]), design.strands)
        self.assertIn(sc.Strand([sc.Substrand(0, True, 8, 16)]), design.strands)
        self.assertIn(sc.Strand([sc.Substrand(0, True, 16, 24)]), design.strands)

    def test_add_nick__small_design_no_nicks_added(self):
        """
        0        8        16
        ACGTACGA AACCGGTA
    0   [------- ------->
        <------- -------]
        TGCATGCT TTGGCCAT

        AAACCCGG TTTGGGCC
    1   [------- ------->
        <------- -------]
        TTTGGGCC AAACCCGG
        """
        self.assertEqual(4, len(self.small_design.strands))
        self.assertIn(sc.Strand([sc.Substrand(0, True, 0, 16)]), self.small_design.strands)
        self.assertIn(sc.Strand([sc.Substrand(0, False, 0, 16)]), self.small_design.strands)
        self.assertIn(sc.Strand([sc.Substrand(1, True, 0, 16)]), self.small_design.strands)
        self.assertIn(sc.Strand([sc.Substrand(1, False, 0, 16)]), self.small_design.strands)
        # DNA
        strand = strand_matching(self.small_design.strands, 0, True, 0, 16)
        self.assertEqual(remove_whitespace('ACGTACGA AACCGGTA'), strand.dna_sequence)
        strand = strand_matching(self.small_design.strands, 0, False, 0, 16)
        self.assertEqual(remove_whitespace('TACCGGTT TCGTACGT'), strand.dna_sequence)
        strand = strand_matching(self.small_design.strands, 1, True, 0, 16)
        self.assertEqual(remove_whitespace('AAACCCGG TTTGGGCC'), strand.dna_sequence)
        strand = strand_matching(self.small_design.strands, 1, False, 0, 16)
        self.assertEqual(remove_whitespace('GGCCCAAA CCGGGTTT'), strand.dna_sequence)

    def test_add_nick__small_design_H0_forward(self):
        """
        0        8        16
        ACGTACGA AACCGGTA
    0   [------> [------>
        <------- -------]
        TGCATGCT TTGGCCAT

        AAACCCGG TTTGGGCC
    1   [------- ------->
        <------- -------]
        TTTGGGCC AAACCCGG
        """
        self.small_design.add_nick(helix=0, offset=8, forward=True)
        self.assertEqual(5, len(self.small_design.strands))
        # two new Strands
        self.assertIn(sc.Strand([sc.Substrand(0, True, 0, 8)]), self.small_design.strands)
        self.assertIn(sc.Strand([sc.Substrand(0, True, 8, 16)]), self.small_design.strands)
        # existing Strands
        self.assertIn(sc.Strand([sc.Substrand(0, False, 0, 16)]), self.small_design.strands)
        self.assertIn(sc.Strand([sc.Substrand(1, True, 0, 16)]), self.small_design.strands)
        self.assertIn(sc.Strand([sc.Substrand(1, False, 0, 16)]), self.small_design.strands)
        # DNA
        strand = strand_matching(self.small_design.strands, 0, True, 0, 8)
        self.assertEqual(remove_whitespace('ACGTACGA'), strand.dna_sequence)
        strand = strand_matching(self.small_design.strands, 0, True, 8, 16)
        self.assertEqual(remove_whitespace('AACCGGTA'), strand.dna_sequence)
        strand = strand_matching(self.small_design.strands, 0, False, 0, 16)
        self.assertEqual(remove_whitespace('TACCGGTT TCGTACGT'), strand.dna_sequence)
        strand = strand_matching(self.small_design.strands, 1, True, 0, 16)
        self.assertEqual(remove_whitespace('AAACCCGG TTTGGGCC'), strand.dna_sequence)
        strand = strand_matching(self.small_design.strands, 1, False, 0, 16)
        self.assertEqual(remove_whitespace('GGCCCAAA CCGGGTTT'), strand.dna_sequence)

    def test_add_nick__small_design_H0_reverse(self):
        """
        0        8        16
        ACGTACGA AACCGGTA
    0   [------- ------->
        <------] <------]
        TGCATGCT TTGGCCAT

        AAACCCGG TTTGGGCC
    1   [------- ------->
        <------- -------]
        TTTGGGCC AAACCCGG
        """
        self.small_design.add_nick(helix=0, offset=8, forward=False)
        self.assertEqual(5, len(self.small_design.strands))
        # two new Strands
        self.assertIn(sc.Strand([sc.Substrand(0, False, 0, 8)]), self.small_design.strands)
        self.assertIn(sc.Strand([sc.Substrand(0, False, 8, 16)]), self.small_design.strands)
        # existing Strands
        self.assertIn(sc.Strand([sc.Substrand(0, True, 0, 16)]), self.small_design.strands)
        self.assertIn(sc.Strand([sc.Substrand(1, True, 0, 16)]), self.small_design.strands)
        self.assertIn(sc.Strand([sc.Substrand(1, False, 0, 16)]), self.small_design.strands)
        # DNA
        strand = strand_matching(self.small_design.strands, 0, True, 0, 16)
        self.assertEqual(remove_whitespace('ACGTACGA AACCGGTA'), strand.dna_sequence)
        strand = strand_matching(self.small_design.strands, 0, False, 8, 16)
        self.assertEqual(remove_whitespace('TACCGGTT'), strand.dna_sequence)
        strand = strand_matching(self.small_design.strands, 0, False, 0, 8)
        self.assertEqual(remove_whitespace('TCGTACGT'), strand.dna_sequence)
        strand = strand_matching(self.small_design.strands, 1, True, 0, 16)
        self.assertEqual(remove_whitespace('AAACCCGG TTTGGGCC'), strand.dna_sequence)
        strand = strand_matching(self.small_design.strands, 1, False, 0, 16)
        self.assertEqual(remove_whitespace('GGCCCAAA CCGGGTTT'), strand.dna_sequence)

    def test_add_nick__small_design_H1_forward(self):
        """
        0        8        16
        ACGTACGA AACCGGTA
    0   [------- ------->
        <------- -------]
        TGCATGCT TTGGCCAT

        AAACCCGG TTTGGGCC
    1   [------> [------>
        <------- -------]
        TTTGGGCC AAACCCGG
        """
        self.small_design.add_nick(helix=1, offset=8, forward=True)
        self.assertEqual(5, len(self.small_design.strands))
        # two new Strands
        self.assertIn(sc.Strand([sc.Substrand(1, True, 0, 8)]), self.small_design.strands)
        self.assertIn(sc.Strand([sc.Substrand(1, True, 8, 16)]), self.small_design.strands)
        # existing Strands
        self.assertIn(sc.Strand([sc.Substrand(1, False, 0, 16)]), self.small_design.strands)
        self.assertIn(sc.Strand([sc.Substrand(0, True, 0, 16)]), self.small_design.strands)
        self.assertIn(sc.Strand([sc.Substrand(0, False, 0, 16)]), self.small_design.strands)
        # DNA
        strand = strand_matching(self.small_design.strands, 0, True, 0, 16)
        self.assertEqual(remove_whitespace('ACGTACGA AACCGGTA'), strand.dna_sequence)
        strand = strand_matching(self.small_design.strands, 0, False, 0, 16)
        self.assertEqual(remove_whitespace('TACCGGTT TCGTACGT'), strand.dna_sequence)
        strand = strand_matching(self.small_design.strands, 1, True, 0, 8)
        self.assertEqual(remove_whitespace('AAACCCGG'), strand.dna_sequence)
        strand = strand_matching(self.small_design.strands, 1, True, 8, 16)
        self.assertEqual(remove_whitespace('TTTGGGCC'), strand.dna_sequence)
        strand = strand_matching(self.small_design.strands, 1, False, 0, 16)
        self.assertEqual(remove_whitespace('GGCCCAAA CCGGGTTT'), strand.dna_sequence)

    def test_add_nick__small_design_H1_reverse(self):
        """
        0        8        16
        ACGTACGA AACCGGTA
    0   [------- ------->
        <------- -------]
        TGCATGCT TTGGCCAT

        AAACCCGG TTTGGGCC
    1   [------- ------->
        <------] <------]
        TTTGGGCC AAACCCGG
        """
        self.small_design.add_nick(helix=1, offset=8, forward=False)
        self.assertEqual(5, len(self.small_design.strands))
        # two new Strands
        self.assertIn(sc.Strand([sc.Substrand(1, False, 0, 8)]), self.small_design.strands)
        self.assertIn(sc.Strand([sc.Substrand(1, False, 8, 16)]), self.small_design.strands)
        # existing Strands
        self.assertIn(sc.Strand([sc.Substrand(1, True, 0, 16)]), self.small_design.strands)
        self.assertIn(sc.Strand([sc.Substrand(0, True, 0, 16)]), self.small_design.strands)
        self.assertIn(sc.Strand([sc.Substrand(0, False, 0, 16)]), self.small_design.strands)
        # DNA
        strand = strand_matching(self.small_design.strands, 0, True, 0, 16)
        self.assertEqual(remove_whitespace('ACGTACGA AACCGGTA'), strand.dna_sequence)
        strand = strand_matching(self.small_design.strands, 0, False, 0, 16)
        self.assertEqual(remove_whitespace('TACCGGTT TCGTACGT'), strand.dna_sequence)
        strand = strand_matching(self.small_design.strands, 1, True, 0, 16)
        self.assertEqual(remove_whitespace('AAACCCGG TTTGGGCC'), strand.dna_sequence)
        strand = strand_matching(self.small_design.strands, 1, False, 0, 8)
        self.assertEqual(remove_whitespace('CCGGGTTT'), strand.dna_sequence)
        strand = strand_matching(self.small_design.strands, 1, False, 8, 16)
        self.assertEqual(remove_whitespace('GGCCCAAA'), strand.dna_sequence)

    def test_add_full_crossover__small_design_H0_forward(self):
        """
        0        8        16
        ACGTACGA AACCGGTA
    0   [------+ +------>
        <------- -------]
        TGCATGCT TTGGCCAT
               | |
        AAACCCGG TTTGGGCC
    1   [------- ------->
        <------+ +------]
        TTTGGGCC AAACCCGG
        """
        self.small_design.add_full_crossover(helix1=0, helix2=1, offset1=8, forward1=True)
        self.assertEqual(4, len(self.small_design.strands))
        # two new Strands
        self.assertIn(sc.Strand([
            sc.Substrand(0, True, 0, 8),
            sc.Substrand(1, False, 0, 8),
        ]), self.small_design.strands)
        self.assertIn(sc.Strand([
            sc.Substrand(1, False, 8, 16),
            sc.Substrand(0, True, 8, 16),
        ]), self.small_design.strands)
        # existing Strands
        self.assertIn(sc.Strand([sc.Substrand(0, False, 0, 16)]), self.small_design.strands)
        self.assertIn(sc.Strand([sc.Substrand(1, True, 0, 16)]), self.small_design.strands)
        # DNA
        strand = strand_matching(self.small_design.strands, 0, False, 0, 16)
        self.assertEqual(remove_whitespace('TACCGGTT TCGTACGT'), strand.dna_sequence)
        strand = strand_matching(self.small_design.strands, 1, True, 0, 16)
        self.assertEqual(remove_whitespace('AAACCCGG TTTGGGCC'), strand.dna_sequence)
        strand = strand_matching(self.small_design.strands, 0, True, 0, 8)
        self.assertEqual(remove_whitespace('ACGTACGA CCGGGTTT'), strand.dna_sequence)
        strand = strand_matching(self.small_design.strands, 1, False, 8, 16)
        self.assertEqual(remove_whitespace('GGCCCAAA AACCGGTA'), strand.dna_sequence)

    def test_add_full_crossover__small_design_H0_reverse(self):
        """
        0        8        16
        ACGTACGA AACCGGTA
    0   [------- ------->
        <------+ +------]
        TGCATGCT TTGGCCAT
               | |
        AAACCCGG TTTGGGCC
    1   [------+ +------>
        <------- -------]
        TTTGGGCC AAACCCGG
        """
        self.small_design.add_full_crossover(helix1=0, helix2=1, offset1=8, forward1=False)
        self.assertEqual(4, len(self.small_design.strands))
        # two new Strands
        self.assertIn(sc.Strand([
            sc.Substrand(1, True, 0, 8),
            sc.Substrand(0, False, 0, 8),
        ]), self.small_design.strands)
        self.assertIn(sc.Strand([
            sc.Substrand(0, False, 8, 16),
            sc.Substrand(1, True, 8, 16),
        ]), self.small_design.strands)
        # existing Strands
        self.assertIn(sc.Strand([sc.Substrand(0, True, 0, 16)]), self.small_design.strands)
        self.assertIn(sc.Strand([sc.Substrand(1, False, 0, 16)]), self.small_design.strands)
        # DNA
        strand = strand_matching(self.small_design.strands, 0, True, 0, 16)
        self.assertEqual(remove_whitespace('ACGTACGA AACCGGTA'), strand.dna_sequence)
        strand = strand_matching(self.small_design.strands, 0, False, 8, 16)
        self.assertEqual(remove_whitespace('TACCGGTT TTTGGGCC'), strand.dna_sequence)
        strand = strand_matching(self.small_design.strands, 1, True, 0, 8)
        self.assertEqual(remove_whitespace('AAACCCGG TCGTACGT'), strand.dna_sequence)
        strand = strand_matching(self.small_design.strands, 1, False, 0, 16)
        self.assertEqual(remove_whitespace('GGCCCAAA CCGGGTTT'), strand.dna_sequence)

    def test_add_half_crossover__small_design_H0_reverse_8(self):
        """
        0        8        16
        ACGTACGA AACCGGTA
    0   [------- ------->
        <------] +------]
        TGCATGCT TTGGCCAT
                 |
        AAACCCGG TTTGGGCC
    1   [------> +------>
        <------- -------]
        TTTGGGCC AAACCCGG
        """
        self.small_design.add_nick(helix=0, offset=8, forward=False)
        self.small_design.add_nick(helix=1, offset=8, forward=True)
        self.small_design.add_half_crossover(helix1=0, helix2=1, offset1=8, forward1=False)
        self.assertEqual(5, len(self.small_design.strands))
        # three new Strands
        self.assertIn(sc.Strand([
            sc.Substrand(1, True, 0, 8),
        ]), self.small_design.strands)
        self.assertIn(sc.Strand([
            sc.Substrand(0, False, 0, 8),
        ]), self.small_design.strands)
        self.assertIn(sc.Strand([
            sc.Substrand(0, False, 8, 16),
            sc.Substrand(1, True, 8, 16),
        ]), self.small_design.strands)
        # existing Strands
        self.assertIn(sc.Strand([sc.Substrand(0, True, 0, 16)]), self.small_design.strands)
        self.assertIn(sc.Strand([sc.Substrand(1, False, 0, 16)]), self.small_design.strands)
        # DNA
        strand = strand_matching(self.small_design.strands, 0, True, 0, 16)
        self.assertEqual(remove_whitespace('ACGTACGA AACCGGTA'), strand.dna_sequence)
        strand = strand_matching(self.small_design.strands, 0, False, 8, 16)
        self.assertEqual(remove_whitespace('TACCGGTT TTTGGGCC'), strand.dna_sequence)
        strand = strand_matching(self.small_design.strands, 0, False, 0, 8)
        self.assertEqual(remove_whitespace('TCGTACGT'), strand.dna_sequence)
        strand = strand_matching(self.small_design.strands, 1, True, 0, 8)
        self.assertEqual(remove_whitespace('AAACCCGG'), strand.dna_sequence)
        strand = strand_matching(self.small_design.strands, 1, False, 0, 16)
        self.assertEqual(remove_whitespace('GGCCCAAA CCGGGTTT'), strand.dna_sequence)

    def test_add_half_crossover__small_design_H0_reverse_0(self):
        """
        0        8        16
        ACGTACGA AACCGGTA
    0   [------- ------->
        +------- -------]
        TGCATGCT TTGGCCAT
        |
        AAACCCGG TTTGGGCC
    1   -------- ------->
        <------- -------]
        TTTGGGCC AAACCCGG
        """
        self.small_design.add_half_crossover(helix1=0, helix2=1, offset1=0, forward1=False)
        self.assertEqual(3, len(self.small_design.strands))
        # one new Strand
        self.assertIn(sc.Strand([
            sc.Substrand(0, False, 0, 16),
            sc.Substrand(1, True, 0, 16),
        ]), self.small_design.strands)
        # existing Strands
        self.assertIn(sc.Strand([sc.Substrand(0, True, 0, 16)]), self.small_design.strands)
        self.assertIn(sc.Strand([sc.Substrand(1, False, 0, 16)]), self.small_design.strands)
        # DNA
        strand = strand_matching(self.small_design.strands, 0, True, 0, 16)
        self.assertEqual(remove_whitespace('ACGTACGA AACCGGTA'), strand.dna_sequence)
        strand = strand_matching(self.small_design.strands, 0, False, 0, 16)
        self.assertEqual(remove_whitespace('TACCGGTT TCGTACGT AAACCCGG TTTGGGCC'), strand.dna_sequence)
        strand = strand_matching(self.small_design.strands, 1, False, 0, 16)
        self.assertEqual(remove_whitespace('GGCCCAAA CCGGGTTT'), strand.dna_sequence)

    def test_add_half_crossover__small_design_H0_reverse_15(self):
        """
        0        8        16
        ACGTACGA AACCGGTA
    0   [------- ------->
        <------- -------+
        TGCATGCT TTGGCCAT
                        |
        AAACCCGG TTTGGGCC
    1   [------- -------+
        <------- -------]
        TTTGGGCC AAACCCGG
        """
        self.small_design.add_half_crossover(helix1=0, helix2=1, offset1=15, forward1=False)
        self.assertEqual(3, len(self.small_design.strands))
        # one new Strand
        self.assertIn(sc.Strand([
            sc.Substrand(1, True, 0, 16),
            sc.Substrand(0, False, 0, 16),
        ]), self.small_design.strands)
        # existing Strands
        self.assertIn(sc.Strand([sc.Substrand(0, True, 0, 16)]), self.small_design.strands)
        self.assertIn(sc.Strand([sc.Substrand(1, False, 0, 16)]), self.small_design.strands)
        # DNA
        strand = strand_matching(self.small_design.strands, 0, True, 0, 16)
        self.assertEqual(remove_whitespace('ACGTACGA AACCGGTA'), strand.dna_sequence)
        strand = strand_matching(self.small_design.strands, 1, True, 0, 16)
        self.assertEqual(remove_whitespace('AAACCCGG TTTGGGCC TACCGGTT TCGTACGT'), strand.dna_sequence)
        strand = strand_matching(self.small_design.strands, 1, False, 0, 16)
        self.assertEqual(remove_whitespace('GGCCCAAA CCGGGTTT'), strand.dna_sequence)

    def test_add_half_crossover__small_design_illegal(self):
        """
        0        8        16
    0   [------- ------->
        <------- -------- ?
                          |
    1   [------- -------- ?
        <------- -------]
        """
        with self.assertRaises(sc.IllegalDNADesignError):
            self.small_design.add_half_crossover(helix1=0, helix2=1, offset1=16, forward1=False)

    def test_add_full_crossover__small_design_illegal(self):
        """
        0        8        16
    0   [------- ------->
        <------- -------+ ?
                        | |
    1   [------- -------+ ?
        <------- -------]
        """
        with self.assertRaises(sc.IllegalDNADesignError):
            self.small_design.add_full_crossover(helix1=0, helix2=1, offset1=16, forward1=False)

    def test_add_full_crossover__small_design_illegal_only_one_helix_has_substrand(self):
        """
        0        8        16
    0   [------- ------->
        <------+ +------]
               | |
    1   [--->  ? ?
        <---]
        """
        design = sc.DNADesign(strands=[
            sc.Strand([sc.Substrand(0, True, 0, 16)]),
            sc.Strand([sc.Substrand(0, False, 0, 16)]),
            sc.Strand([sc.Substrand(1, True, 0, 5)]),
            sc.Strand([sc.Substrand(1, False, 0, 5)]),
        ], grid=sc.square)
        with self.assertRaises(sc.IllegalDNADesignError):
            design.add_full_crossover(helix1=0, helix2=1, offset1=10, forward1=False)

    r"""
    0        8        16       24       32       40       48       56       64       72       80       88       96
0   [------- -------- -------- -------- -------- -------- -------- -------- -------- -------- -------- ------->
    <------- -------- -------- -------- -------] <------- -------- -------- -------] <------- -------- -------]

1   [------- -------- -------> [------- -------- -------- -------> [------- -------- -------- -------- ------->
    <------- -------- -------- -------- -------- -------- -------- -------- -------- -------- -------- -------]

2   [------- -------- -------- -------- -------- -------- -------- -------- -------- -------- -------- ------->
    <------- -------- -------- -------- -------] <------- -------- -------- -------] <------- -------- -------]

3   [------- -------- -------- -------- -------- -------- -------- -------- -------- -------- -------- ------->
    <------- -------- -------> [------- -------- -------- -------> [------- -------- -------- -------- -------]

4   [------- -------- -------- -------- -------- -------- -------- -------- -------- -------- -------- ------->
    <------- -------- -------- -------- -------] <------- -------- -------- -------] <------- -------- -------]

5   [------- -------- -------> [------- -------- -------- -------> [------- -------- -------- -------- ------->
    <------- -------- -------- -------- -------- -------] <------- -------- -------- -------- -------- -------]
    """

    def add_nicks(self, design: sc.DNADesign):
        design.add_nick(helix=5, offset=48, forward=False)
        design.add_nick(helix=0, offset=40, forward=False)
        design.add_nick(helix=0, offset=72, forward=False)
        design.add_nick(helix=2, offset=40, forward=False)
        design.add_nick(helix=2, offset=72, forward=False)
        design.add_nick(helix=4, offset=40, forward=False)
        design.add_nick(helix=4, offset=72, forward=False)
        design.add_nick(helix=1, offset=24, forward=True)
        design.add_nick(helix=1, offset=56, forward=True)
        design.add_nick(helix=3, offset=24, forward=True)
        design.add_nick(helix=3, offset=56, forward=True)
        design.add_nick(helix=5, offset=24, forward=True)
        design.add_nick(helix=5, offset=56, forward=True)

    def test_add_nick__6_helix_rectangle(self):
        self.add_nicks(self.design)
        self.assertEqual(25, len(self.design.strands))
        for helix in range(0, len(self.design.helices), 2):
            # even helix
            self.assertIn(sc.Strand([sc.Substrand(helix, True, 0, 96)]), self.design.strands)
            self.assertIn(sc.Strand([sc.Substrand(helix, False, 0, 40)]), self.design.strands)
            self.assertIn(sc.Strand([sc.Substrand(helix, False, 40, 72)]), self.design.strands)
            self.assertIn(sc.Strand([sc.Substrand(helix, False, 72, 96)]), self.design.strands)
            # odd helix
            if helix + 1 < len(self.design.helices) - 1:
                self.assertIn(sc.Strand([sc.Substrand(helix + 1, False, 0, 96)]), self.design.strands)
            else:
                # nick in scaffold on bottom helix
                self.assertIn(sc.Strand([sc.Substrand(helix + 1, False, 0, 48)]), self.design.strands)
                self.assertIn(sc.Strand([sc.Substrand(helix + 1, False, 48, 96)]), self.design.strands)
            self.assertIn(sc.Strand([sc.Substrand(helix + 1, True, 0, 24)]), self.design.strands)
            self.assertIn(sc.Strand([sc.Substrand(helix + 1, True, 24, 56)]), self.design.strands)
            self.assertIn(sc.Strand([sc.Substrand(helix + 1, True, 56, 96)]), self.design.strands)

    # TODO: re-write this test after support for circular Strands is added and test making crossovers first
    r"""
    Crossovers needed to be added after nicks so that there will be no circular strands 
    0        8        16       24       32       40       48       56       64       72       80       88       96
0   +------- -------- -------- -------- -------- -------- -------- -------- -------- -------- -------- -------+
   /<------- -------+ +------- -------- -------] <------- -------- -------- -------] <------+ +------- -------]\
  (                 | |                                                                     | |                 )
1  \[------- -------+ +------> [------+ +------- -------- -------> [------+ +------- -------+ +------- ------->/
    +------- -------- -------- -------- -------- -------+ +------- -------- -------- -------- -------- -------+
                                      | |               | |               | |
2   +------- -------- -------- -------- -------- -------+ +------- -------- -------- -------- -------- -------+
   /<------- -------+ +------- -------+ +------] <------- -------- -------+ +------] <------+ +------- -------]\
  (                 | |                                                                     | |                 )
3  \[------- -------+ +------> [------+ +------- -------- -------> [------+ +------- -------+ +------- ------->/
    +------- -------- -------- -------- -------- -------+ +------- -------- -------- -------- -------- -------+
                                      | |               | |               | |
4   +------- -------- -------- -------- -------- -------+ +------- -------- -------- -------- -------- -------+
   /<------- -------+ +------- -------+ +------] <------- -------- -------+ +------] <------+ +------- -------]\
  (                 | |                                                                     | |                 )
5  \[------- -------+ +------> [------- -------- -------- -------> [------- -------- -------+ +------- ------->/
    +------- -------- -------- -------- -------- -------] <------- -------- -------- -------- -------- -------+

    """

    def add_crossovers_after_nicks(self, design: sc.DNADesign):
        # scaffold seam crossovers
        design.add_full_crossover(helix1=1, helix2=2, offset1=48, forward1=False)
        design.add_full_crossover(helix1=3, helix2=4, offset1=48, forward1=False)

        # staple crossovers
        design.add_full_crossover(helix1=0, helix2=1, offset1=16, forward1=False)
        design.add_full_crossover(helix1=0, helix2=1, offset1=80, forward1=False)
        design.add_full_crossover(helix1=1, helix2=2, offset1=32, forward1=True)
        design.add_full_crossover(helix1=1, helix2=2, offset1=64, forward1=True)
        design.add_full_crossover(helix1=2, helix2=3, offset1=16, forward1=False)
        design.add_full_crossover(helix1=2, helix2=3, offset1=80, forward1=False)
        design.add_full_crossover(helix1=3, helix2=4, offset1=32, forward1=True)
        design.add_full_crossover(helix1=3, helix2=4, offset1=64, forward1=True)
        design.add_full_crossover(helix1=4, helix2=5, offset1=16, forward1=False)
        design.add_full_crossover(helix1=4, helix2=5, offset1=80, forward1=False)

        # The left and right edge crossovers need to be added last to ensure the Strands remain
        # non-circular during all intermediate stages.

        # scaffold left crossovers
        design.add_half_crossover(helix1=0, helix2=1, offset1=0, forward1=True)
        design.add_half_crossover(helix1=2, helix2=3, offset1=0, forward1=True)
        design.add_half_crossover(helix1=4, helix2=5, offset1=0, forward1=True)

        # scaffold right crossovers
        design.add_half_crossover(helix1=0, helix2=1, offset1=95, forward1=True)
        design.add_half_crossover(helix1=2, helix2=3, offset1=95, forward1=True)
        design.add_half_crossover(helix1=4, helix2=5, offset1=95, forward1=True)

    def test_add_nick_then_add_crossovers__6_helix_rectangle(self):
        self.add_nicks(self.design)
        self.add_crossovers_after_nicks(self.design)

        self.assertEqual(19, len(self.design.strands))

        # staples left edge
        # {"helix": 1, "forward": true, "start": 0, "end": 16},
        # {"helix": 0, "forward": false, "start": 0, "end": 16}
        stap = sc.Strand([
            sc.Substrand(1, True, 0, 16),
            sc.Substrand(0, False, 0, 16),
        ])
        self.assertIn(stap, self.design.strands)

        # {"helix": 3, "forward": true, "start": 0, "end": 16},
        # {"helix": 2, "forward": false, "start": 0, "end": 16}
        stap = sc.Strand([
            sc.Substrand(3, True, 0, 16),
            sc.Substrand(2, False, 0, 16),
        ])
        self.assertIn(stap, self.design.strands)

        # {"helix": 5, "forward": true, "start": 0, "end": 16},
        # {"helix": 4, "forward": false, "start": 0, "end": 16}
        stap3 = sc.Strand([
            sc.Substrand(5, True, 0, 16),
            sc.Substrand(4, False, 0, 16),
        ])
        self.assertIn(stap, self.design.strands)

        # staples right edge
        # {"helix": 0, "forward": false, "start": 80, "end": 96},
        # {"helix": 1, "forward": true, "start": 80, "end": 96}
        stap = sc.Strand([
            sc.Substrand(0, False, 80, 96),
            sc.Substrand(1, True, 80, 96),
        ])
        self.assertIn(stap, self.design.strands)

        # {"helix": 2, "forward": false, "start": 80, "end": 96},
        # {"helix": 3, "forward": true, "start": 80, "end": 96}
        stap = sc.Strand([
            sc.Substrand(2, False, 80, 96),
            sc.Substrand(3, True, 80, 96),
        ])
        self.assertIn(stap, self.design.strands)

        # {"helix": 4, "forward": false, "start": 80, "end": 96},
        # {"helix": 5, "forward": true, "start": 80, "end": 96}
        stap = sc.Strand([
            sc.Substrand(4, False, 80, 96),
            sc.Substrand(5, True, 80, 96),
        ])
        self.assertIn(stap, self.design.strands)

        # staples remainder
        # {"helix": 0, "forward": false, "start": 40, "end": 72}
        stap = sc.Strand([sc.Substrand(0, False, 40, 72)])
        self.assertIn(stap, self.design.strands)

        # {"helix": 2, "forward": false, "start": 32, "end": 40},
        # {"helix": 1, "forward": true, "start": 32, "end": 56}
        stap = sc.Strand([
            sc.Substrand(2, False, 32, 40),
            sc.Substrand(1, True, 32, 56),
        ])
        self.assertIn(stap, self.design.strands)

        # {"helix": 1, "forward": true, "start": 56, "end": 64},
        # {"helix": 2, "forward": false, "start": 40, "end": 64}
        stap = sc.Strand([
            sc.Substrand(1, True, 56, 64),
            sc.Substrand(2, False, 40, 64),
        ])
        self.assertIn(stap, self.design.strands)

        # {"helix": 4, "forward": false, "start": 32, "end": 40},
        # {"helix": 3, "forward": true, "start": 32, "end": 56}
        stap = sc.Strand([
            sc.Substrand(4, False, 32, 40),
            sc.Substrand(3, True, 32, 56),
        ])
        self.assertIn(stap, self.design.strands)

        # {"helix": 3, "forward": true, "start": 56, "end": 64},
        # {"helix": 4, "forward": false, "start": 40, "end": 64}
        stap = sc.Strand([
            sc.Substrand(3, True, 56, 64),
            sc.Substrand(4, False, 40, 64),
        ])
        self.assertIn(stap, self.design.strands)

        # {"helix": 5, "forward": true, "start": 24, "end": 56}
        stap = sc.Strand([sc.Substrand(5, True, 24, 56)])
        self.assertIn(stap, self.design.strands)

        # {"helix": 0, "forward": false, "start": 16, "end": 40},
        # {"helix": 1, "forward": true, "start": 16, "end": 24}
        stap = sc.Strand([
            sc.Substrand(0, False, 16, 40),
            sc.Substrand(1, True, 16, 24),
        ])
        self.assertIn(stap, self.design.strands)

        # {"helix": 1, "forward": true, "start": 24, "end": 32},
        # {"helix": 2, "forward": false, "start": 16, "end": 32},
        # {"helix": 3, "forward": true, "start": 16, "end": 24}
        stap = sc.Strand([
            sc.Substrand(1, True, 24, 32),
            sc.Substrand(2, False, 16, 32),
            sc.Substrand(3, True, 16, 24),
        ])
        self.assertIn(stap, self.design.strands)

        # {"helix": 3, "forward": true, "start": 24, "end": 32},
        # {"helix": 4, "forward": false, "start": 16, "end": 32},
        # {"helix": 5, "forward": true, "start": 16, "end": 24}
        stap = sc.Strand([
            sc.Substrand(3, True, 24, 32),
            sc.Substrand(4, False, 16, 32),
            sc.Substrand(5, True, 16, 24),
        ])
        self.assertIn(stap, self.design.strands)

        # {"helix": 5, "forward": true, "start": 56, "end": 80},
        # {"helix": 4, "forward": false, "start": 72, "end": 80}
        stap = sc.Strand([
            sc.Substrand(5, True, 56, 80),
            sc.Substrand(4, False, 72, 80),
        ])
        self.assertIn(stap, self.design.strands)

        # {"helix": 2, "forward": false, "start": 64, "end": 72},
        # {"helix": 1, "forward": true, "start": 64, "end": 80},
        # {"helix": 0, "forward": false, "start": 72, "end": 80}
        stap = sc.Strand([
            sc.Substrand(2, False, 64, 72),
            sc.Substrand(1, True, 64, 80),
            sc.Substrand(0, False, 72, 80),
        ])
        self.assertIn(stap, self.design.strands)

        # {"helix": 4, "forward": false, "start": 64, "end": 72},
        # {"helix": 3, "forward": true, "start": 64, "end": 80},
        # {"helix": 2, "forward": false, "start": 72, "end": 80}
        stap = sc.Strand([
            sc.Substrand(4, False, 64, 72),
            sc.Substrand(3, True, 64, 80),
            sc.Substrand(2, False, 72, 80),
        ])
        self.assertIn(stap, self.design.strands)

        # scaffold
        #     {"helix": 5, "forward": false, "start": 0, "end": 48},
        #     {"helix": 4, "forward": true, "start": 0, "end": 48},
        #     {"helix": 3, "forward": false, "start": 0, "end": 48},
        #     {"helix": 2, "forward": true, "start": 0, "end": 48},
        #     {"helix": 1, "forward": false, "start": 0, "end": 48},
        #     {"helix": 0, "forward": true, "start": 0, "end": 96},
        #     {"helix": 1, "forward": false, "start": 48, "end": 96},
        #     {"helix": 2, "forward": true, "start": 48, "end": 96},
        #     {"helix": 3, "forward": false, "start": 48, "end": 96},
        #     {"helix": 4, "forward": true, "start": 48, "end": 96},
        #     {"helix": 5, "forward": false, "start": 48, "end": 96}
        scaf = sc.Strand([
            sc.Substrand(5, False, 0, 48),
            sc.Substrand(4, True, 0, 48),
            sc.Substrand(3, False, 0, 48),
            sc.Substrand(2, True, 0, 48),
            sc.Substrand(1, False, 0, 48),
            sc.Substrand(0, True, 0, 96),
            sc.Substrand(1, False, 48, 96),
            sc.Substrand(2, True, 48, 96),
            sc.Substrand(3, False, 48, 96),
            sc.Substrand(4, True, 48, 96),
            sc.Substrand(5, False, 48, 96),
        ])
        self.assertIn(scaf, self.design.strands)


class TestAutocalculatedData(unittest.TestCase):

    def test_helix_min_max_offsets_illegal_explicitly_specified(self):
        helices = [sc.Helix(min_offset=5, max_offset=5)]
        with self.assertRaises(sc.IllegalDNADesignError):
            design = sc.DNADesign(helices=helices, strands=[], grid=sc.square)

    def test_helix_min_max_offsets_illegal_autocalculated(self):
        helices = [sc.Helix(min_offset=5)]
        ss = sc.Substrand(0, True, 0, 4)
        strand = sc.Strand([ss])
        with self.assertRaises(sc.IllegalDNADesignError):
            design = sc.DNADesign(helices=helices, strands=[strand], grid=sc.square)

    def test_helix_min_max_offsets(self):
        helices = [sc.Helix(), sc.Helix(min_offset=-5), sc.Helix(max_offset=5),
                   sc.Helix(min_offset=5, max_offset=10)]
        ss_0 = sc.Substrand(helix=0, forward=True, start=20, end=25)
        ss_1 = sc.Substrand(helix=1, forward=False, start=-5, end=30)
        ss_2 = sc.Substrand(helix=2, forward=True, start=0, end=5)
        ss_3 = sc.Substrand(helix=3, forward=False, start=5, end=10)
        strand = sc.Strand([ss_0, ss_1, ss_2, ss_3])
        design = sc.DNADesign(helices=helices, strands=[strand], grid=sc.square)
        self.assertEqual(0, design.helices[0].min_offset)
        self.assertEqual(25, design.helices[0].max_offset)
        self.assertEqual(-5, design.helices[1].min_offset)
        self.assertEqual(30, design.helices[1].max_offset)
        self.assertEqual(0, design.helices[2].min_offset)
        self.assertEqual(5, design.helices[2].max_offset)
        self.assertEqual(5, design.helices[3].min_offset)
        self.assertEqual(10, design.helices[3].max_offset)

    def test_helix_max_offset(self):
        helices = [sc.Helix(), sc.Helix(max_offset=8), sc.Helix()]
        ss_0 = sc.Substrand(helix=0, forward=True, start=5, end=10)
        ss_1 = sc.Substrand(helix=1, forward=False, start=2, end=6)
        ss_2 = sc.Substrand(helix=2, forward=True, start=0, end=5)
        strand = sc.Strand([ss_0, ss_1, ss_2])
        design = sc.DNADesign(helices=helices, strands=[strand], grid=sc.square)
        self.assertEqual(10, design.helices[0].max_offset)
        self.assertEqual(8, design.helices[1].max_offset)
        self.assertEqual(5, design.helices[2].max_offset)


class TestJSON(unittest.TestCase):
    def test_to_json__hairpin(self):
        """
        01234
        AAACC    # helix 0 going forward
             \
              T  # loopout
              G  # loopout
              C  # loopout
              A  # loopout
              C  # loopout
             /
        TTTGG    # helix 0 going reverse
        """
        ss_f = sc.Substrand(helix=0, forward=True, start=0, end=5)
        loop = sc.Loopout(length=5)
        ss_r = sc.Substrand(helix=0, forward=False, start=0, end=5)
        strand_forward = sc.Strand([ss_f, loop, ss_r])
        design = sc.DNADesign(strands=[strand_forward], grid=sc.square)
        design.assign_dna(strand_forward, 'AAACC TGCAC')
        json = design.to_json()
        # should be no error getting here

    def test_to_json__rotation(self):
        helix = sc.Helix(rotation=math.pi / 2, rotation_anchor=31)
        ss_f = sc.Substrand(helix=0, forward=True, start=0, end=5)
        ss_r = sc.Substrand(helix=0, forward=False, start=0, end=5)
        strand_f = sc.Strand([ss_f])
        strand_r = sc.Strand([ss_r])
        design = sc.DNADesign(helices=[helix], strands=[strand_f, strand_r], grid=sc.square)
        json = design.to_json()
        # should be no error getting here


class TestIllegalStructuresPrevented(unittest.TestCase):

    def test_consecutive_substrands_loopout(self):
        helices = [sc.Helix(max_offset=10)]
        ss1 = sc.Substrand(0, sc.forward, 0, 3)
        ss2 = sc.Loopout(4)
        ss3 = sc.Loopout(4)
        with self.assertRaises(sc.IllegalDNADesignError):
            strand = sc.Strand([ss1, ss2, ss3])

        strand = sc.Strand([ss1, ss2])
        strand.substrands.append(ss3)
        with self.assertRaises(sc.IllegalDNADesignError):
            design = sc.DNADesign(helices=helices, strands=[strand], grid=sc.square)

    def test_singleton_loopout(self):
        helices = [sc.Helix(max_offset=10)]
        ss1 = sc.Loopout(4)
        with self.assertRaises(sc.IllegalDNADesignError):
            strand = sc.Strand([ss1])

        strand = sc.Strand([])
        strand.substrands.append(ss1)
        with self.assertRaises(sc.IllegalDNADesignError):
            design = sc.DNADesign(helices=helices, strands=[strand], grid=sc.square)

    def test_strand_offset_beyond_maxbases(self):
        helices = [sc.Helix(max_offset=10)]
        ss1 = sc.Substrand(0, sc.forward, 0, 20)
        strands = [sc.Strand([ss1])]
        with self.assertRaises(sc.StrandError):
            design = sc.DNADesign(helices=helices, strands=strands)

    def test_to_idt_bulk_input_format__duplicate_names_same_sequence(self):
        length = 8
        helices = [sc.Helix(max_offset=length)]
        ss1_r = sc.Substrand(0, sc.forward, 0, 4)
        ss2_r = sc.Substrand(0, sc.forward, 4, 8)
        ss_l = sc.Substrand(0, sc.reverse, 0, 4)

        s1_r = sc.Strand([ss1_r], idt=sc.IDTFields('s1_r'))
        s2_r = sc.Strand([ss2_r], idt=sc.IDTFields('s1_r'))
        s_l = sc.Strand([ss_l], idt=sc.IDTFields('s_l'))

        strands = [s1_r, s2_r, s_l]

        design = sc.DNADesign(helices=helices, strands=strands, grid=sc.square)

        design.assign_dna(s_l, 'AGTT')
        design.assign_dna(s2_r, 'AACT')

        # should not raise exception
        idt_str = design.to_idt_bulk_input_format()

    def test_to_idt_bulk_input_format__duplicate_names_different_sequences(self):
        ss1_r = sc.Substrand(0, sc.forward, 0, 4)
        ss2_r = sc.Substrand(0, sc.forward, 4, 8)
        ss_l = sc.Substrand(0, sc.reverse, 0, 4)

        s1_r = sc.Strand([ss1_r], idt=sc.IDTFields('s1_r'))
        s2_r = sc.Strand([ss2_r], idt=sc.IDTFields('s1_r'))
        s_l = sc.Strand([ss_l], idt=sc.IDTFields('s_l'))

        strands = [s1_r, s2_r, s_l]

        design = sc.DNADesign(strands=strands, grid=sc.square)

        design.assign_dna(s_l, 'AGTT')
        design.assign_dna(s2_r, 'GGGG')

        with self.assertRaises(sc.IllegalDNADesignError):
            idt_str = design.to_idt_bulk_input_format()

    def test_to_idt_bulk_input_format__duplicate_names_different_scales(self):
        ss1_r = sc.Substrand(0, sc.forward, 0, 4)
        ss2_r = sc.Substrand(0, sc.forward, 4, 8)
        ss_l = sc.Substrand(0, sc.reverse, 0, 4)

        s1_r = sc.Strand([ss1_r], idt=sc.IDTFields('s1_r', scale='25nm'))
        s2_r = sc.Strand([ss2_r], idt=sc.IDTFields('s1_r', scale='100nm'))
        s_l = sc.Strand([ss_l], idt=sc.IDTFields('s_l'))

        strands = [s1_r, s2_r, s_l]

        design = sc.DNADesign(strands=strands, grid=sc.square)

        design.assign_dna(s_l, 'AGTT')
        design.assign_dna(s2_r, 'AACT')

        with self.assertRaises(sc.IllegalDNADesignError):
            idt_str = design.to_idt_bulk_input_format()

    def test_to_idt_bulk_input_format__duplicate_names_different_purifications(self):
        length = 8
        ss1_r = sc.Substrand(0, sc.forward, 0, 4)
        ss2_r = sc.Substrand(0, sc.forward, 4, 8)
        ss_l = sc.Substrand(0, sc.reverse, 0, 4)

        s1_r = sc.Strand([ss1_r], idt=sc.IDTFields('s1_r', purification='STD'))
        s2_r = sc.Strand([ss2_r], idt=sc.IDTFields('s1_r', purification='HPLC'))
        s_l = sc.Strand([ss_l], idt=sc.IDTFields('s_l'))

        strands = [s1_r, s2_r, s_l]

        design = sc.DNADesign(strands=strands, grid=sc.square)

        design.assign_dna(s_l, 'AGTT')
        design.assign_dna(s2_r, 'AACT')

        with self.assertRaises(sc.IllegalDNADesignError):
            idt_str = design.to_idt_bulk_input_format()

    def test_assign_dna__conflicting_sequences_directly_assigned(self):
        ss_right = sc.Substrand(0, sc.forward, 0, 5)
        ss_left = sc.Substrand(0, sc.reverse, 0, 5)
        strand_right = sc.Strand([ss_right])
        strand_left = sc.Strand([ss_left])
        design = sc.DNADesign(strands=[strand_left, strand_right])
        design.assign_dna(strand_right, 'ACGTT')
        with self.assertRaises(sc.IllegalDNADesignError):
            design.assign_dna(strand_right, 'TTTTT')

    def test_assign_dna__conflicting_sequences_indirectly_assigned(self):
        ss_right = sc.Substrand(0, sc.forward, 0, 5)
        ss_left = sc.Substrand(0, sc.reverse, 0, 5)
        strand_right = sc.Strand([ss_right])
        strand_left = sc.Strand([ss_left])
        design = sc.DNADesign(strands=[strand_left, strand_right])
        design.assign_dna(strand_right, 'ACGTT')
        with self.assertRaises(sc.IllegalDNADesignError):
            design.assign_dna(strand_left, 'GGGGG')

    def test_overlapping_caught_in_strange_counterexample(self):
        # found this counterexample as a simplified version of something caught in practice
        s1_left_ss0 = sc.Substrand(0, sc.reverse, 0, 5)
        s1_ss1 = sc.Substrand(0, sc.forward, 0, 15)
        s1_right_ss0 = sc.Substrand(0, sc.reverse, 5, 15)
        s1 = sc.Strand([s1_left_ss0, s1_ss1, s1_right_ss0])

        s2_ss1 = sc.Substrand(0, sc.forward, 10, 20)
        s2_ss0 = sc.Substrand(0, sc.reverse, 10, 20)
        s2 = sc.Strand([s2_ss1, s2_ss0])

        strands = [s1, s2]

        with self.assertRaises(sc.IllegalDNADesignError):
            design = sc.DNADesign(strands=strands, grid=sc.square)

    def test_major_tick_outside_range(self):
        with self.assertRaises(sc.IllegalDNADesignError):
            helix = sc.Helix(max_offset=9, major_ticks=[2, 5, 10])

    def test_major_tick_just_inside_range(self):
        helix = sc.Helix(max_offset=9, major_ticks=[0, 5, 9])

    def test_two_illegally_overlapping_strands(self):
        ss_bot = sc.Substrand(helix=0, forward=sc.reverse, start=0, end=9)
        ss_top = sc.Substrand(helix=0, forward=sc.reverse, start=0, end=9)
        strand_bot = sc.Strand(substrands=[ss_bot])
        strand_top = sc.Strand(substrands=[ss_top])
        strands = [strand_bot, strand_top]
        with self.assertRaises(sc.IllegalDNADesignError):
            sc.DNADesign(grid=sc.square, strands=strands)

    def test_two_nonconsecutive_illegally_overlapping_strands(self):
        ss_top1 = sc.Substrand(helix=0, forward=sc.reverse, start=0, end=5)
        ss_bot = sc.Substrand(helix=0, forward=sc.forward, start=2, end=9)
        ss_top2 = sc.Substrand(helix=0, forward=sc.reverse, start=4, end=8)
        strand_bot = sc.Strand(substrands=[ss_bot])
        strand_top1 = sc.Strand(substrands=[ss_top1])
        strand_top2 = sc.Strand(substrands=[ss_top2])
        strands = [strand_bot, strand_top1, strand_top2]
        with self.assertRaises(sc.IllegalDNADesignError):
            sc.DNADesign(grid=sc.square, strands=strands)

    def test_four_legally_leapfrogging_strands(self):
        ss_top1 = sc.Substrand(helix=0, forward=sc.reverse, start=0, end=20)
        ss_bot1 = sc.Substrand(helix=0, forward=sc.forward, start=10, end=30)
        ss_top2 = sc.Substrand(helix=0, forward=sc.reverse, start=20, end=40)
        ss_bot2 = sc.Substrand(helix=0, forward=sc.forward, start=30, end=50)
        strand_bot1 = sc.Strand(substrands=[ss_bot1])
        strand_bot2 = sc.Strand(substrands=[ss_bot2])
        strand_top1 = sc.Strand(substrands=[ss_top1])
        strand_top2 = sc.Strand(substrands=[ss_top2])
        strands = [strand_bot1, strand_bot2, strand_top1, strand_top2]
        sc.DNADesign(grid=sc.square, strands=strands)

    def test_strand_references_nonexistent_helix(self):
        h1 = sc.Helix(max_offset=9)
        h2 = sc.Helix(max_offset=9)
        ss_bot = sc.Substrand(helix=2, forward=sc.reverse, start=0, end=9)
        ss_top = sc.Substrand(helix=3, forward=sc.reverse, start=0, end=9)
        strand_bot = sc.Strand(substrands=[ss_bot])
        strand_top = sc.Strand(substrands=[ss_top])
        strands = [strand_bot, strand_top]
        with self.assertRaises(sc.IllegalDNADesignError):
            sc.DNADesign(grid=sc.square, helices=[h1, h2], strands=strands)


class TestAddStrand(unittest.TestCase):

    def test_add_strand__with_loopout(self):
        helices = [sc.Helix(max_offset=10), sc.Helix(max_offset=10)]
        design = sc.DNADesign(helices=helices, strands=[])

        ss1 = sc.Substrand(0, True, 0, 10)
        loop = sc.Loopout(4)
        ss2 = sc.Substrand(1, False, 0, 10)
        strand = sc.Strand([ss1, loop, ss2])

        design.add_strand(strand)

        self.assertEqual(1, len(design.strands))
        self.assertEqual(strand, design.strands[0])
        self.assertEqual(ss1, design.substrand_at(0, 0, True))
        self.assertEqual(ss2, design.substrand_at(1, 0, False))


class TestAssignDNA(unittest.TestCase):

    def test_assign_dna__hairpin(self):
        """
        01234
        AAACC    # helix 0 going forward
             \
              T  # loopout
              G  # loopout
              C  # loopout
              A  # loopout
              C  # loopout
             /
        TTTGG    # helix 0 going reverse
        """
        ss_f = sc.Substrand(helix=0, forward=True, start=0, end=5)
        loop = sc.Loopout(length=5)
        ss_r = sc.Substrand(helix=0, forward=False, start=0, end=5)
        strand_forward = sc.Strand([ss_f, loop, ss_r])
        design = sc.DNADesign(strands=[strand_forward], grid=sc.square)
        design.assign_dna(strand_forward, 'AAACC TGCAC')
        self.assertEqual('AAACC TGCAC GGTTT'.replace(' ', ''), strand_forward.dna_sequence)

    def test_assign_dna__from_strand_with_loopout(self):
        """
          01234
        <-TTTGG-]
        [-AAACC-    # helix 0
                \
                 T  # loopout
                 G  # loopout
                 C  # loopout
                 A  # loopout
                 C  # loopout
                /
        <-GCTTA-    # helix 1
        [-CGAAT->
        """
        ss_f = sc.Substrand(helix=0, forward=True, start=0, end=5)
        loop = sc.Loopout(length=5)
        ss_r = sc.Substrand(helix=1, forward=False, start=0, end=5)
        strand_multi = sc.Strand([ss_f, loop, ss_r])

        ss_single0 = sc.Substrand(helix=0, forward=False, start=0, end=5)
        strand_single0 = sc.Strand([ss_single0])

        ss_single1 = sc.Substrand(helix=1, forward=True, start=0, end=5)
        strand_single1 = sc.Strand([ss_single1])

        design = sc.DNADesign(strands=[strand_multi, strand_single0, strand_single1], grid=sc.square)

        design.assign_dna(strand_multi, 'AAACC TGCAC ATTCG')

        self.assertEqual('AAACC TGCAC ATTCG'.replace(' ', ''), strand_multi.dna_sequence)
        self.assertEqual('GGTTT'.replace(' ', ''), strand_single0.dna_sequence)
        self.assertEqual('CGAAT'.replace(' ', ''), strand_single1.dna_sequence)

    def test_assign_dna__to_strand_with_loopout(self):
        """
          01234
        <-TTTGG-]
        [-AAACC-    # helix 0
                \
                 ?  # loopout
                 ?  # loopout
                 ?  # loopout
                 ?  # loopout
                 ?  # loopout
                /
        <-GCTTA-    # helix 1
        [-CGAAT->
        """
        ss_f = sc.Substrand(helix=0, forward=True, start=0, end=5)
        loop = sc.Loopout(length=5)
        ss_r = sc.Substrand(helix=1, forward=False, start=0, end=5)
        strand_multi = sc.Strand([ss_f, loop, ss_r])

        ss_single0 = sc.Substrand(helix=0, forward=False, start=0, end=5)
        strand_single0 = sc.Strand([ss_single0])

        ss_single1 = sc.Substrand(helix=1, forward=True, start=0, end=5)
        strand_single1 = sc.Strand([ss_single1])

        design = sc.DNADesign(strands=[strand_multi, strand_single0, strand_single1], grid=sc.square)

        design.assign_dna(strand_single0, 'GGTTT')

        self.assertEqual('AAACC ????? ?????'.replace(' ', ''), strand_multi.dna_sequence)
        self.assertEqual('GGTTT'.replace(' ', ''), strand_single0.dna_sequence)

        design.assign_dna(strand_single1, 'CGAAT')

        self.assertEqual('AAACC ????? ATTCG'.replace(' ', ''), strand_multi.dna_sequence)
        self.assertEqual('GGTTT'.replace(' ', ''), strand_single0.dna_sequence)
        self.assertEqual('CGAAT'.replace(' ', ''), strand_single1.dna_sequence)

    def test_assign_dna__assign_from_strand_multi_other_single(self):
        """
          01234567
        <-TTTG----GACA-]
        +-AAAC->[-CTGT-+   # helix 0
        |              |
        +-GCTT----AGTA-+   # helix 1
        """
        ss_f_left = sc.Substrand(helix=0, forward=True, start=0, end=4)
        ss_f_right = sc.Substrand(helix=0, forward=True, start=4, end=8)
        ss_h1 = sc.Substrand(helix=1, forward=False, start=0, end=8)
        strand_multi = sc.Strand([ss_f_right, ss_h1, ss_f_left])

        ss_single = sc.Substrand(helix=0, forward=False, start=0, end=8)
        strand_single = sc.Strand([ss_single])

        design = sc.DNADesign(strands=[strand_multi, strand_single], grid=sc.square)

        design.assign_dna(strand_multi, 'CTGT ATGA TTCG AAAC')

        self.assertEqual('CTGT ATGA TTCG AAAC'.replace(' ', ''), strand_multi.dna_sequence)
        self.assertEqual('ACAG GTTT'.replace(' ', ''), strand_single.dna_sequence)

    def test_assign_dna__assign_to_strand_multi_other_single(self):
        """
          01234567
        <-TTTG----GACA-]
        +-AAAC->[-CTGT-+   # helix 0
        |              |
        +-????----????-+   # helix 1
        """
        ss_f_left = sc.Substrand(helix=0, forward=True, start=0, end=4)
        ss_f_right = sc.Substrand(helix=0, forward=True, start=4, end=8)
        ss_h1 = sc.Substrand(helix=1, forward=False, start=0, end=8)
        strand_multi = sc.Strand([ss_f_right, ss_h1, ss_f_left])

        ss_single = sc.Substrand(helix=0, forward=False, start=0, end=8)
        strand_single = sc.Strand([ss_single])

        design = sc.DNADesign(strands=[strand_multi, strand_single], grid=sc.square)

        design.assign_dna(strand_single, 'ACAG GTTT')

        self.assertEqual('CTGT ???? ???? AAAC'.replace(' ', ''), strand_multi.dna_sequence)
        self.assertEqual('ACAG GTTT'.replace(' ', ''), strand_single.dna_sequence)

    def test_assign_dna__other_strand_fully_defined_already(self):
        """
        01234567
        [------>
        CAAAGTCG
        GTTT
        <--]
        """
        ss_r = sc.Substrand(helix=0, forward=sc.forward, start=0, end=8)
        ss_l = sc.Substrand(helix=0, forward=sc.reverse, start=0, end=4)
        strand_r = sc.Strand(substrands=[ss_r])
        strand_l = sc.Strand(substrands=[ss_l])
        design = sc.DNADesign(grid=sc.square, strands=[strand_r, strand_l])
        design.assign_dna(strand_r, 'CAAAGTCG')
        design.assign_dna(strand_l, 'TTTG')
        # should not have an error by this point

    def test_assign_dna__other_strand_fully_defined_already_and_other_extends_beyond(self):
        """
        01234567
        [------>
        CAAAGTCG
          TTCA
          <--]
        """
        ss_r = sc.Substrand(helix=0, forward=sc.forward, start=0, end=8)
        ss_l = sc.Substrand(helix=0, forward=sc.reverse, start=2, end=6)
        strand_r = sc.Strand(substrands=[ss_r])
        strand_l = sc.Strand(substrands=[ss_l])
        design = sc.DNADesign(grid=sc.square, strands=[strand_r, strand_l])
        design.assign_dna(strand_r, 'CAAAGTCG')
        design.assign_dna(strand_l, 'ACTT')
        # should not have an error by this point

    def test_assign_dna__other_strand_fully_defined_already_and_self_extends_beyond(self):
        """
        01234567
        [------>
        CAAAGTCG
          TTCA
          <--]
        """
        ss_r = sc.Substrand(helix=0, forward=sc.forward, start=0, end=8)
        ss_l = sc.Substrand(helix=0, forward=sc.reverse, start=2, end=6)
        strand_r = sc.Strand(substrands=[ss_r])
        strand_l = sc.Strand(substrands=[ss_l])
        design = sc.DNADesign(grid=sc.square, strands=[strand_r, strand_l])
        design.assign_dna(strand_l, 'ACTT')
        design.assign_dna(strand_r, 'CAAAGTCG')
        # should not have an error by this point

    def test_assign_dna__two_equal_length_strands_on_one_helix(self):
        """
        01234
        <---]
        CAAAA
        GTTTT
        [--->
        """
        ss_r = sc.Substrand(helix=0, forward=sc.forward, start=0, end=5)
        ss_l = sc.Substrand(helix=0, forward=sc.reverse, start=0, end=5)
        strand_r = sc.Strand(substrands=[ss_r])
        strand_l = sc.Strand(substrands=[ss_l])
        design = sc.DNADesign(grid=sc.square, strands=[strand_r, strand_l])
        design.assign_dna(strand_l, 'AAAAC')
        self.assertEqual('GTTTT', strand_r.dna_sequence)

    def test_assign_dna__assign_seq_with_wildcards(self):
        """
        01234
        <---]
        C??AA
        G??TT
        [--->
        """
        ss_bot = sc.Substrand(helix=0, forward=sc.forward, start=0, end=5)
        ss_top = sc.Substrand(helix=0, forward=sc.reverse, start=0, end=5)
        strand_bot = sc.Strand(substrands=[ss_bot])
        strand_top = sc.Strand(substrands=[ss_top])
        strands = [strand_bot, strand_top]
        design = sc.DNADesign(grid=sc.square, strands=strands)
        design.assign_dna(strand_top, 'AA??C')
        self.assertEqual('G??TT', strand_bot.dna_sequence)

    def test_assign_dna__one_strand_assigned_by_complement_from_two_other_strands(self):
        """
          0123     4567
        <-AAAC-] <-GGGA-]
        [-TTTG-----CCCT->
        """
        ss_top_left = sc.Substrand(0, sc.reverse, 0, 4)
        ss_top_right = sc.Substrand(0, sc.reverse, 4, 8)
        ss_bot = sc.Substrand(0, sc.forward, 0, 8)
        st_top_left = sc.Strand([ss_top_left])
        st_top_right = sc.Strand([ss_top_right])
        st_bot = sc.Strand([ss_bot])
        design = sc.DNADesign(strands=[st_bot, st_top_left, st_top_right])
        design.assign_dna(st_top_left, 'CAAA')
        self.assertEqual('TTTG????', st_bot.dna_sequence)
        design.assign_dna(st_top_right, 'AGGG')
        self.assertEqual('TTTGCCCT', st_bot.dna_sequence)

    def test_assign_dna__adapter_assigned_from_scaffold_and_tiles(self):
        """
        XXX: it appears the behavior this tests (which the other tests miss) is assigning DNA to
        tile0 first, then to tile1, and adap is connected to each of them on different helices.
        This means that when tile1 is assigned, we need to ensure when assigning to adap that we
        keep the old information and don't discard it by simply padding the shorter portion of it on
        helix 1 with ?'s, but remember the old DNA sequence.

               01 2345     6789  01
        adap    [-TTTC-----CATT-------+
        scaf <-GT-AAAG-+ <-GTAA--AA-] |
                       |              |
             [-AA-TTTG-+ [-TGCC--GG-> |
                <-AAAC-----ACGG-------+
        """
        scaf0_ss = sc.Substrand(0, sc.reverse, 0, 6)
        scaf1_ss = sc.Substrand(1, sc.forward, 0, 6)
        tile1_ss = sc.Substrand(1, sc.forward, 6, 12)
        tile0_ss = sc.Substrand(0, sc.reverse, 6, 12)
        adap0_ss = sc.Substrand(0, sc.forward, 2, 10)
        adap1_ss = sc.Substrand(1, sc.reverse, 2, 10)
        scaf = sc.Strand([scaf1_ss, scaf0_ss])
        adap = sc.Strand([adap0_ss, adap1_ss])
        tile0 = sc.Strand([tile0_ss])
        tile1 = sc.Strand([tile1_ss])

        design = sc.DNADesign(strands=[scaf, adap, tile0, tile1])

        design.assign_dna(tile0, 'AA AATG')
        self.assertEqual('???? CATT ???? ????'.replace(' ', ''), adap.dna_sequence)

        design.assign_dna(tile1, 'TGCC GG')
        self.assertEqual('???? CATT GGCA ????'.replace(' ', ''), adap.dna_sequence)

        design.assign_dna(scaf, 'AA TTTG GAAA TG')
        self.assertEqual('TTTC CATT GGCA CAAA'.replace(' ', ''), adap.dna_sequence)

    def test_assign_dna__adapter_assigned_from_scaffold_and_tiles_with_deletions(self):
        """
        XXX: it appears the behavior this tests (which the other tests miss) is assigning DNA to
        tile0 first, then to tile1, and adap is connected to each of them on different helices.
        This means that when tile1 is assigned, we need to ensure when assigning to adap that we
        keep the old information and don't discard it by simply padding the shorter portion of it on
        helix 1 with ?'s, but remember the old DNA sequence.

               01 2345     6789  01
                   X         X           deletions
        adap    [-T TC-----CA T-------+
        scaf <-GT-A AG-+ <-GT A--AA-] |
                       |              |
             [-AA-TTTG-+ [-TG C--GG-> |
                <-AAAC-----AC G-------+
                             X           deletions
        """
        scaf0_ss = sc.Substrand(0, sc.reverse, 0, 6)
        scaf1_ss = sc.Substrand(1, sc.forward, 0, 6)
        tile1_ss = sc.Substrand(1, sc.forward, 6, 12)
        tile0_ss = sc.Substrand(0, sc.reverse, 6, 12)
        adap0_ss = sc.Substrand(0, sc.forward, 2, 10)
        adap1_ss = sc.Substrand(1, sc.reverse, 2, 10)
        scaf = sc.Strand([scaf1_ss, scaf0_ss])
        adap = sc.Strand([adap0_ss, adap1_ss])
        tile0 = sc.Strand([tile0_ss])
        tile1 = sc.Strand([tile1_ss])

        design = sc.DNADesign(strands=[scaf, adap, tile0, tile1])
        design.add_deletion(0, 3)
        design.add_deletion(0, 8)
        design.add_deletion(1, 8)

        design.assign_dna(tile0, 'AA ATG')
        self.assertEqual('??? CAT ??? ????'.replace(' ', ''), adap.dna_sequence)

        design.assign_dna(tile1, 'TGC GG')
        self.assertEqual('??? CAT GCA ????'.replace(' ', ''), adap.dna_sequence)

        design.assign_dna(scaf, 'AA TTTG GAA TG')
        self.assertEqual('TTC CAT GCA CAAA'.replace(' ', ''), adap.dna_sequence)

    def test_assign_dna__adapter_assigned_from_scaffold_and_tiles_with_insertions(self):
        """
        XXX: it appears the behavior this tests (which the other tests miss) is assigning DNA to
        tile0 first, then to tile1, and adap is connected to each of them on different helices.
        This means that when tile1 is assigned, we need to ensure when assigning to adap that we
        keep the old information and don't discard it by simply padding the shorter portion of it on
        helix 1 with ?'s, but remember the old DNA sequence.

               01 2345     678I9  01
                              I          insertions
        adap    [-TTTC-----CATTT-------+
        scaf <-GT-AAAG-+ <-GTAAA--AA-] |
                       |              |
             [-AA-TTTG-+ [-TGCCC--GG-> |
                <-AAAC-----ACGGG-------+
                              I          insertions
        """
        scaf0_ss = sc.Substrand(0, sc.reverse, 0, 6)
        scaf1_ss = sc.Substrand(1, sc.forward, 0, 6)
        tile1_ss = sc.Substrand(1, sc.forward, 6, 12)
        tile0_ss = sc.Substrand(0, sc.reverse, 6, 12)
        adap0_ss = sc.Substrand(0, sc.forward, 2, 10)
        adap1_ss = sc.Substrand(1, sc.reverse, 2, 10)
        scaf = sc.Strand([scaf1_ss, scaf0_ss])
        adap = sc.Strand([adap0_ss, adap1_ss])
        tile0 = sc.Strand([tile0_ss])
        tile1 = sc.Strand([tile1_ss])

        design = sc.DNADesign(strands=[scaf, adap, tile0, tile1])
        design.add_insertion(0, 8, 1)
        design.add_insertion(1, 8, 1)

        design.assign_dna(tile0, 'AA AAATG')
        self.assertEqual('???? CATTT ????? ????'.replace(' ', ''), adap.dna_sequence)

        design.assign_dna(tile1, 'TGCCC GG')
        self.assertEqual('???? CATTT GGGCA ????'.replace(' ', ''), adap.dna_sequence)

        design.assign_dna(scaf, 'AA TTTG GAAA TG')
        self.assertEqual('TTTC CATTT GGGCA CAAA'.replace(' ', ''), adap.dna_sequence)

    def test_assign_dna__dna_sequence_shorter_than_complementary_strand_right_strand_longer(self):
        """
        <---]
        CAAAA
        GTTTT?????
        [-------->
        """
        ss_long = sc.Substrand(helix=0, forward=sc.forward, start=0, end=10)
        ss_short = sc.Substrand(helix=0, forward=sc.reverse, start=0, end=5)
        strand_long = sc.Strand(substrands=[ss_long])
        strand_short = sc.Strand(substrands=[ss_short])
        strands = [strand_long, strand_short]
        design = sc.DNADesign(grid=sc.square, strands=strands)
        design.assign_dna(strand_short, 'AAAAC')
        self.assertEqual('GTTTT?????', strand_long.dna_sequence)

    def test_assign_dna__dna_sequence_shorter_than_complementary_strand_left_strand_longer(self):
        """
        [--->
        AAAAC
        TTTTG?????
        <--------]
        """
        ss_long = sc.Substrand(helix=0, forward=sc.reverse, start=0, end=10)
        ss_short = sc.Substrand(helix=0, forward=sc.forward, start=0, end=5)
        strand_long = sc.Strand(substrands=[ss_long])
        strand_short = sc.Strand(substrands=[ss_short])
        strands = [strand_long, strand_short]
        design = sc.DNADesign(grid=sc.square, strands=strands)
        design.assign_dna(strand_short, 'AAAAC')
        self.assertEqual('?????GTTTT', strand_long.dna_sequence)

    def test_assign_dna__dna_sequence_with_uncomplemented_substrand_on_different_helix(self):
        """
        <---]
        CAAAA
        GTTTT?????
        [--------+
                 |
               <-+
               ???
        """
        ss_long = sc.Substrand(helix=0, forward=sc.forward, start=0, end=10)
        ss_long_h1 = sc.Substrand(helix=0, forward=sc.reverse, start=7, end=10)
        ss_short = sc.Substrand(helix=0, forward=sc.reverse, start=0, end=5)
        strand_long = sc.Strand(substrands=[ss_long, ss_long_h1])
        strand_short = sc.Strand(substrands=[ss_short])
        strands = [strand_long, strand_short]
        design = sc.DNADesign(grid=sc.square, strands=strands)
        design.assign_dna(strand_short, 'AAAAC')
        self.assertEqual('GTTTT????????', strand_long.dna_sequence)

    def test_assign_dna__dna_sequence_with_uncomplemented_substrand_on_different_helix_wildcards_both_ends(
            self):
        """
             <---]
             CAAAA
        ?????GTTTT
        [--------+
                 |
               <-+
               ???
        """
        ss_long_h0 = sc.Substrand(helix=0, forward=sc.forward, start=0, end=10)
        ss_long_h1 = sc.Substrand(helix=1, forward=sc.reverse, start=7, end=10)
        ss_short_h0 = sc.Substrand(helix=0, forward=sc.reverse, start=5, end=10)
        strand_long = sc.Strand(substrands=[ss_long_h0, ss_long_h1])
        strand_short = sc.Strand(substrands=[ss_short_h0])
        strands = [strand_long, strand_short]
        design = sc.DNADesign(grid=sc.square, strands=strands)
        design.assign_dna(strand_short, 'AAAAC')
        self.assertEqual('?????GTTTT???', strand_long.dna_sequence)

    def test_assign_dna__one_helix_with_one_bottom_strand_and_three_top_strands(self):
        """
         012   345   678
        -TTT> -GGG> -CCC>
        <AAA---CCC---GGG-
         876   543   210
        """
        ss_bot = sc.Substrand(helix=0, forward=sc.reverse, start=0, end=9)
        ss_top1 = sc.Substrand(helix=0, forward=sc.forward, start=0, end=3)
        ss_top2 = sc.Substrand(helix=0, forward=sc.forward, start=3, end=6)
        ss_top3 = sc.Substrand(helix=0, forward=sc.forward, start=6, end=9)
        strand_bot = sc.Strand(substrands=[ss_bot])
        strand_top1 = sc.Strand(substrands=[ss_top1])
        strand_top2 = sc.Strand(substrands=[ss_top2])
        strand_top3 = sc.Strand(substrands=[ss_top3])
        strands = [strand_bot, strand_top1, strand_top2, strand_top3]
        design = sc.DNADesign(grid=sc.square, strands=strands)
        design.assign_dna(strand_bot, 'AAACCCGGG')
        self.assertEqual('CCC', strand_top1.dna_sequence)
        self.assertEqual('GGG', strand_top2.dna_sequence)
        self.assertEqual('TTT', strand_top3.dna_sequence)

    def test_assign_dna__two_helices_with_multiple_substrand_intersections(self):
        """
                012    345   678    901
        M13   [-ACC----TAA---GAA----AAC---+
              +-TGG-]<-ATT-+ CTT----TTG-+ |
              |            | |          | |
              +-GAT----TTC-+ ATG->[-AGT-+ |
              <-CTA----AAG---TAC----TCA---+
        """
        scaf0_ss = sc.Substrand(helix=0, forward=sc.forward, start=0, end=12)
        scaf1_ss = sc.Substrand(helix=1, forward=sc.reverse, start=0, end=12)
        scaf = sc.Strand(substrands=[scaf0_ss, scaf1_ss])

        first_stap0_left_ss = sc.Substrand(helix=0, forward=sc.reverse, start=0, end=3)
        first_stap1_ss = sc.Substrand(helix=1, forward=sc.forward, start=0, end=6)
        first_stap0_right_ss = sc.Substrand(helix=0, forward=sc.reverse, start=3, end=6)
        first_stap = sc.Strand(substrands=[first_stap0_left_ss, first_stap1_ss, first_stap0_right_ss])

        second_stap1_right_ss = sc.Substrand(helix=1, forward=sc.forward, start=9, end=12)
        second_stap0_ss = sc.Substrand(helix=0, forward=sc.reverse, start=6, end=12)
        second_stap1_left_ss = sc.Substrand(helix=1, forward=sc.forward, start=6, end=9)
        second_stap = sc.Strand(substrands=[second_stap1_right_ss, second_stap0_ss, second_stap1_left_ss])

        strands = [scaf, first_stap, second_stap]
        design = sc.DNADesign(grid=sc.square, strands=strands)
        design.assign_dna(scaf, 'ACC TAA GAA AAC ACT CAT GAA ATC'.replace(' ', ''))
        self.assertEqual('GGT GAT TTC TTA'.replace(' ', ''), first_stap.dna_sequence)
        self.assertEqual('AGT GTT TTC ATG'.replace(' ', ''), second_stap.dna_sequence)

    def test_assign_dna__upper_left_edge_staple_of_16H_origami_rectangle(self):
        """
        staple <ACATAAGAAAACGGAG--+
        M13   +-TGTATTCTTTTGCCTC> |
              |                   |
              +-GATTTTGTGAGTAGAA- |
               -CTAAAACACTCATCTT--+
        """
        scaf0_ss = sc.Substrand(helix=0, forward=sc.forward, start=0, end=16)
        scaf1_ss = sc.Substrand(helix=1, forward=sc.reverse, start=0, end=16)
        stap0_ss = sc.Substrand(helix=0, forward=sc.reverse, start=0, end=16)
        stap1_ss = sc.Substrand(helix=1, forward=sc.forward, start=0, end=16)
        scaf = sc.Strand(substrands=[scaf1_ss, scaf0_ss])
        stap = sc.Strand(substrands=[stap1_ss, stap0_ss])
        strands = [scaf, stap]
        design = sc.DNADesign(grid=sc.square, strands=strands)

        seq_m13_upper_left = 'AAGATGAGTGTTTTAGTGTATTCTTTTGCCTC'
        design.assign_dna(scaf, seq_m13_upper_left)
        expected_seq_stap_upperleft = 'CTAAAACACTCATCTTGAGGCAAAAGAATACA'
        self.assertEqual(expected_seq_stap_upperleft, stap.dna_sequence)

    def test_assign_dna__2helix_with_deletions(self):
        """
        scaf index: 2     3  4     5
        offset:     0 D1  2  3 D4  5
                    +     -  -     +
                   /C     A  T     C\
                  | G     T  A     G |
        helix 0   | <     +  +     ] |
                  |       |  |       |
        helix 1   | [     +  +     > |
                  | T     T  A     C |
                   \A     A  T     G/
                    +     ]  <     +
        offset:     0 D1  2  3 D4  5
        scaf index: 1     0  7     6
        """
        width = 6
        width_h = width // 2
        helices = [sc.Helix(max_offset=width), sc.Helix(max_offset=width)]
        stap_left_ss1 = sc.Substrand(1, sc.forward, 0, width_h)
        stap_left_ss0 = sc.Substrand(0, sc.reverse, 0, width_h)
        stap_right_ss0 = sc.Substrand(0, sc.reverse, width_h, width)
        stap_right_ss1 = sc.Substrand(1, sc.forward, width_h, width)
        scaf_ss1_left = sc.Substrand(1, sc.reverse, 0, width_h)
        scaf_ss0 = sc.Substrand(0, sc.forward, 0, width)
        scaf_ss1_right = sc.Substrand(1, sc.reverse, width_h, width)
        stap_left = sc.Strand([stap_left_ss1, stap_left_ss0])
        stap_right = sc.Strand([stap_right_ss0, stap_right_ss1])
        scaf = sc.Strand([scaf_ss1_left, scaf_ss0, scaf_ss1_right], color=sc.default_scaffold_color)
        strands = [stap_left, stap_right, scaf]
        design = sc.DNADesign(helices=helices, strands=strands, grid=sc.square)
        design.add_deletion(helix=0, offset=1)
        design.add_deletion(helix=0, offset=4)
        design.add_deletion(helix=1, offset=1)
        design.add_deletion(helix=1, offset=4)
        design.assign_dna(scaf, 'AACATCGT')
        self.assertEqual("AACATCGT", scaf.dna_sequence)
        self.assertEqual("TTTG", stap_left.dna_sequence)
        self.assertEqual("GAAC", stap_right.dna_sequence)

    def test_assign_dna__wildcards_simple(self):
        """
         012   345   678
        -TTC> -GGA> -CCT>
        <AAG---CCT---GGA-
         876   543   210
        """
        ss_bot = sc.Substrand(helix=0, forward=sc.reverse, start=0, end=9)
        ss_top1 = sc.Substrand(helix=0, forward=sc.forward, start=0, end=3)
        ss_top2 = sc.Substrand(helix=0, forward=sc.forward, start=3, end=6)
        ss_top3 = sc.Substrand(helix=0, forward=sc.forward, start=6, end=9)
        strand_bot = sc.Strand(substrands=[ss_bot])
        strand_top1 = sc.Strand(substrands=[ss_top1])
        strand_top2 = sc.Strand(substrands=[ss_top2])
        strand_top3 = sc.Strand(substrands=[ss_top3])
        strands = [strand_bot, strand_top1, strand_top2, strand_top3]
        design = sc.DNADesign(grid=sc.square, strands=strands)

        design.assign_dna(strand_top1, 'TTC')
        self.assertEqual('??????GAA', strand_bot.dna_sequence)

        design.assign_dna(strand_top3, 'CCT')
        self.assertEqual('AGG???GAA', strand_bot.dna_sequence)

        design.assign_dna(strand_top2, 'GGA')
        self.assertEqual('AGGTCCGAA', strand_bot.dna_sequence)

    def test_assign_dna__wildcards_multiple_overlaps(self):
        """
         012   345   678   901   234   567
              +---------------+
              |               |
              |         +-----|-------+
              |         |     |       |
        -ACG> +TTC> -GGA+ -AAC+ -TGC> +TTG>
        <TGC---AAG---CCT---TTG---ACG---AAC---???-
         098   765   432   109   876   543   210
        """
        ss_bot = sc.Substrand(helix=0, forward=sc.reverse, start=0, end=21)
        ss_top0 = sc.Substrand(helix=0, forward=sc.forward, start=0, end=3)
        ss_top3 = sc.Substrand(helix=0, forward=sc.forward, start=3, end=6)
        ss_top6 = sc.Substrand(helix=0, forward=sc.forward, start=6, end=9)
        ss_top9 = sc.Substrand(helix=0, forward=sc.forward, start=9, end=12)
        ss_top12 = sc.Substrand(helix=0, forward=sc.forward, start=12, end=15)
        ss_top15 = sc.Substrand(helix=0, forward=sc.forward, start=15, end=18)
        strand_bot = sc.Strand(substrands=[ss_bot])
        strand_top_small0 = sc.Strand(substrands=[ss_top0])
        strand_top_small12 = sc.Strand(substrands=[ss_top12])
        strand_top_big9 = sc.Strand(substrands=[ss_top9, ss_top3])
        strand_top_big6 = sc.Strand(substrands=[ss_top6, ss_top15])
        strands = [strand_bot, strand_top_small0, strand_top_small12, strand_top_big9, strand_top_big6]
        design = sc.DNADesign(grid=sc.square, strands=strands)

        design.assign_dna(strand_top_big9, 'AACTTC')
        self.assertEqual('?????????GTT???GAA???', strand_bot.dna_sequence)

        design.assign_dna(strand_top_small12, 'TGC')
        self.assertEqual('??????GCAGTT???GAA???', strand_bot.dna_sequence)

        design.assign_dna(strand_top_small0, 'ACG')
        self.assertEqual('??????GCAGTT???GAACGT', strand_bot.dna_sequence)

        design.assign_dna(strand_top_big6, 'GGATTG')
        self.assertEqual('???CAAGCAGTTTCCGAACGT', strand_bot.dna_sequence)


TEST_OFFSETS_AT_DELETION_INSERTIONS = False


class TestSubstrandDNASequenceIn(unittest.TestCase):

    def test_dna_sequence_in__right_then_left(self):
        ss0 = sc.Substrand(0, sc.forward, 0, 10)
        ss1 = sc.Substrand(1, sc.reverse, 0, 10)
        strand = sc.Strand([ss0, ss1])
        strand.dna_sequence = "AAAACCCCGGGGTTTTACGT"
        # offset: 0  1  2  3  4  5  6  7  8  9
        # index:  0  1  2  3  4  5  6  7  8  9
        #         A  A  A  A  C  C  C  C  G  G
        # helix 0 [  -  -  -  -  -  -  -  -  +
        #                                    |
        # helix 1 <  -  -  -  -  -  -  -  -  +
        #         T  G  C  A  T  T  T  T  G  G
        # offset: 0  1  2  3  4  5  6  7  8  9
        # index: 19 18 17 16 15 14 13 12 11 10
        self.assertEqual("A", ss0.dna_sequence_in(0, 0))
        self.assertEqual("AA", ss0.dna_sequence_in(0, 1))
        self.assertEqual("AAA", ss0.dna_sequence_in(0, 2))
        self.assertEqual("AAAA", ss0.dna_sequence_in(0, 3))
        self.assertEqual("AAAAC", ss0.dna_sequence_in(0, 4))
        self.assertEqual("AAAACC", ss0.dna_sequence_in(0, 5))
        self.assertEqual("AAAACCC", ss0.dna_sequence_in(0, 6))
        self.assertEqual("AAAACCCC", ss0.dna_sequence_in(0, 7))
        self.assertEqual("AAAACCCCG", ss0.dna_sequence_in(0, 8))
        self.assertEqual("AAAACCCCGG", ss0.dna_sequence_in(0, 9))
        #
        self.assertEqual("G", ss1.dna_sequence_in(9, 9))
        self.assertEqual("GG", ss1.dna_sequence_in(8, 9))
        self.assertEqual("GGT", ss1.dna_sequence_in(7, 9))
        self.assertEqual("GGTT", ss1.dna_sequence_in(6, 9))
        self.assertEqual("GGTTT", ss1.dna_sequence_in(5, 9))
        self.assertEqual("GGTTTT", ss1.dna_sequence_in(4, 9))
        self.assertEqual("GGTTTTA", ss1.dna_sequence_in(3, 9))
        self.assertEqual("GGTTTTAC", ss1.dna_sequence_in(2, 9))
        self.assertEqual("GGTTTTACG", ss1.dna_sequence_in(1, 9))
        self.assertEqual("GGTTTTACGT", ss1.dna_sequence_in(0, 9))

    def test_dna_sequence_in__right_then_left_deletions(self):
        ss0 = sc.Substrand(0, sc.forward, 0, 10, deletions=[2, 5, 6])
        ss1 = sc.Substrand(1, sc.reverse, 0, 10, deletions=[2, 6, 7])
        strand = sc.Strand([ss0, ss1])
        strand.dna_sequence = "AAACCGGGGTTAGT"
        # offset: 0  1 D2  3  4 D5 D6  7  8  9
        # index:  0  1     2  3        4  5  6
        #         A  A     A  C        C  G  G
        # helix 0 [  -  -  -  -  -  -  -  -  +
        #                                    |
        # helix 1 <  -  -  -  -  -  -  -  -  +
        #         T  G     A  T  T        G  G
        # offset: 0  1 D2  3  4  5 D6 D7  8  9
        # index: 13 12    11 10  9        9  7
        self.assertEqual("A", ss0.dna_sequence_in(0, 0))
        self.assertEqual("AA", ss0.dna_sequence_in(0, 1))
        self.assertEqual("AA", ss0.dna_sequence_in(0, 2))
        self.assertEqual("AAA", ss0.dna_sequence_in(0, 3))
        self.assertEqual("AAAC", ss0.dna_sequence_in(0, 4))
        self.assertEqual("AAAC", ss0.dna_sequence_in(0, 5))
        self.assertEqual("AAAC", ss0.dna_sequence_in(0, 6))
        self.assertEqual("AAACC", ss0.dna_sequence_in(0, 7))
        self.assertEqual("AAACCG", ss0.dna_sequence_in(0, 8))
        self.assertEqual("AAACCGG", ss0.dna_sequence_in(0, 9))
        #
        self.assertEqual("G", ss1.dna_sequence_in(9, 9))
        self.assertEqual("GG", ss1.dna_sequence_in(8, 9))
        self.assertEqual("GG", ss1.dna_sequence_in(7, 9))
        self.assertEqual("GG", ss1.dna_sequence_in(6, 9))
        self.assertEqual("GGT", ss1.dna_sequence_in(5, 9))
        self.assertEqual("GGTT", ss1.dna_sequence_in(4, 9))
        self.assertEqual("GGTTA", ss1.dna_sequence_in(3, 9))
        self.assertEqual("GGTTA", ss1.dna_sequence_in(2, 9))
        self.assertEqual("GGTTAG", ss1.dna_sequence_in(1, 9))
        self.assertEqual("GGTTAGT", ss1.dna_sequence_in(0, 9))

        # if TEST_OFFSETS_AT_DELETION_INSERTIONS:
        #     self.assertEqual("AA", ss0.dna_sequence_in(0, 3))
        #     self.assertEqual("AAACC", ss0.dna_sequence_in(0, 7))
        #     self.assertEqual("GGT", ss1.dna_sequence_in(6, 10))
        #     self.assertEqual("GGTTTA", ss1.dna_sequence_in(2, 10))

    def test_dna_sequence_in__right_then_left_insertions(self):
        ss0 = sc.Substrand(0, sc.forward, 0, 10, insertions=[(2, 1), (6, 2)])
        ss1 = sc.Substrand(1, sc.reverse, 0, 10, insertions=[(2, 1), (6, 2)])
        strand = sc.Strand([ss0, ss1])
        strand.dna_sequence = "AAAACCCCGGGGTTTTACGTACGTAC"
        # offset: 0  1  2  I  3  4  5  6  I  I  7  8  9
        # index:  0  1  2  3  4  5  6  7  8  9 10 11 12
        #         A  A  A  A  C  C  C  C  G  G  G  G  T
        # helix 0 [  -  -  -  -  -  -  -  -  -  -  -  +
        #                                             |
        # helix 1 <  -  -  -  -  -  -  -  -  -  -  -  +
        #         C  A  T  G  C  A  T  G  C  A  T  T  T
        # offset: 0  1  2  I  3  4  5  6  I  I  7  8  9
        # index: 25 24 23 22 21 20 19 18 17 16 15 14 13
        self.assertEqual("A", ss0.dna_sequence_in(0, 0))
        self.assertEqual("AA", ss0.dna_sequence_in(0, 1))
        self.assertEqual("AAAA", ss0.dna_sequence_in(0, 2))
        self.assertEqual("AAAAC", ss0.dna_sequence_in(0, 3))
        self.assertEqual("AAAACC", ss0.dna_sequence_in(0, 4))
        self.assertEqual("AAAACCC", ss0.dna_sequence_in(0, 5))
        self.assertEqual("AAAACCCCGG", ss0.dna_sequence_in(0, 6))
        self.assertEqual("AAAACCCCGGG", ss0.dna_sequence_in(0, 7))
        self.assertEqual("AAAACCCCGGGG", ss0.dna_sequence_in(0, 8))
        self.assertEqual("AAAACCCCGGGGT", ss0.dna_sequence_in(0, 9))
        #
        self.assertEqual("T", ss1.dna_sequence_in(9, 9))
        self.assertEqual("TT", ss1.dna_sequence_in(8, 9))
        self.assertEqual("TTT", ss1.dna_sequence_in(7, 9))
        self.assertEqual("TTTACG", ss1.dna_sequence_in(6, 9))
        self.assertEqual("TTTACGT", ss1.dna_sequence_in(5, 9))
        self.assertEqual("TTTACGTA", ss1.dna_sequence_in(4, 9))
        self.assertEqual("TTTACGTAC", ss1.dna_sequence_in(3, 9))
        self.assertEqual("TTTACGTACGT", ss1.dna_sequence_in(2, 9))
        self.assertEqual("TTTACGTACGTA", ss1.dna_sequence_in(1, 9))
        self.assertEqual("TTTACGTACGTAC", ss1.dna_sequence_in(0, 9))

        # if TEST_OFFSETS_AT_DELETION_INSERTIONS:
        #     self.assertEqual("AAAA", ss0.dna_sequence_in(0, 3))
        #     self.assertEqual("AAAACCCCGG", ss0.dna_sequence_in(0, 7))
        #     self.assertEqual("TTTACG", ss1.dna_sequence_in(6, 10))
        #     self.assertEqual("TTTACGTACGT", ss1.dna_sequence_in(2, 10))

    def test_dna_sequence_in__right_then_left_deletions_and_insertions(self):
        ss0 = sc.Substrand(0, sc.forward, 0, 10, deletions=[4], insertions=[(2, 1), (6, 2)])
        ss1 = sc.Substrand(1, sc.reverse, 0, 10, deletions=[4], insertions=[(2, 1), (6, 2)])
        strand = sc.Strand([ss0, ss1])
        strand.dna_sequence = "AAAACCCCGGGGTTTTACGTACGTAC"
        # offset: 0  1  2  I  3 D4  5  6  I  I  7  8  9
        # index:  0  1  2  3  4     5  6  7  8  9 10 11
        #         A  A  A  A  C     C  C  C  G  G  G  G
        # helix 0 [  -  -  -  -  -  -  -  -  -  -  -  +
        #                                             |
        # helix 1 <  -  -  -  -  -  -  -  -  -  -  -  +
        #         T  G  C  A  T     G  C  A  T  T  T  T
        # offset: 0  1  2  I  3 D4  5  6  I  I  7  8  9
        # index: 23 22 21 20 19    18 17 16 15 14 13 12
        self.assertEqual("AA", ss0.dna_sequence_in(2, 2))
        self.assertEqual("CCG", ss0.dna_sequence_in(6, 6))
        self.assertEqual("TAC", ss1.dna_sequence_in(6, 6))
        self.assertEqual("AC", ss1.dna_sequence_in(2, 2))
        #
        self.assertEqual("A           ".strip(), ss0.dna_sequence_in(0, 0))
        self.assertEqual("AA          ".strip(), ss0.dna_sequence_in(0, 1))
        self.assertEqual("AAAA        ".strip(), ss0.dna_sequence_in(0, 2))
        self.assertEqual("AAAAC       ".strip(), ss0.dna_sequence_in(0, 3))
        self.assertEqual("AAAAC       ".strip(), ss0.dna_sequence_in(0, 4))
        self.assertEqual("AAAACC      ".strip(), ss0.dna_sequence_in(0, 5))
        self.assertEqual("AAAACCCCG   ".strip(), ss0.dna_sequence_in(0, 6))
        self.assertEqual("AAAACCCCGG  ".strip(), ss0.dna_sequence_in(0, 7))
        self.assertEqual("AAAACCCCGGG ".strip(), ss0.dna_sequence_in(0, 8))
        self.assertEqual("AAAACCCCGGGG".strip(), ss0.dna_sequence_in(0, 9))
        self.assertEqual(" AAACCCCGGGG".strip(), ss0.dna_sequence_in(1, 9))
        self.assertEqual("  AACCCCGGGG".strip(), ss0.dna_sequence_in(2, 9))
        self.assertEqual("    CCCCGGGG".strip(), ss0.dna_sequence_in(3, 9))
        self.assertEqual("     CCCGGGG".strip(), ss0.dna_sequence_in(4, 9))
        self.assertEqual("     CCCGGGG".strip(), ss0.dna_sequence_in(5, 9))
        self.assertEqual("      CCGGGG".strip(), ss0.dna_sequence_in(6, 9))
        self.assertEqual("         GGG".strip(), ss0.dna_sequence_in(7, 9))
        self.assertEqual("          GG".strip(), ss0.dna_sequence_in(8, 9))
        self.assertEqual("           G".strip(), ss0.dna_sequence_in(9, 9))
        # strip() below is so auto-formatting preserves nice lineup of characters
        self.assertEqual("T           ".strip(), ss1.dna_sequence_in(9, 9))
        self.assertEqual("TT          ".strip(), ss1.dna_sequence_in(8, 9))
        self.assertEqual("TTT         ".strip(), ss1.dna_sequence_in(7, 9))
        self.assertEqual("TTTTAC      ".strip(), ss1.dna_sequence_in(6, 9))
        self.assertEqual("TTTTACG     ".strip(), ss1.dna_sequence_in(5, 9))
        self.assertEqual("TTTTACG     ".strip(), ss1.dna_sequence_in(4, 9))
        self.assertEqual("TTTTACGT    ".strip(), ss1.dna_sequence_in(3, 9))
        self.assertEqual("TTTTACGTAC  ".strip(), ss1.dna_sequence_in(2, 9))
        self.assertEqual("TTTTACGTACG ".strip(), ss1.dna_sequence_in(1, 9))
        self.assertEqual("TTTTACGTACGT".strip(), ss1.dna_sequence_in(0, 9))
        self.assertEqual(" TTTACGTACGT".strip(), ss1.dna_sequence_in(0, 8))
        self.assertEqual("  TTACGTACGT".strip(), ss1.dna_sequence_in(0, 7))
        self.assertEqual("   TACGTACGT".strip(), ss1.dna_sequence_in(0, 6))
        self.assertEqual("      GTACGT".strip(), ss1.dna_sequence_in(0, 5))
        self.assertEqual("       TACGT".strip(), ss1.dna_sequence_in(0, 4))
        self.assertEqual("       TACGT".strip(), ss1.dna_sequence_in(0, 3))
        self.assertEqual("        ACGT".strip(), ss1.dna_sequence_in(0, 2))
        self.assertEqual("          GT".strip(), ss1.dna_sequence_in(0, 1))
        self.assertEqual("           T".strip(), ss1.dna_sequence_in(0, 0))

        # if TEST_OFFSETS_AT_DELETION_INSERTIONS:
        #     self.assertEqual("AAAA", ss0.dna_sequence_in(0, 3))
        #     self.assertEqual("AAAAC", ss0.dna_sequence_in(0, 5))
        #     self.assertEqual("AAAACCCCGG", ss0.dna_sequence_in(0, 7))
        #     self.assertEqual("TTTACG", ss1.dna_sequence_in(6, 10))
        #     self.assertEqual("TTTTACG", ss1.dna_sequence_in(4, 10))
        #     self.assertEqual("TTTACGTACGT", ss1.dna_sequence_in(2, 10))
