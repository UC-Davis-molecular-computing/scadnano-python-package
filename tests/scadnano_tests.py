import os
import unittest
import re
import json
import math
from typing import Iterable, Union, Dict, Any

import xlrd

import scadnano as sc
import scadnano.origami_rectangle as rect
import scadnano.modifications as mod


def strand_matching(strands: Iterable[sc.Strand], helix: int, forward: bool, start: int,
                    end: int) -> sc.Strand:
    """
    Finds strand whose first bound domain matches the given parameters.
    """
    return next(s for s in strands if
                s.first_bound_domain().helix == helix and
                s.first_bound_domain().forward == forward and
                s.first_bound_domain().start == start and
                s.first_bound_domain().end == end)


def remove_whitespace(sequence: str) -> str:
    sequence = re.sub(r'\s*', '', sequence)
    return sequence


class TestCreateStrandChainedMethods(unittest.TestCase):
    # tests methods for creating strands using chained method notation as in this issue:
    # https://github.com/UC-Davis-molecular-computing/scadnano-python-package/issues/85

    def setUp(self) -> None:
        helices = [sc.Helix(max_offset=100) for _ in range(6)]
        self.design_6helix: sc.Design = sc.Design(helices=helices, strands=[], grid=sc.square)

    def test_strand__loopouts_with_labels(self) -> None:
        design = self.design_6helix
        sb = design.strand(0, 0)
        sb.to(10)
        sb.loopout(1, 8)
        sb.with_domain_label('loop0')
        sb.to(5)
        sb.with_domain_label('dom1')
        sb.cross(2)
        sb.to(10)
        sb.with_domain_label('dom2')
        sb.loopout(3, 12)
        sb.with_domain_label('loop1')
        sb.to(5)
        expected_strand = sc.Strand([
            sc.Domain(0, True, 0, 10),
            sc.Loopout(8, label='loop0'),
            sc.Domain(1, False, 5, 10, label='dom1'),
            sc.Domain(2, True, 5, 10, label='dom2'),
            sc.Loopout(12, label='loop1'),
            sc.Domain(3, False, 5, 10),
        ])
        self.assertEqual(1, len(design.strands))
        self.assertEqual(expected_strand, design.strands[0])

    def test_strand__loopouts_with_labels_to_json(self) -> None:
        design = self.design_6helix
        sb = design.strand(0, 0)
        sb.to(10)
        sb.loopout(1, 8)
        sb.with_domain_label('loop0')
        sb.to(5)
        sb.with_domain_label('dom1')
        sb.cross(2)
        sb.to(10)
        sb.with_domain_label('dom2')
        sb.loopout(3, 12)
        sb.with_domain_label('loop1')
        sb.to(5)
        design_json_map = design.to_json_serializable(suppress_indent=False)
        design_from_json = sc.Design.from_scadnano_json_map(design_json_map)
        expected_strand = sc.Strand([
            sc.Domain(0, True, 0, 10),
            sc.Loopout(8, label='loop0'),
            sc.Domain(1, False, 5, 10, label='dom1'),
            sc.Domain(2, True, 5, 10, label='dom2'),
            sc.Loopout(12, label='loop1'),
            sc.Domain(3, False, 5, 10),
        ])
        self.assertEqual(1, len(design_from_json.strands))
        self.assertEqual(expected_strand, design_from_json.strands[0])

    def test_strand__0_0_to_10_cross_1_to_5(self) -> None:
        design = self.design_6helix
        sb = design.strand(0, 0)
        sb.to(10)
        sb.cross(1)
        sb.to(5)
        expected_strand = sc.Strand([
            sc.Domain(0, True, 0, 10),
            sc.Domain(1, False, 5, 10),
        ])
        self.assertEqual(1, len(design.strands))
        self.assertEqual(expected_strand, design.strands[0])
        self.assertEqual(1, len(design.helices[0].domains))
        self.assertEqual(1, len(design.helices[1].domains))
        self.assertEqual(0, len(design.helices[2].domains))
        self.assertEqual(0, len(design.helices[3].domains))
        self.assertEqual(0, len(design.helices[4].domains))
        self.assertEqual(0, len(design.helices[5].domains))

    def test_strand__0_0_to_10_cross_1_to_5__reverse(self) -> None:
        design = self.design_6helix
        design.strand(1, 5).to(10).cross(0).to(0)
        expected_strand = sc.Strand([
            sc.Domain(1, True, 5, 10),
            sc.Domain(0, False, 0, 10),
        ])
        self.assertEqual(1, len(design.strands))
        self.assertEqual(expected_strand, design.strands[0])
        self.assertEqual(1, len(design.helices[0].domains))
        self.assertEqual(1, len(design.helices[1].domains))
        self.assertEqual(0, len(design.helices[2].domains))
        self.assertEqual(0, len(design.helices[3].domains))
        self.assertEqual(0, len(design.helices[4].domains))
        self.assertEqual(0, len(design.helices[5].domains))

    def test_strand__h0_off0_to_off10_cross_h1_to_off5_loopout_length3_h2_to_off15(self) -> None:
        design = self.design_6helix
        sb = design.strand(0, 0)
        sb.to(10)
        sb.cross(1)
        sb.to(5)
        sb.loopout(2, 3)
        sb.to(15)
        expected_strand = sc.Strand([
            sc.Domain(0, True, 0, 10),
            sc.Domain(1, False, 5, 10),
            sc.Loopout(3),
            sc.Domain(2, True, 5, 15),
        ])
        self.assertEqual(1, len(design.strands))
        self.assertEqual(expected_strand, design.strands[0])
        self.assertEqual(1, len(design.helices[0].domains))
        self.assertEqual(1, len(design.helices[1].domains))
        self.assertEqual(1, len(design.helices[2].domains))
        self.assertEqual(0, len(design.helices[3].domains))
        self.assertEqual(0, len(design.helices[4].domains))
        self.assertEqual(0, len(design.helices[5].domains))

    def test_strand__two_forward_paranemic_crossovers(self) -> None:
        design = self.design_6helix
        design.strand(0, 0).to(10).cross(1).to(15).cross(2).to(20)
        expected_strand = sc.Strand([
            sc.Domain(0, True, 0, 10),
            sc.Domain(1, True, 10, 15),
            sc.Domain(2, True, 15, 20),
        ])
        self.assertEqual(1, len(design.strands))
        self.assertEqual(expected_strand, design.strands[0])
        self.assertEqual(1, len(design.helices[0].domains))
        self.assertEqual(1, len(design.helices[1].domains))
        self.assertEqual(1, len(design.helices[2].domains))
        self.assertEqual(0, len(design.helices[3].domains))
        self.assertEqual(0, len(design.helices[4].domains))
        self.assertEqual(0, len(design.helices[5].domains))

    def test_strand__two_reverse_paranemic_crossovers(self) -> None:
        design = self.design_6helix
        design.strand(0, 20).to(10).cross(1).to(5).cross(2).to(0)
        expected_strand = sc.Strand([
            sc.Domain(0, False, 10, 20),
            sc.Domain(1, False, 5, 10),
            sc.Domain(2, False, 0, 5),
        ])
        self.assertEqual(1, len(design.strands))
        self.assertEqual(expected_strand, design.strands[0])
        self.assertEqual(1, len(design.helices[0].domains))
        self.assertEqual(1, len(design.helices[1].domains))
        self.assertEqual(1, len(design.helices[2].domains))
        self.assertEqual(0, len(design.helices[3].domains))
        self.assertEqual(0, len(design.helices[4].domains))
        self.assertEqual(0, len(design.helices[5].domains))

    def test_strand__multiple_strands(self) -> None:
        design = self.design_6helix
        design.strand(0, 0).to(10).cross(1).to(0)
        design.strand(0, 20).to(10).cross(1).to(20)
        expected_strand0 = sc.Strand([
            sc.Domain(0, True, 0, 10),
            sc.Domain(1, False, 0, 10),
        ])
        expected_strand1 = sc.Strand([
            sc.Domain(0, False, 10, 20),
            sc.Domain(1, True, 10, 20),
        ])
        self.assertEqual(2, len(design.strands))
        self.assertEqual(expected_strand0, design.strands[0])
        self.assertEqual(expected_strand1, design.strands[1])
        self.assertEqual(2, len(design.helices[0].domains))
        self.assertEqual(2, len(design.helices[1].domains))
        self.assertEqual(0, len(design.helices[2].domains))
        self.assertEqual(0, len(design.helices[3].domains))
        self.assertEqual(0, len(design.helices[4].domains))
        self.assertEqual(0, len(design.helices[5].domains))

    def test_strand__multiple_strands_other_order(self) -> None:
        design = self.design_6helix
        design.strand(0, 20).to(10).cross(1).to(20)
        design.strand(0, 0).to(10).cross(1).to(0)
        expected_strand0 = sc.Strand([
            sc.Domain(0, False, 10, 20),
            sc.Domain(1, True, 10, 20),
        ])
        expected_strand1 = sc.Strand([
            sc.Domain(0, True, 0, 10),
            sc.Domain(1, False, 0, 10),
        ])
        self.assertEqual(2, len(design.strands))
        self.assertEqual(expected_strand0, design.strands[0])
        self.assertEqual(expected_strand1, design.strands[1])
        self.assertEqual(2, len(design.helices[0].domains))
        self.assertEqual(2, len(design.helices[1].domains))
        self.assertEqual(0, len(design.helices[2].domains))
        self.assertEqual(0, len(design.helices[3].domains))
        self.assertEqual(0, len(design.helices[4].domains))
        self.assertEqual(0, len(design.helices[5].domains))

    def test_strand__multiple_strands_overlap_no_error(self) -> None:
        design = self.design_6helix
        design.strand(0, 0).to(10).cross(1).to(0) \
            .as_scaffold() \
            .with_modification_internal(5, mod.cy3_int, warn_on_no_dna=False)
        design.strand(0, 10).to(0).cross(1).to(10).with_modification_5p(mod.biotin_5p)
        expected_strand0 = sc.Strand([
            sc.Domain(0, True, 0, 10),
            sc.Domain(1, False, 0, 10),
        ], is_scaffold=True)
        expected_strand1 = sc.Strand([
            sc.Domain(0, False, 0, 10),
            sc.Domain(1, True, 0, 10),
        ])

        expected_strand0.set_modification_internal(5, mod.cy3_int, warn_on_no_dna=False)
        expected_strand1.set_modification_5p(mod.biotin_5p)

        self.assertEqual(2, len(design.strands))

        self.assertEqual(expected_strand0, design.strands[0])
        self.assertEqual(None, design.strands[0].modification_5p)
        self.assertEqual(None, design.strands[0].modification_3p)
        self.assertDictEqual({5: mod.cy3_int}, design.strands[0].modifications_int)

        self.assertEqual(expected_strand1, design.strands[1])
        self.assertEqual(mod.biotin_5p, design.strands[1].modification_5p)
        self.assertEqual(None, design.strands[1].modification_3p)
        self.assertDictEqual({}, design.strands[1].modifications_int)

        self.assertEqual(2, len(design.helices[0].domains))
        self.assertEqual(2, len(design.helices[1].domains))
        self.assertEqual(0, len(design.helices[2].domains))
        self.assertEqual(0, len(design.helices[3].domains))
        self.assertEqual(0, len(design.helices[4].domains))
        self.assertEqual(0, len(design.helices[5].domains))

    def test_strand__multiple_strands_overlap_error(self) -> None:
        design = self.design_6helix
        design.strand(0, 0).to(10).cross(1).to(0)
        with self.assertRaises(sc.IllegalDesignError):
            design.strand(0, 2).to(8)

    def test_strand__call_to_twice_legally(self) -> None:
        design = self.design_6helix
        sb = design.strand(0, 0)
        sb.to(10)
        sb.cross(1)
        sb.to(5)
        sb.to(0)
        expected_strand = sc.Strand([
            sc.Domain(0, True, 0, 10),
            sc.Domain(1, False, 5, 10),
            sc.Domain(1, False, 0, 5),
        ])
        self.assertEqual(1, len(design.strands))
        self.assertEqual(expected_strand, design.strands[0])
        self.assertEqual(1, len(design.helices[0].domains))
        self.assertEqual(2, len(design.helices[1].domains))
        self.assertEqual(0, len(design.helices[2].domains))
        self.assertEqual(0, len(design.helices[3].domains))
        self.assertEqual(0, len(design.helices[4].domains))
        self.assertEqual(0, len(design.helices[5].domains))

    def test_strand__call_update_to_twice_legally(self) -> None:
        design = self.design_6helix
        sb = design.strand(0, 0)
        sb.to(10)
        sb.cross(1)
        sb.update_to(5)
        sb.update_to(0)
        expected_strand = sc.Strand([
            sc.Domain(0, True, 0, 10),
            sc.Domain(1, False, 0, 10),
        ])
        self.assertEqual(1, len(design.strands))
        self.assertEqual(expected_strand, design.strands[0])
        self.assertEqual(1, len(design.helices[0].domains))
        self.assertEqual(1, len(design.helices[1].domains))
        self.assertEqual(0, len(design.helices[2].domains))
        self.assertEqual(0, len(design.helices[3].domains))
        self.assertEqual(0, len(design.helices[4].domains))
        self.assertEqual(0, len(design.helices[5].domains))

    def test_strand__call_to_then_update_to_legally(self) -> None:
        design = self.design_6helix
        sb = design.strand(0, 0)
        sb.to(10)
        sb.cross(1)
        sb.to(5)
        sb.update_to(0)
        expected_strand = sc.Strand([
            sc.Domain(0, True, 0, 10),
            sc.Domain(1, False, 0, 10),
        ])
        self.assertEqual(1, len(design.strands))
        self.assertEqual(expected_strand, design.strands[0])
        self.assertEqual(1, len(design.helices[0].domains))
        self.assertEqual(1, len(design.helices[1].domains))
        self.assertEqual(0, len(design.helices[2].domains))
        self.assertEqual(0, len(design.helices[3].domains))
        self.assertEqual(0, len(design.helices[4].domains))
        self.assertEqual(0, len(design.helices[5].domains))

    def test_strand__call_to_twice_increase_decrease_forward(self) -> None:
        design = self.design_6helix
        sb = design.strand(0, 0)
        sb.to(10)
        with self.assertRaises(sc.IllegalDesignError):
            sb.to(5)

    def test_strand__call_to_twice_decrease_increase_reverse(self) -> None:
        design = self.design_6helix
        sb = design.strand(0, 10)
        sb.to(0)
        with self.assertRaises(sc.IllegalDesignError):
            sb.to(5)


class TestCreateHelix(unittest.TestCase):

    def test_helix_constructor_no_max_offset_with_major_ticks(self) -> None:
        # tests bug where an exception is raised if major ticks is defined but not max_offset
        sc.Helix(major_ticks=[0, 5, 10])


class TestM13(unittest.TestCase):

    def test_p7249(self) -> None:
        p7249 = sc.m13()
        self.assertEqual('TTCCCTTCCTTTCTCG', p7249[:16])
        self.assertEqual(7249, len(p7249))
        p7249 = sc.m13(rotation=0)
        self.assertEqual('AATGCTACTACTATTA', p7249[:16])
        self.assertEqual(7249, len(p7249))

    def test_p7560(self) -> None:
        p7560 = sc.m13(rotation=0, variant=sc.M13Variant.p7560)
        self.assertEqual('AGCTTGGCACTGGCCG', p7560[:16])
        self.assertEqual(7560, len(p7560))

    def test_p8064(self) -> None:
        p8064 = sc.m13(rotation=0, variant=sc.M13Variant.p8064)
        self.assertEqual('GGCAATGACCTGATAG', p8064[:16])
        self.assertEqual(8064, len(p8064))


class TestModifications(unittest.TestCase):

    def test_to_json__names_unique_for_modifications_raises_no_error(self) -> None:
        helices = [sc.Helix(max_offset=100)]
        design: sc.Design = sc.Design(helices=helices, strands=[], grid=sc.square)
        name = 'mod_name'
        design.strand(0, 0).move(5).with_modification_5p(sc.Modification5Prime(display_text=name, id=name))
        design.strand(0, 5).move(5).with_modification_3p(
            sc.Modification3Prime(display_text=name, id=name + '3'))
        design.to_json(True)

    def test_to_json__names_not_unique_for_modifications_raises_error(self) -> None:
        helices = [sc.Helix(max_offset=100)]
        design: sc.Design = sc.Design(helices=helices, strands=[], grid=sc.square)
        name = 'mod_name'
        design.strand(0, 0).move(5).with_modification_5p(sc.Modification5Prime(display_text=name, id=name))
        design.strand(0, 5).move(5).with_modification_3p(sc.Modification3Prime(display_text=name, id=name))
        with self.assertRaises(sc.IllegalDesignError):
            design.to_json(True)

    def test_mod_illegal_exceptions_raised(self) -> None:
        strand = sc.Strand(domains=[sc.Domain(0, True, 0, 5)], dna_sequence='AATGC')
        strand.set_modification_internal(2, mod.biotin_int)
        with self.assertRaises(sc.IllegalDesignError):
            strand.set_modification_internal(1, mod.biotin_int)

        # biotin3_1 = mod.Biotin(location=sc.ModLocation.prime3)
        # biotin3_2 = mod.Biotin(location=sc.ModLocation.prime3)
        # with self.assertRaises(sc.IllegalDesignError):
        #     strand = sc.Strand(domains=[sc.Substrand(0, True, 0, 5)], dna_sequence='AATGC',
        #                        modifications=[biotin3_1, biotin3_2])
        #
        # biotin_i_1 = mod.Biotin(location=sc.ModLocation.internal, offset=2)
        # biotin_i_2 = mod.Biotin(location=sc.ModLocation.internal, offset=2)
        # with self.assertRaises(sc.IllegalDesignError):
        #     strand = sc.Strand(domains=[sc.Substrand(0, True, 0, 5)], dna_sequence='AATGC',
        #                        modifications=[biotin_i_1, biotin_i_2])
        #
        # biotin_i_small = mod.Biotin(location=sc.ModLocation.internal, offset=-1)
        # with self.assertRaises(sc.IllegalDesignError):
        #     strand = sc.Strand(domains=[sc.Substrand(0, True, 0, 5)], dna_sequence='AATGC',
        #                        modifications=[biotin_i_small])
        #     seq = strand.idt_dna_sequence()
        #
        # biotin_i_large = mod.Biotin(location=sc.ModLocation.internal, offset=10)
        # with self.assertRaises(sc.IllegalDesignError):
        #     strand = sc.Strand(domains=[sc.Substrand(0, True, 0, 5)], dna_sequence='AATGC',
        #                        modifications=[biotin_i_small])
        #     seq = strand.idt_dna_sequence()
        #
        # biotin_i_offset_not_T = mod.Biotin(location=sc.ModLocation.internal, offset=0)
        # with self.assertRaises(sc.IllegalDesignError):
        #     strand = sc.Strand(domains=[sc.Substrand(0, True, 0, 5)], dna_sequence='AATGC',
        #                        modifications=[biotin_i_offset_not_T])
        #     seq = strand.idt_dna_sequence()
        #
        # cy3I_offset_off_end = mod.Cy3(location=sc.ModLocation.internal, offset=4)
        # with self.assertRaises(sc.IllegalDesignError):
        #     strand = sc.Strand(domains=[sc.Substrand(0, True, 0, 5)], dna_sequence='AATGC',
        #                        modifications=[cy3I_offset_off_end])
        #     seq = strand.idt_dna_sequence()

    def test_Cy3(self) -> None:
        cy3_5 = mod.cy3_5p
        self.assertEqual(r'/5Cy3/', cy3_5.idt_text)
        self.assertEqual(r'/5Cy3/', cy3_5.id)
        self.assertEqual('Cy3', cy3_5.display_text)
        cy3_3 = mod.cy3_3p
        self.assertEqual(r'/3Cy3Sp/', cy3_3.idt_text)
        self.assertEqual(r'/3Cy3Sp/', cy3_3.id)
        self.assertEqual('Cy3', cy3_3.display_text)
        # cy3_i1 = mod.Cy3(location=sc.ModLocation.internal, offset=1)
        cy3_i1 = mod.cy3_int
        self.assertEqual(r'/iCy3/', cy3_i1.idt_text)
        self.assertEqual(r'/iCy3/', cy3_i1.id)
        self.assertEqual('Cy3', cy3_i1.display_text)
        # cy3_i2 = mod.Cy3(location=sc.ModLocation.internal, offset=3)
        cy3_i2 = mod.cy3_int
        self.assertEqual(r'/iCy3/', cy3_i2.idt_text)
        self.assertEqual(r'/iCy3/', cy3_i2.id)
        self.assertEqual('Cy3', cy3_i2.display_text)

        strand5 = sc.Strand(domains=[sc.Domain(0, True, 0, 5)], dna_sequence='ATTGC',
                            modification_5p=cy3_5)
        strand3 = sc.Strand(domains=[sc.Domain(1, True, 0, 5)], dna_sequence='ATTGC',
                            modification_3p=cy3_3)
        strand_i = sc.Strand(domains=[sc.Domain(2, True, 0, 5)], dna_sequence='ATTGC',
                             modifications_int={1: cy3_i1, 3: cy3_i2})
        strand53 = sc.Strand(domains=[sc.Domain(3, True, 0, 5)], dna_sequence='ATTGC',
                             modification_5p=cy3_5, modification_3p=cy3_3)
        strand53_i = sc.Strand(domains=[sc.Domain(4, True, 0, 5)], dna_sequence='ATTGC',
                               modification_5p=cy3_5, modification_3p=cy3_3,
                               modifications_int={1: cy3_i1, 3: cy3_i2})

        self.assertEqual(r'/5Cy3/ATTGC', strand5.idt_dna_sequence())
        self.assertEqual(r'ATTGC/3Cy3Sp/', strand3.idt_dna_sequence())
        self.assertEqual(r'/5Cy3/ATTGC/3Cy3Sp/', strand53.idt_dna_sequence())
        self.assertEqual(r'AT/iCy3/TG/iCy3/C', strand_i.idt_dna_sequence())
        self.assertEqual(r'/5Cy3/AT/iCy3/TG/iCy3/C/3Cy3Sp/', strand53_i.idt_dna_sequence())

    def test_biotin(self) -> None:
        biotin5 = mod.biotin_5p
        self.assertEqual(r'/5Biosg/', biotin5.idt_text)
        self.assertEqual(r'/5Biosg/', biotin5.id)
        self.assertEqual('B', biotin5.display_text)
        biotin3 = mod.biotin_3p
        self.assertEqual(r'/3Bio/', biotin3.idt_text)
        self.assertEqual(r'/3Bio/', biotin3.id)
        self.assertEqual('B', biotin3.display_text)
        # biotin_i_1 = mod.Biotin(location=sc.ModLocation.internal, offset=1)
        biotin_i_1 = mod.biotin_int
        self.assertEqual(r'/iBiodT/', biotin_i_1.idt_text)
        self.assertEqual(r'/iBiodT/', biotin_i_1.id)
        self.assertEqual('B', biotin_i_1.display_text)
        # biotin_i_2 = mod.Biotin(location=sc.ModLocation.internal, offset=2)
        biotin_i_2 = mod.biotin_int
        self.assertEqual(r'/iBiodT/', biotin_i_2.idt_text)
        self.assertEqual(r'/iBiodT/', biotin_i_2.id)
        self.assertEqual('B', biotin_i_2.display_text)

        strand5 = sc.Strand(domains=[sc.Domain(0, True, 0, 5)], dna_sequence='ATTGC',
                            modification_5p=biotin5)
        strand3 = sc.Strand(domains=[sc.Domain(1, True, 0, 5)], dna_sequence='ATTGC',
                            modification_3p=biotin3)
        strand_i = sc.Strand(domains=[sc.Domain(2, True, 0, 5)], dna_sequence='ATTGC',
                             modifications_int={1: biotin_i_1, 2: biotin_i_2})
        strand53 = sc.Strand(domains=[sc.Domain(3, True, 0, 5)], dna_sequence='ATTGC',
                             modification_5p=biotin5, modification_3p=biotin3)
        strand53i = sc.Strand(domains=[sc.Domain(4, True, 0, 5)], dna_sequence='ATTGC',
                              modification_5p=biotin5, modification_3p=biotin3,
                              modifications_int={1: biotin_i_1, 2: biotin_i_2})
        self.assertEqual(r'/5Biosg/ATTGC', strand5.idt_dna_sequence())
        self.assertEqual(r'ATTGC/3Bio/', strand3.idt_dna_sequence())
        self.assertEqual(r'A/iBiodT//iBiodT/GC', strand_i.idt_dna_sequence())
        self.assertEqual(r'/5Biosg/ATTGC/3Bio/', strand53.idt_dna_sequence())
        self.assertEqual(r'/5Biosg/A/iBiodT//iBiodT/GC/3Bio/', strand53i.idt_dna_sequence())

    def test_to_json_serializable(self) -> None:
        biotin5 = mod.biotin_5p
        self.assertEqual(r'/5Biosg/', biotin5.idt_text)
        self.assertEqual(r'/5Biosg/', biotin5.id)
        self.assertEqual('B', biotin5.display_text)
        biotin3 = mod.biotin_3p
        self.assertEqual(r'/3Bio/', biotin3.idt_text)
        self.assertEqual(r'/3Bio/', biotin3.id)
        self.assertEqual('B', biotin3.display_text)
        # biotin_i_1 = mod.Biotin(location=sc.ModLocation.internal, offset=1)
        biotin_i_1 = mod.biotin_int
        self.assertEqual(r'/iBiodT/', biotin_i_1.idt_text)
        self.assertEqual(r'/iBiodT/', biotin_i_1.id)
        self.assertEqual('B', biotin_i_1.display_text)
        # biotin_i_2 = mod.Biotin(location=sc.ModLocation.internal, offset=2)
        biotin_i_2 = mod.biotin_int
        self.assertEqual(r'/iBiodT/', biotin_i_2.idt_text)
        self.assertEqual(r'/iBiodT/', biotin_i_2.id)
        self.assertEqual('B', biotin_i_2.display_text)

        strand5 = sc.Strand(domains=[sc.Domain(0, True, 0, 5)], dna_sequence='ATTGC',
                            modification_5p=biotin5)
        strand3 = sc.Strand(domains=[sc.Domain(1, True, 0, 5)], dna_sequence='ATTGC',
                            modification_3p=biotin3)
        strand_i = sc.Strand(domains=[sc.Domain(2, True, 0, 5)], dna_sequence='ATTGC',
                             modifications_int={1: biotin_i_1, 2: biotin_i_2})
        strand53 = sc.Strand(domains=[sc.Domain(3, True, 0, 5)], dna_sequence='ATTGC',
                             modification_5p=biotin5, modification_3p=biotin3)
        strand53_i = sc.Strand(domains=[sc.Domain(4, True, 0, 5)], dna_sequence='ATTGC',
                               modification_5p=biotin5, modification_3p=biotin3,
                               modifications_int={1: biotin_i_1, 2: biotin_i_2})

        strands = [strand5, strand3, strand_i, strand53, strand53_i]
        design = sc.Design(strands=strands, grid=sc.square)

        # print(design.to_json())

        json_dict = design.to_json_serializable(suppress_indent=False)
        self.assertTrue(sc.design_modifications_key in json_dict)
        mods_dict = json_dict[sc.design_modifications_key]
        self.assertTrue(r'/5Biosg/' in mods_dict)
        self.assertTrue(r'/3Bio/' in mods_dict)
        self.assertTrue(r'/iBiodT/' in mods_dict)

        strand5_mod5_json = json_dict[sc.strands_key][0][sc.modification_5p_key]
        strand3_mod3_json = json_dict[sc.strands_key][1][sc.modification_3p_key]
        self.assertEqual("/5Biosg/", strand5_mod5_json)
        self.assertEqual("/3Bio/", strand3_mod3_json)

        strand_i_mods_int_json = json_dict[sc.strands_key][2][sc.modifications_int_key]
        self.assertDictEqual({"1": "/iBiodT/", "2": "/iBiodT/"}, strand_i_mods_int_json)

        strand53_mod5_json = json_dict[sc.strands_key][3][sc.modification_5p_key]
        strand53_mod3_json = json_dict[sc.strands_key][3][sc.modification_3p_key]
        self.assertEqual("/5Biosg/", strand53_mod5_json)
        self.assertEqual("/3Bio/", strand53_mod3_json)

        strand53_i_mod5_json = json_dict[sc.strands_key][4][sc.modification_5p_key]
        strand53_i_mod3_json = json_dict[sc.strands_key][4][sc.modification_3p_key]
        strand53_i_mods_int_json = json_dict[sc.strands_key][4][sc.modifications_int_key]
        self.assertEqual("/5Biosg/", strand53_i_mod5_json)
        self.assertEqual("/3Bio/", strand53_i_mod3_json)
        self.assertDictEqual({"1": "/iBiodT/", "2": "/iBiodT/"}, strand53_i_mods_int_json)


class TestImportCadnanoV2(unittest.TestCase):
    """
    Tests the import feature to cadnano v2 (see misc/cadnano-format-specs/v2.txt).
    """
    folder = "cadnano_v2_import"
    input_path = os.path.join('tests_inputs', folder)
    output_path = os.path.join('tests_outputs', folder)

    def test_32_helix_rectangle(self) -> None:
        design = sc.Design.from_cadnano_v2(directory=self.input_path,
                                           filename='test_32_helix_rectangle.json')
        design.write_scadnano_file(directory=self.output_path,
                                   filename=f'test_32_helix_rectangle.{sc.default_scadnano_file_extension}')

    def test_helices_order(self) -> None:
        design = sc.Design.from_cadnano_v2(directory=self.input_path,
                                           filename='test_helices_order.json')
        design.write_scadnano_file(directory=self.output_path,
                                   filename=f'test_helices_order.{sc.default_scadnano_file_extension}')

    def test_huge_hex(self) -> None:
        design = sc.Design.from_cadnano_v2(directory=self.input_path,
                                           filename='test_huge_hex.json')
        design.write_scadnano_file(directory=self.output_path,
                                   filename=f'test_huge_hex.{sc.default_scadnano_file_extension}')

    def test_Science09_prot120_98_v3(self) -> None:
        file_name = "test_Science09_prot120_98_v3"
        design = sc.Design.from_cadnano_v2(directory=self.input_path,
                                           filename=file_name + ".json")
        design.write_scadnano_file(directory=self.output_path,
                                   filename=f'{file_name}.{sc.default_scadnano_file_extension}')

    def test_Nature09_monolith(self) -> None:
        file_name = "test_Nature09_monolith"
        design = sc.Design.from_cadnano_v2(directory=self.input_path,
                                           filename=file_name + ".json")
        design.write_scadnano_file(directory=self.output_path,
                                   filename=f'{file_name}.{sc.default_scadnano_file_extension}')

    def test_circular_auto_staple(self) -> None:
        file_name = "test_circular_auto_staple"
        design = sc.Design.from_cadnano_v2(directory=self.input_path,
                                           filename=file_name + ".json")
        design.write_scadnano_file(directory=self.output_path,
                                   filename=f'{file_name}.{sc.default_scadnano_file_extension}')

    def test_circular_auto_staple_hex(self) -> None:
        file_name = "test_circular_auto_staple_hex"
        design = sc.Design.from_cadnano_v2(directory=self.input_path,
                                           filename=file_name + ".json")
        design.write_scadnano_file(directory=self.output_path,
                                   filename=f'{file_name}.{sc.default_scadnano_file_extension}')

    def test_paranemic_crossover(self) -> None:
        file_name = "test_paranemic_crossover"
        design = sc.Design.from_cadnano_v2(directory=self.input_path,
                                           filename=file_name + ".json")
        design.write_scadnano_file(directory=self.output_path,
                                   filename=f'{file_name}.{sc.default_scadnano_file_extension}')


class TestExportDNASequences(unittest.TestCase):

    def setUp(self) -> None:
        r""" Removing scaffold from this design:
    0        8        16       24       32       40       48       56       64       72       80       88      96
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

gives this (removed scaffold to make it easier for me to visually track the strands)

    0        8        16       24       32       40       48       56       64       72       80       88      96
0   <-------A-------+ +-------D--------D-------] <-------G--------G--------G-------] <------+ +-------P-------]
                    | |                                                                     | |
1   [-------A-------+ +------> [------+ +-------H--------H-------> [------+ +-------M-------+ +-------P------->
                                      | |                                 | |
2   <-------B-------+ +-------E-------+ +------] <-------I--------I-------+ +------] <------+ +-------Q-------]
                    | |                                                                     | |
3   [-------B-------+ +------> [------+ +-------J--------J-------> [------+ +-------N-------+ +-------Q------->
                                      | |                                 | |
4   <-------C-------+ +-------F-------+ +------] <-------K--------K-------+ +------] <------+ +-------R-------]
                    | |                                                                     | |
5   [-------C-------+ +------> [-------L--------L--------L-------> [-------O--------O-------+ +-------R------->

We give the strands simple single-letter names to help test they appear in the correct order in the output.

row major 5':                    DGPAEIHMQBFKJNRCLO
col major 5':                    ABCEFLDHJIKOGMNPQR
row major 3':                    AGMDHPBINEJQCKOFLR
col major 3':                    ABCDEFGIKHJLMNOPQR
row major 5' or 3':              ADGMPEHIBNQFJKCORL
col major 5' or 3':              ABCDEFLHJGIKOMNPQR
row major top-left domain start: ADGMPEHIBNQFJKCORL
col major top-left domain start: ABCDEFLHJGIKMNOPQR
        """
        helices = [sc.Helix(max_offset=100) for _ in range(6)]
        self.design_6h: sc.Design = sc.Design(helices=helices, strands=[], grid=sc.square)
        d = self.design_6h

        d.strand(1, 0).move(16).cross(0).move(-16).with_name('A')
        d.strand(3, 0).move(16).cross(2).move(-16).with_name('B')
        d.strand(5, 0).move(16).cross(4).move(-16).with_name('C')

        d.strand(0, 40).move(-24).cross(1).move(8).with_name('D')

        d.strand(1, 24).move(8).cross(2).move(-16).cross(3).move(8).with_name('E')
        d.strand(3, 24).move(8).cross(4).move(-16).cross(5).move(8).with_name('F')

        d.strand(0, 72).move(-32).with_name('G')

        d.strand(2, 40).move(-8).cross(1).move(24).with_name('H')
        d.strand(1, 56).move(8).cross(2).move(-24).with_name('I')

        d.strand(4, 40).move(-8).cross(3).move(24).with_name('J')
        d.strand(3, 56).move(8).cross(4).move(-24).with_name('K')

        d.strand(5, 24).move(32).with_name('L')

        d.strand(2, 72).move(-8).cross(1).move(16).cross(0).move(-8).with_name('M')
        d.strand(4, 72).move(-8).cross(3).move(16).cross(2).move(-8).with_name('N')

        d.strand(5, 56).move(24).cross(4).move(-8).with_name('O')

        d.strand(0, 96).move(-16).cross(1).move(16).with_name('P')
        d.strand(2, 96).move(-16).cross(3).move(16).with_name('Q')
        d.strand(4, 96).move(-16).cross(5).move(16).with_name('R')

        for strand in d.strands:
            d.assign_dna(strand, 'A' * 32, assign_complement=False)

    def _get_names_idt(self, design: sc.Design, key: sc.KeyFunction[sc.Strand]) -> str:
        # call design.to_idt_bulk_input_format with given key functions,
        # get IDT names of strands exported, and return them joined into a single string
        idt_str = design.to_idt_bulk_input_format(key=key)
        idt_lines = idt_str.split('\n')
        names = []
        for line in idt_lines:
            name = line.split(',')[0]
            names.append(name)
        names_joined = ''.join(names)
        return names_joined

    def test_to_idt_bulk_input_format__row_major_5p(self) -> None:
        key = sc.strand_order_key_function(column_major=False, strand_order=sc.StrandOrder.five_prime)
        names_joined = self._get_names_idt(self.design_6h, key)
        self.assertEqual('DGPAEIHMQBFKJNRCLO', names_joined)

    def test_to_idt_bulk_input_format__col_major_5p(self) -> None:
        key = sc.strand_order_key_function(column_major=True, strand_order=sc.StrandOrder.five_prime)
        names_joined = self._get_names_idt(self.design_6h, key)
        self.assertEqual('ABCEFLDHJIKOGMNPQR', names_joined)

    def test_to_idt_bulk_input_format__row_major_3p(self) -> None:
        key = sc.strand_order_key_function(column_major=False, strand_order=sc.StrandOrder.three_prime)
        names_joined = self._get_names_idt(self.design_6h, key)
        self.assertEqual('AGMDHPBINEJQCKOFLR', names_joined)

    def test_to_idt_bulk_input_format__col_major_3p(self) -> None:
        key = sc.strand_order_key_function(column_major=True, strand_order=sc.StrandOrder.three_prime)
        names_joined = self._get_names_idt(self.design_6h, key)
        self.assertEqual('ABCDEFGIKHJLMNOPQR', names_joined)

    def test_to_idt_bulk_input_format__row_major_5p_or_3p(self) -> None:
        key = sc.strand_order_key_function(column_major=False,
                                           strand_order=sc.StrandOrder.five_or_three_prime)
        names_joined = self._get_names_idt(self.design_6h, key)
        self.assertEqual('ADGMPEHIBNQFJKCORL', names_joined)

    def test_to_idt_bulk_input_format__col_major_5p_or_3p(self) -> None:
        key = sc.strand_order_key_function(column_major=True, strand_order=sc.StrandOrder.five_or_three_prime)
        names_joined = self._get_names_idt(self.design_6h, key)
        self.assertEqual('ABCDEFLHJGIKOMNPQR', names_joined)

    def test_to_idt_bulk_input_format__row_major_top_left_domain_start(self) -> None:
        key = sc.strand_order_key_function(column_major=False, strand_order=sc.StrandOrder.top_left_domain)
        names_joined = self._get_names_idt(self.design_6h, key)
        self.assertEqual('ADGMPEHIBNQFJKCORL', names_joined)

    def test_to_idt_bulk_input_format__col_major_top_left_domain_start(self) -> None:
        key = sc.strand_order_key_function(column_major=True, strand_order=sc.StrandOrder.top_left_domain)
        names_joined = self._get_names_idt(self.design_6h, key)
        self.assertEqual('ABCDEFLHJGIKMNOPQR', names_joined)

    def test_to_idt_bulk_input_format__duplicate_names_same_sequence(self) -> None:
        length = 8
        helices = [sc.Helix(max_offset=length)]
        ss1_r = sc.Domain(0, True, 0, 4)
        ss2_r = sc.Domain(0, True, 4, 8)
        ss_l = sc.Domain(0, False, 0, 4)

        s1_r = sc.Strand([ss1_r], idt=sc.IDTFields(), name='s1_r')
        s2_r = sc.Strand([ss2_r], idt=sc.IDTFields(), name='s1_r')
        s_l = sc.Strand([ss_l], idt=sc.IDTFields(), name='s_l')

        strands = [s1_r, s2_r, s_l]

        design = sc.Design(helices=helices, strands=strands, grid=sc.square)

        design.assign_dna(s_l, 'AGTT')
        design.assign_dna(s2_r, 'AACT')

        # should not raise exception
        design.to_idt_bulk_input_format()

    def test_to_idt_bulk_input_format__duplicate_names_different_sequences(self) -> None:
        ss1_r = sc.Domain(0, True, 0, 4)
        ss2_r = sc.Domain(0, True, 4, 8)
        ss_l = sc.Domain(0, False, 0, 4)

        s1_r = sc.Strand([ss1_r], idt=sc.IDTFields(), name='s1_r')
        s2_r = sc.Strand([ss2_r], idt=sc.IDTFields(), name='s1_r')
        s_l = sc.Strand([ss_l], idt=sc.IDTFields(), name='s_l')

        strands = [s1_r, s2_r, s_l]

        design = sc.Design(strands=strands, grid=sc.square)

        design.assign_dna(s_l, 'AGTT')
        design.assign_dna(s2_r, 'GGGG')

        with self.assertRaises(sc.IllegalDesignError):
            design.to_idt_bulk_input_format()

    def test_to_idt_bulk_input_format__duplicate_names_different_scales(self) -> None:
        ss1_r = sc.Domain(0, True, 0, 4)
        ss2_r = sc.Domain(0, True, 4, 8)
        ss_l = sc.Domain(0, False, 0, 4)

        s1_r = sc.Strand([ss1_r], idt=sc.IDTFields(scale='25nm'), name='s1_r')
        s2_r = sc.Strand([ss2_r], idt=sc.IDTFields(scale='100nm'), name='s1_r')
        s_l = sc.Strand([ss_l], idt=sc.IDTFields(), name='s_l')

        strands = [s1_r, s2_r, s_l]

        design = sc.Design(strands=strands, grid=sc.square)

        design.assign_dna(s_l, 'AGTT')
        design.assign_dna(s2_r, 'AACT')

        with self.assertRaises(sc.IllegalDesignError):
            design.to_idt_bulk_input_format()

    # https://stackoverflow.com/a/46307456/5339430
    # @mock.patch('sys.stdout', new_callable=io.StringIO)
    def test_to_idt_bulk_input_format__duplicate_names_different_purifications(self) -> None:
        ss1_r = sc.Domain(0, True, 0, 4)
        ss2_r = sc.Domain(0, True, 4, 8)
        ss_l = sc.Domain(0, False, 0, 4)

        s1_r = sc.Strand([ss1_r], idt=sc.IDTFields(purification='STD'), name='s1_r')
        s2_r = sc.Strand([ss2_r], idt=sc.IDTFields(purification='HPLC'), name='s1_r')
        s_l = sc.Strand([ss_l], idt=sc.IDTFields(), name='s_l')

        strands = [s1_r, s2_r, s_l]

        design = sc.Design(strands=strands, grid=sc.square)

        design.assign_dna(s_l, 'AGTT')
        design.assign_dna(s2_r, 'AACT')

        with self.assertRaises(sc.IllegalDesignError):
            design.to_idt_bulk_input_format(warn_duplicate_name=True)
        # printed = mock_print.getvalue()
        # self.assertIn('two strands with same IDT name', printed)
        # self.assertIn('s1_r', printed)

    def test_write_idt_plate_excel_file(self) -> None:
        strand_len = 10

        # add 10 strands in excess of 3 plates
        for plate_type in [sc.PlateType.wells96, sc.PlateType.wells384]:
            filename = f'test_excel_export_{plate_type.num_wells_per_plate()}.xls'
            helices = [sc.Helix(max_offset=100) for _ in range(1)]
            design = sc.Design(helices=helices, strands=[], grid=sc.square)
            for strand_idx in range(3 * plate_type.num_wells_per_plate() + 10):
                design.strand(0, strand_len * strand_idx).move(strand_len).with_name(f's{strand_idx}')
                design.strands[-1].set_dna_sequence('T' * strand_len)

            design.write_idt_plate_excel_file(filename=filename, plate_type=plate_type)

            book = xlrd.open_workbook(filename)
            self.assertEqual(4, book.nsheets)
            for plate in range(4):
                sheet = book.sheet_by_index(plate)
                self.assertEqual(3, sheet.ncols)

                if plate == 2:  # penultimate plate
                    expected_wells = plate_type.num_wells_per_plate() - plate_type.min_wells_per_plate() + 10
                elif plate == 3:  # last plate
                    expected_wells = plate_type.min_wells_per_plate()
                else:
                    expected_wells = plate_type.num_wells_per_plate()

                self.assertEqual(expected_wells + 1, sheet.nrows)

            os.remove(filename)


class TestExportCadnanoV2(unittest.TestCase):
    """
    Tests the export feature to cadnano v2 (see misc/cadnano-format-specs/v2.txt).
    """
    folder = "cadnano_v2_export"
    input_path = os.path.join('tests_inputs', folder)
    output_path = os.path.join('tests_outputs', folder)
    ext = sc.default_scadnano_file_extension

    def test_export_design_with_helix_group(self):
        e = 'east'
        s = 'south'
        helices = [
            sc.Helix(max_offset=24, group=s),
            sc.Helix(max_offset=25, group=s),
        ]
        helices.extend([
            sc.Helix(max_offset=22, group=e),
            sc.Helix(max_offset=23, group=e),
        ])

        group_south = sc.HelixGroup(position=sc.Position3D(x=0, y=10, z=0),
                                    grid=sc.square)
        group_east = sc.HelixGroup(position=sc.Position3D(x=10, y=0, z=0), grid=sc.square)

        groups = {
            e: group_east,
            s: group_south
        }

        design = sc.Design(helices=helices, groups=groups, strands=[])
        design.write_scadnano_file(directory=self.input_path,
                                   filename=f'test_export_design_with_helix_group.{self.ext}')
        design.export_cadnano_v2(directory=self.output_path,
                                 filename='test_export_design_with_helix_group.json')

    def test_export_design_with_helix_group_not_same_grid(self):
        e = 'east'
        s = 'south'
        helices = [
            sc.Helix(max_offset=24, group=s),
            sc.Helix(max_offset=25, group=s),
        ]
        helices.extend([
            sc.Helix(max_offset=22, group=e),
            sc.Helix(max_offset=23, group=e),
        ])

        group_south = sc.HelixGroup(position=sc.Position3D(x=0, y=10, z=0),
                                    grid=sc.square)
        group_east = sc.HelixGroup(position=sc.Position3D(x=10, y=0, z=0), grid=sc.honeycomb)

        groups = {
            e: group_east,
            s: group_south
        }

        design = sc.Design(helices=helices, groups=groups, strands=[])
        design.write_scadnano_file(directory=self.input_path,
                                   filename=f'test_export_design_with_helix_group_not_same_grid.{self.ext}')

        with self.assertRaises(ValueError) as context:
            design.export_cadnano_v2(directory=self.output_path,
                                     filename='test_export_design_with_helix_group_not_same_grid.json')
        self.assertTrue('helix groups' in context.exception.args[0])

    def test_2_staple_2_helix_origami_extremely_simple(self):
        helices = [sc.Helix(max_offset=32), sc.Helix(max_offset=32)]
        scaf_part = sc.Domain(helix=0, forward=True, start=0, end=32)
        scaf = sc.Strand(domains=[scaf_part], is_scaffold=True)
        design = sc.Design(helices=helices, strands=[scaf], grid=sc.square)
        design.write_scadnano_file(directory=self.input_path,
                                   filename=f'test_2_stape_2_helix_origami_extremely_simple.{self.ext}')
        design.export_cadnano_v2(directory=self.output_path,
                                 filename='test_2_stape_2_helix_origami_extremely_simple.json')

    def test_2_staple_2_helix_origami_extremely_simple_2(self) -> None:
        helices = [sc.Helix(max_offset=32), sc.Helix(max_offset=32)]
        scaf_part1 = sc.Domain(helix=0, forward=True, start=0, end=32)
        scaf_part2 = sc.Domain(helix=1, forward=False, start=0, end=32)
        scaf = sc.Strand(domains=[scaf_part1, scaf_part2], is_scaffold=True)
        design = sc.Design(helices=helices, strands=[scaf], grid=sc.square)
        design.write_scadnano_file(directory=self.input_path,
                                   filename=f'test_2_stape_2_helix_origami_extremely_simple_2.{self.ext}')
        design.export_cadnano_v2(directory=self.output_path,
                                 filename='test_2_stape_2_helix_origami_extremely_simple_2.json')

    def test_2_staple_2_helix_origami_deletions_insertions(self) -> None:
        # left staple
        stap_left_ss1 = sc.Domain(helix=1, forward=True, start=0, end=16)
        stap_left_ss0 = sc.Domain(helix=0, forward=False, start=0, end=16)
        stap_left = sc.Strand(domains=[stap_left_ss1, stap_left_ss0])

        # right staple
        stap_right_ss0 = sc.Domain(helix=0, forward=False, start=16, end=32)
        stap_right_ss1 = sc.Domain(helix=1, forward=True, start=16, end=32)
        stap_right = sc.Strand(domains=[stap_right_ss0, stap_right_ss1])

        # scaffold
        scaf_ss1_left = sc.Domain(helix=1, forward=False, start=0, end=16)
        scaf_ss0 = sc.Domain(helix=0, forward=True, start=0, end=32)
        # loopout = sc.Loopout(length=3) No loopout in cadnano
        scaf_ss1_right = sc.Domain(helix=1, forward=False, start=16, end=32)
        scaf = sc.Strand(domains=[scaf_ss1_left, scaf_ss0, scaf_ss1_right], is_scaffold=True)

        # whole design
        design = sc.Design(strands=[scaf, stap_left, stap_right], grid=sc.square)

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
        design.write_scadnano_file(directory=self.input_path,
                                   filename=f'test_2_stape_2_helix_origami_deletions_insertions.{self.ext}')
        design.export_cadnano_v2(directory=self.output_path,
                                 filename='test_2_stape_2_helix_origami_deletions_insertions.json')

    def test_6_helix_origami_rectangle(self) -> None:
        design = rect.create(num_helices=6, num_cols=10, nick_pattern=rect.staggered,
                             twist_correction_deletion_spacing=3)
        design.write_scadnano_file(directory=self.input_path,
                                   filename=f'test_6_helix_origami_rectangle.{self.ext}')
        design.export_cadnano_v2(directory=self.output_path,
                                 filename='test_6_helix_origami_rectangle.json')

    def test_6_helix_bundle_honeycomb(self) -> None:
        design = sc.Design.from_scadnano_file(
            os.path.join(self.input_path, f'test_6_helix_bundle_honeycomb.{self.ext}'))
        design.export_cadnano_v2(directory=self.output_path,
                                 filename='test_6_helix_bundle_honeycomb.json')

    def test_16_helix_origami_rectangle_no_twist(self) -> None:
        design = rect.create(num_helices=16, num_cols=26, assign_seq=True,
                             twist_correction_deletion_spacing=3)
        design.write_scadnano_file(directory=self.input_path,
                                   filename=f'test_16_helix_origami_rectangle_no_twist.{self.ext}')
        design.export_cadnano_v2(directory=self.output_path,
                                 filename='test_16_helix_origami_rectangle_no_twist.json')

    def test_circular_strand(self) -> None:
        helices = [sc.Helix(max_offset=24) for _ in range(2)]
        design = sc.Design(helices=helices, grid=sc.square)

        design.strand(1, 0).move(8).cross(0).move(-8).as_circular()
        design.write_scadnano_file(directory=self.input_path,
                                   filename=f'test_circular_strand.{self.ext}')
        design.export_cadnano_v2(directory=self.output_path,
                                 filename='test_circular_strand.json')

    def test_big_circular_staples_hex(self) -> None:
        design = sc.Design.from_scadnano_file(
            os.path.join(self.input_path, f'test_big_circular_staples_hex.{self.ext}'))
        design.export_cadnano_v2(directory=self.output_path,
                                 filename='test_big_circular_staples_hex.json')

    def test_big_circular_staples(self) -> None:
        design = sc.Design.from_scadnano_file(
            os.path.join(self.input_path, f'test_big_circular_staples.{self.ext}'))
        design.export_cadnano_v2(directory=self.output_path,
                                 filename='test_big_circular_staples.json')

    def test_paranemic_crossover(self) -> None:
        design = sc.Design.from_scadnano_file(
            os.path.join(self.input_path, f'test_paranemic_crossover.{self.ext}'))
        design.export_cadnano_v2(directory=self.output_path,
                                 filename='test_paranemic_crossover.json')

    def test_bad_cases(self) -> None:
        """ We do not handle Loopouts and design where the parity of the helix
        does not correspond to the direction.
        """

        # Bad case one: parity issue in design (see cadnano v2 format spec, v2.txt)
        helices = [sc.Helix(max_offset=32), sc.Helix(max_offset=32)]
        scaf_part = sc.Domain(helix=1, forward=True, start=0, end=32)
        scaf = sc.Strand(domains=[scaf_part], is_scaffold=True)
        design = sc.Design(helices=helices, strands=[scaf], grid=sc.square)

        with self.assertRaises(ValueError) as context:
            design.export_cadnano_v2(directory=self.output_path,
                                     filename='test_parity_issue.json')
        self.assertTrue('forward' in context.exception.args[0])

        # Bad case two: Loopouts
        helices = [sc.Helix(max_offset=48), sc.Helix(max_offset=48)]

        # left staple
        stap_left_ss1 = sc.Domain(helix=1, forward=True, start=8, end=24)
        stap_left_ss0 = sc.Domain(helix=0, forward=False, start=8, end=24)
        stap_left = sc.Strand(domains=[stap_left_ss1, stap_left_ss0])

        # right staple
        stap_right_ss0 = sc.Domain(helix=0, forward=False, start=24, end=40)
        stap_right_ss1 = sc.Domain(helix=1, forward=True, start=24, end=40)
        stap_right = sc.Strand(domains=[stap_right_ss0, stap_right_ss1])

        # scaffold
        scaf_ss1_left = sc.Domain(helix=1, forward=False, start=8, end=24)
        scaf_ss0 = sc.Domain(helix=0, forward=True, start=8, end=40)
        loopout = sc.Loopout(length=3)
        scaf_ss1_right = sc.Domain(helix=1, forward=False, start=24, end=40)
        scaf = sc.Strand(domains=[scaf_ss1_left, scaf_ss0, loopout, scaf_ss1_right], is_scaffold=True)

        # whole design
        design = sc.Design(helices=helices, strands=[scaf, stap_left, stap_right], grid=sc.square)

        # deletions and insertions added to design are added to both strands on a helix
        design.add_deletion(helix=1, offset=20)
        design.add_insertion(helix=0, offset=14, length=1)
        design.add_insertion(helix=0, offset=26, length=2)

        with self.assertRaises(ValueError) as context:
            design.export_cadnano_v2(directory=self.output_path,
                                     filename='test_loopout_issue.json')
        self.assertTrue('Loopouts' in context.exception.args[0])


class TestDesignFromJson(unittest.TestCase):
    """
    Tests reading a design from a dict derived from JSON.
    """

    def setUp(self) -> None:
        r"""
                0       8       16
                |       |       |
            0   +--X-----------+
               /<--X---++------]\
               |       ||       loopout(3)
               \[---2--++------>/
            1   +---2--]<------+
                """
        st_l = sc.Strand([
            sc.Domain(1, True, 0, 8, insertions=[(4, 2)]),
            sc.Domain(0, False, 0, 8, deletions=[3]),
        ], modification_5p=mod.biotin_5p)
        st_r = sc.Strand([
            sc.Domain(0, False, 8, 16),
            sc.Domain(1, True, 8, 16),
        ], modification_5p=mod.biotin_5p, modification_3p=mod.cy3_3p, modifications_int={
            1: mod.biotin_int, 2: mod.cy5_int
        })
        scaf = sc.Strand([
            sc.Domain(1, False, 0, 8, insertions=[(4, 2)]),
            sc.Domain(0, True, 0, 16, deletions=[3]),
            sc.Loopout(3),
            sc.Domain(1, False, 8, 16, deletions=[]),
        ], is_scaffold=True)

        self.design_pre_json = sc.Design(strands=[st_l, st_r, scaf], grid=sc.square)
        self.design_pre_json.assign_dna(scaf, 'A' * 36)

    def test_from_json__from_and_to_file_contents(self) -> None:
        json_str = self.design_pre_json.to_json()
        json_map = json.loads(json_str)
        design = sc.Design.from_scadnano_json_map(json_map)
        design.to_json_serializable()

    def test_from_json__three_strands(self) -> None:
        json_str = self.design_pre_json.to_json()
        json_map = json.loads(json_str)
        design = sc.Design.from_scadnano_json_map(json_map)

        self.assertTrue(isinstance(design, sc.Design))

        self.assertEqual(sc.Grid.square, design.grid)

        self.assertEqual(2, len(design.helices))
        helix0 = design.helices[0]
        helix1 = design.helices[1]
        self.assertEqual(0, helix0.idx)
        self.assertEqual(0, helix0.min_offset)
        self.assertEqual(16, helix0.max_offset)
        self.assertEqual((0, 0), helix0.grid_position)
        self.assertEqual(1, helix1.idx)
        self.assertEqual(0, helix1.min_offset)
        self.assertEqual(16, helix1.max_offset)
        self.assertEqual((0, 1), helix1.grid_position)

        self.assertEqual(3, len(design.strands))
        st_l = design.strands[0]
        st_r = design.strands[1]
        scaf = design.strands[2]

        self.assertEqual(scaf, design.scaffold)

        self.assertEqual(2, len(st_l.domains))
        self.assertEqual(2, len(st_r.domains))
        self.assertEqual(4, len(scaf.domains))

        self.assertEqual('A' * 36, scaf.dna_sequence)
        self.assertEqual('T' * 17, st_l.dna_sequence)
        self.assertEqual('T' * 16, st_r.dna_sequence)

        st_l_ss0 = st_l.domains[0]
        st_l_ss1 = st_l.domains[1]
        st_r_ss0 = st_r.domains[0]
        st_r_ss1 = st_r.domains[1]
        scaf_ss0 = scaf.domains[0]
        scaf_ss1 = scaf.domains[1]
        scaf_loop = scaf.domains[2]
        scaf_ss2 = scaf.domains[3]

        self.assertEqual(3, scaf_loop.length)

        self.assertEqual(1, st_l_ss0.helix)
        self.assertEqual(0, st_l_ss1.helix)
        self.assertEqual(0, st_r_ss0.helix)
        self.assertEqual(1, st_r_ss1.helix)
        self.assertEqual(1, scaf_ss0.helix)
        self.assertEqual(0, scaf_ss1.helix)
        self.assertEqual(1, scaf_ss2.helix)

        self.assertEqual(True, st_l_ss0.forward)
        self.assertEqual(False, st_l_ss1.forward)
        self.assertEqual(False, st_r_ss0.forward)
        self.assertEqual(True, st_r_ss1.forward)
        self.assertEqual(False, scaf_ss0.forward)
        self.assertEqual(True, scaf_ss1.forward)
        self.assertEqual(False, scaf_ss2.forward)

        self.assertEqual(0, st_l_ss0.start)
        self.assertEqual(8, st_l_ss0.end)
        self.assertEqual(0, st_l_ss1.start)
        self.assertEqual(8, st_l_ss1.end)
        self.assertEqual(8, st_r_ss0.start)
        self.assertEqual(16, st_r_ss0.end)
        self.assertEqual(8, st_r_ss1.start)
        self.assertEqual(16, st_r_ss1.end)
        self.assertEqual(0, scaf_ss0.start)
        self.assertEqual(8, scaf_ss0.end)
        self.assertEqual(0, scaf_ss1.start)
        self.assertEqual(16, scaf_ss1.end)
        self.assertEqual(8, scaf_ss2.start)
        self.assertEqual(16, scaf_ss2.end)

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

        self.assertEqual(mod.biotin_5p, st_l.modification_5p)
        self.assertEqual(None, st_l.modification_3p)
        self.assertDictEqual({}, st_l.modifications_int)

        self.assertEqual(mod.biotin_5p, st_r.modification_5p)
        self.assertEqual(mod.cy3_3p, st_r.modification_3p)
        self.assertDictEqual({1: mod.biotin_int, 2: mod.cy5_int}, st_r.modifications_int)

        self.assertEqual(None, scaf.modification_5p)
        self.assertEqual(None, scaf.modification_3p)
        self.assertDictEqual({}, scaf.modifications_int)

    def test_from_json__helices_non_default_indices(self) -> None:
        h2 = sc.Helix(idx=2)
        h3 = sc.Helix(idx=3)
        h5 = sc.Helix(idx=5)
        helices = [h2, h3, h5]
        s1 = sc.Strand([
            sc.Domain(2, True, 0, 4),
            sc.Domain(3, False, 0, 4),
        ])
        s2 = sc.Strand([
            sc.Domain(3, True, 4, 8),
            sc.Domain(5, False, 4, 8),
        ])
        self.design_pre_json = sc.Design(helices=helices, strands=[s1, s2], grid=sc.square)

        json_str = self.design_pre_json.to_json()
        json_map = json.loads(json_str)
        design = sc.Design.from_scadnano_json_map(json_map)

        self.assertEqual(3, len(design.helices))

        self.assertEqual(2, design.helices[2].idx)
        self.assertEqual(3, design.helices[3].idx)
        self.assertEqual(5, design.helices[5].idx)

    def test_from_json__helices_non_default_indices_mixed_with_default(self) -> None:
        h2 = sc.Helix(idx=2)
        h3 = sc.Helix()
        h5 = sc.Helix(idx=5)
        helices = [h2, h3, h5]
        s1 = sc.Strand([
            sc.Domain(2, True, 0, 4),
            sc.Domain(1, False, 0, 4),
        ])
        s2 = sc.Strand([
            sc.Domain(1, True, 4, 8),
            sc.Domain(5, False, 4, 8),
        ])
        self.design_pre_json = sc.Design(helices=helices, strands=[s1, s2], grid=sc.square)

        json_str = self.design_pre_json.to_json()

        json_map = json.loads(json_str)
        design = sc.Design.from_scadnano_json_map(json_map)

        self.assertEqual(3, len(design.helices))

        self.assertEqual(2, design.helices[2].idx)
        self.assertEqual(1, design.helices[1].idx)
        self.assertEqual(5, design.helices[5].idx)

    def test_from_json__helices_non_default_error_if_some_have_idx_not_others(self) -> None:
        h2 = sc.Helix(idx=2)
        h3 = sc.Helix(idx=3)
        h3_2 = sc.Helix(idx=3)
        helices = [h2, h3, h3_2]
        with self.assertRaises(sc.IllegalDesignError):
            self.design_pre_json = sc.Design(helices=helices, strands=[], grid=sc.square)

    def test_from_json__helices_non_default_error_if_some_have_idx_not_others_mixed_default(self) -> None:
        h0 = sc.Helix()
        h0_2 = sc.Helix(idx=0)
        helices = [h0, h0_2]
        with self.assertRaises(sc.IllegalDesignError):
            self.design_pre_json = sc.Design(helices=helices, strands=[], grid=sc.square)


class TestStrandReversePolarity(unittest.TestCase):
    """
    Tests reversing polarity of all strands in the system.
    """

    def test_reverse_all__one_strand(self) -> None:
        """
        before
        0       8
        |       |
    0   [---X-->

        after
    0   <---X--]
        """
        design = sc.Design(
            strands=[sc.Strand([sc.Domain(0, True, 0, 8, deletions=[4])])],
            grid=sc.square)
        design.reverse_all()

        self.assertEqual(1, len(design.strands))
        strand = design.strands[0]
        self.assertEqual(1, len(strand.domains))
        self.assertEqual(7, strand.offset_5p())
        self.assertEqual(0, strand.offset_3p())
        ss = strand.domains[0]
        self.assertEqual(0, ss.helix)
        self.assertEqual(False, ss.forward)
        self.assertEqual(0, ss.start)
        self.assertEqual(8, ss.end)

    def test_reverse_all__three_strands(self) -> None:
        r"""
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
        design = sc.Design(
            strands=[
                sc.Strand([
                    sc.Domain(1, True, 0, 8, insertions=[(4, 2)]),
                    sc.Domain(0, False, 0, 8, deletions=[3]),
                ]),
                sc.Strand([
                    sc.Domain(0, False, 8, 16),
                    sc.Domain(1, True, 8, 16),
                ]),
                sc.Strand([
                    sc.Domain(1, False, 0, 8, insertions=[(4, 2)]),
                    sc.Domain(0, True, 0, 16, deletions=[3]),
                    sc.Domain(1, False, 8, 16, deletions=[]),
                ]),
            ],
            grid=sc.square)
        design.reverse_all()

        self.assertEqual(3, len(design.strands))
        stap_l = design.strands[0]
        stap_r = design.strands[1]
        scaf = design.strands[2]

        self.assertEqual(0, stap_l.offset_5p())
        self.assertEqual(0, stap_l.offset_3p())
        self.assertEqual(15, stap_r.offset_5p())
        self.assertEqual(15, stap_r.offset_3p())
        self.assertEqual(8, scaf.offset_5p())
        self.assertEqual(7, scaf.offset_3p())

        self.assertEqual(2, len(stap_l.domains))
        stap_l_ss0 = stap_l.domains[0]
        stap_l_ss1 = stap_l.domains[1]

        self.assertEqual(0, stap_l_ss0.helix)
        self.assertEqual(True, stap_l_ss0.forward)
        self.assertEqual(0, stap_l_ss0.start)
        self.assertEqual(8, stap_l_ss0.end)
        self.assertEqual([3], stap_l_ss0.deletions)
        self.assertEqual([], stap_l_ss0.insertions)

        self.assertEqual(1, stap_l_ss1.helix)
        self.assertEqual(False, stap_l_ss1.forward)
        self.assertEqual(0, stap_l_ss1.start)
        self.assertEqual(8, stap_l_ss1.end)
        self.assertEqual([], stap_l_ss1.deletions)
        self.assertEqual([(4, 2)], stap_l_ss1.insertions)

        self.assertEqual(2, len(stap_r.domains))
        stap_r_ss0 = stap_r.domains[0]
        stap_r_ss1 = stap_r.domains[1]

        self.assertEqual(1, stap_r_ss0.helix)
        self.assertEqual(False, stap_r_ss0.forward)
        self.assertEqual(8, stap_r_ss0.start)
        self.assertEqual(16, stap_r_ss0.end)
        self.assertEqual([], stap_r_ss0.deletions)
        self.assertEqual([], stap_r_ss0.insertions)

        self.assertEqual(0, stap_r_ss1.helix)
        self.assertEqual(True, stap_r_ss1.forward)
        self.assertEqual(8, stap_r_ss1.start)
        self.assertEqual(16, stap_r_ss1.end)
        self.assertEqual([], stap_r_ss1.deletions)
        self.assertEqual([], stap_r_ss1.insertions)

        self.assertEqual(3, len(scaf.domains))
        scaf_ss0 = scaf.domains[0]
        scaf_ss1 = scaf.domains[1]
        scaf_ss2 = scaf.domains[2]

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

    def test_inline_deletions_insertions__one_deletion(self) -> None:
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
        design = sc.Design(
            helices=[sc.Helix(max_offset=24, major_tick_distance=8)],
            strands=[sc.Strand([sc.Domain(0, True, 0, 8, deletions=[4])])],
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
        self.assertEqual(start, strand.domains[0].start)
        self.assertEqual(end, strand.domains[0].end)
        self.assertListEqual([], strand.domains[0].deletions)
        self.assertListEqual([], strand.domains[0].insertions)

    def test_inline_deletions_insertions__two_deletions(self) -> None:
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
        design = sc.Design(
            helices=[sc.Helix(max_offset=24, major_tick_distance=8)],
            strands=[sc.Strand([sc.Domain(0, True, 0, 8, deletions=[2, 4])])],
            grid=sc.square)
        design.inline_deletions_insertions()
        self.helix0_strand0_inlined_test(design, max_offset=22, major_ticks=[0, 6, 14, 22], start=0, end=6)

    def test_inline_deletions_insertions__one_insertion(self) -> None:
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
        design = sc.Design(
            helices=[sc.Helix(max_offset=24, major_tick_distance=8)],
            strands=[sc.Strand([sc.Domain(0, True, 0, 8, insertions=[(4, 1)])])],
            grid=sc.square)
        design.inline_deletions_insertions()
        self.helix0_strand0_inlined_test(design, max_offset=25, major_ticks=[0, 9, 17, 25], start=0, end=9)

    def test_inline_deletions_insertions__two_insertions(self) -> None:
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
        design = sc.Design(
            helices=[sc.Helix(max_offset=24, major_tick_distance=8)],
            strands=[sc.Strand([sc.Domain(0, True, 0, 8, insertions=[(2, 3), (4, 1)])])],
            grid=sc.square)
        design.inline_deletions_insertions()
        self.helix0_strand0_inlined_test(design, max_offset=28, major_ticks=[0, 12, 20, 28], start=0, end=12)

    def test_inline_deletions_insertions__one_deletion_one_insertion(self) -> None:
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
        design = sc.Design(
            helices=[sc.Helix(max_offset=24, major_tick_distance=8)],
            strands=[sc.Strand([sc.Domain(0, True, 0, 8, deletions=[4], insertions=[(2, 3)])])],
            grid=sc.square)
        design.inline_deletions_insertions()
        self.helix0_strand0_inlined_test(design, max_offset=26, major_ticks=[0, 10, 18, 26], start=0, end=10)

    def test_inline_deletions_insertions__one_deletion_right_of_major_tick(self) -> None:
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
        design = sc.Design(
            helices=[sc.Helix(max_offset=24, major_tick_distance=8)],
            strands=[sc.Strand([sc.Domain(0, True, 0, 12, deletions=[9])])],
            grid=sc.square)
        design.inline_deletions_insertions()
        self.helix0_strand0_inlined_test(design, max_offset=23, major_ticks=[0, 8, 15, 23], start=0, end=11)

    def test_inline_deletions_insertions__one_deletion_on_major_tick(self) -> None:
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
        design = sc.Design(
            helices=[sc.Helix(max_offset=24, major_tick_distance=8)],
            strands=[sc.Strand([sc.Domain(0, True, 0, 12, deletions=[8])])],
            grid=sc.square)
        design.inline_deletions_insertions()
        self.helix0_strand0_inlined_test(design, max_offset=23, major_ticks=[0, 8, 15, 23], start=0, end=11)

    def test_inline_deletions_insertions__one_deletion_left_of_major_tick(self) -> None:
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
        design = sc.Design(
            helices=[sc.Helix(max_offset=24, major_tick_distance=8)],
            strands=[sc.Strand([sc.Domain(0, True, 0, 12, deletions=[7])])],
            grid=sc.square)
        design.inline_deletions_insertions()
        self.helix0_strand0_inlined_test(design, max_offset=23, major_ticks=[0, 7, 15, 23], start=0, end=11)

    def test_inline_deletions_insertions__one_insertion_right_of_major_tick(self) -> None:
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
        design = sc.Design(
            helices=[sc.Helix(max_offset=24, major_tick_distance=8)],
            strands=[sc.Strand([sc.Domain(0, True, 0, 12, insertions=[(9, 1)])])],
            grid=sc.square)
        design.inline_deletions_insertions()
        self.helix0_strand0_inlined_test(design, max_offset=25, major_ticks=[0, 8, 17, 25], start=0, end=13)

    def test_inline_deletions_insertions__one_insertion_on_major_tick(self) -> None:
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
        design = sc.Design(
            helices=[sc.Helix(max_offset=24, major_tick_distance=8)],
            strands=[sc.Strand([sc.Domain(0, True, 0, 12, insertions=[(8, 1)])])],
            grid=sc.square)
        design.inline_deletions_insertions()
        self.helix0_strand0_inlined_test(design, max_offset=25, major_ticks=[0, 8, 17, 25], start=0, end=13)

    def test_inline_deletions_insertions__one_insertion_left_of_major_tick(self) -> None:
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
        design = sc.Design(
            helices=[sc.Helix(max_offset=24, major_tick_distance=8)],
            strands=[sc.Strand([sc.Domain(0, True, 0, 12, insertions=[(7, 1)])])],
            grid=sc.square)
        design.inline_deletions_insertions()
        self.helix0_strand0_inlined_test(design, max_offset=25, major_ticks=[0, 9, 17, 25], start=0, end=13)

    def test_inline_deletions_insertions__deletions_insertions_in_multiple_domains(self) -> None:
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
        design = sc.Design(
            helices=[sc.Helix(max_offset=24, major_tick_distance=8)],
            strands=[sc.Strand([sc.Domain(0, True, 0, 24, deletions=[19], insertions=[(5, 2), (11, 1)])])],
            grid=sc.square)
        design.inline_deletions_insertions()
        self.helix0_strand0_inlined_test(design, max_offset=26, major_ticks=[0, 10, 19, 26], start=0, end=26)

    def test_inline_deletions_insertions__deletions_insertions_in_multiple_domains_two_strands(self) -> None:
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
        design = sc.Design(
            helices=[sc.Helix(max_offset=24, major_tick_distance=8)],
            strands=[
                sc.Strand([sc.Domain(0, True, 0, 14, deletions=[2], insertions=[(5, 2), (10, 1)])]),
                sc.Strand([sc.Domain(0, True, 14, 24, deletions=[19])]),
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
        self.assertEqual(0, strand0.domains[0].start)
        self.assertEqual(16, strand0.domains[0].end)
        self.assertEqual(16, strand1.domains[0].start)
        self.assertEqual(25, strand1.domains[0].end)
        self.assertListEqual([], strand0.domains[0].deletions)
        self.assertListEqual([], strand0.domains[0].insertions)
        self.assertListEqual([], strand1.domains[0].deletions)
        self.assertListEqual([], strand1.domains[0].insertions)


class TestNickLigateAndCrossover(unittest.TestCase):
    """
    Tests add_nick(), ligate(), add_*_crossover(), and remove_*_crossover() methods on Design
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
    
    
    small_nicked_design
    0        8        16
    ACGTACGA AACCGGTA
0   [------> [------>
    <------] <------]
    TGCATGCT TTGGCCAT

    AAACCCGG TTTGGGCC
1   [------> [------>
    <------] <------]
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
    
    two_strand_design:    
    0        8        16
    
0   [------- ------->
    
   
1   <------- -------]        
    """

    def setUp(self) -> None:
        # small design
        strands_small_design = [
            sc.Strand([sc.Domain(0, True, 0, 16)]),
            sc.Strand([sc.Domain(0, False, 0, 16)]),
            sc.Strand([sc.Domain(1, True, 0, 16)]),
            sc.Strand([sc.Domain(1, False, 0, 16)]),
        ]
        self.small_design = sc.Design(strands=strands_small_design, grid=sc.square)
        self.small_design.assign_dna(strands_small_design[0], "ACGTACGA AACCGGTA")
        self.small_design.assign_dna(strands_small_design[2], "AAACCCGG TTTGGGCC")

        # small nicked design
        small_nicked_helices = [sc.Helix(max_offset=100) for _ in range(2)]
        self.small_nicked_design = sc.Design(helices=small_nicked_helices, grid=sc.square)
        # forward strands
        self.small_nicked_design.strand(0, 0).move(8)
        self.small_nicked_design.strand(0, 8).move(8)
        self.small_nicked_design.strand(1, 0).move(8)
        self.small_nicked_design.strand(1, 8).move(8)
        # reverse strands
        self.small_nicked_design.strand(0, 8).move(-8)
        self.small_nicked_design.strand(0, 16).move(-8)
        self.small_nicked_design.strand(1, 8).move(-8)
        self.small_nicked_design.strand(1, 16).move(-8)

        self.small_nicked_design.assign_dna(self.small_nicked_design.strands[0], "ACGTACGA")
        self.small_nicked_design.assign_dna(self.small_nicked_design.strands[1], "AACCGGTA")
        self.small_nicked_design.assign_dna(self.small_nicked_design.strands[2], "AAACCCGG")
        self.small_nicked_design.assign_dna(self.small_nicked_design.strands[3], "TTTGGGCC")

        # origami
        self.max_offset: int = 8 * 12
        scafs = []
        staps = []
        for helix in range(6):
            scaf_ss = sc.Domain(helix, helix % 2 == 0, 0, self.max_offset)
            stap_ss = sc.Domain(helix, helix % 2 == 1, 0, self.max_offset)
            scaf = sc.Strand([scaf_ss])
            stap = sc.Strand([stap_ss])
            scafs.append(scaf)
            staps.append(stap)
        self.origami: sc.Design = sc.Design(strands=scafs + staps, grid=sc.square)

        self.two_strand_design = sc.Design(helices=[sc.Helix(max_offset=100) for _ in range(2)])
        self.two_strand_design.strand(0, 0).move(16)
        self.two_strand_design.strand(1, 16).move(-16)

    def test_add_nick__twice_on_same_domain(self) -> None:
        """
        before
        0        8        16       24
    0   [------- -------- ------->

        after
        0        8        16       24
    0   [------> [------> [------>
        """
        design = sc.Design(strands=[
            sc.Strand([sc.Domain(0, True, 0, 24)]),
        ], grid=sc.square)
        design.add_nick(helix=0, offset=8, forward=True)
        design.add_nick(helix=0, offset=16, forward=True)
        self.assertEqual(3, len(design.strands))
        self.assertIn(sc.Strand([sc.Domain(0, True, 0, 8)]), design.strands)
        self.assertIn(sc.Strand([sc.Domain(0, True, 8, 16)]), design.strands)
        self.assertIn(sc.Strand([sc.Domain(0, True, 16, 24)]), design.strands)

    def test_ligate__twice_on_same_domain(self) -> None:
        """
        before
        0        8        16       24
    0   [------> [------> [------>

        after
        0        8        16       24
    0   [------- -------- ------->
        """
        design = sc.Design(helices=[sc.Helix(max_offset=24)], grid=sc.square)
        design.strand(0, 0).move(8)
        design.strand(0, 8).move(8)
        design.strand(0, 16).move(8)
        design.ligate(helix=0, offset=8, forward=True)
        design.ligate(helix=0, offset=16, forward=True)
        self.assertEqual(1, len(design.strands))
        self.assertIn(sc.Strand([sc.Domain(0, True, 0, 24)]), design.strands)

    def test_add_nick__small_design_no_nicks_added_yet(self) -> None:
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
        self.assertIn(sc.Strand([sc.Domain(0, True, 0, 16)]), self.small_design.strands)
        self.assertIn(sc.Strand([sc.Domain(0, False, 0, 16)]), self.small_design.strands)
        self.assertIn(sc.Strand([sc.Domain(1, True, 0, 16)]), self.small_design.strands)
        self.assertIn(sc.Strand([sc.Domain(1, False, 0, 16)]), self.small_design.strands)
        # DNA
        strand = strand_matching(self.small_design.strands, 0, True, 0, 16)
        self.assertEqual(remove_whitespace('ACGTACGA AACCGGTA'), strand.dna_sequence)
        strand = strand_matching(self.small_design.strands, 0, False, 0, 16)
        self.assertEqual(remove_whitespace('TACCGGTT TCGTACGT'), strand.dna_sequence)
        strand = strand_matching(self.small_design.strands, 1, True, 0, 16)
        self.assertEqual(remove_whitespace('AAACCCGG TTTGGGCC'), strand.dna_sequence)
        strand = strand_matching(self.small_design.strands, 1, False, 0, 16)
        self.assertEqual(remove_whitespace('GGCCCAAA CCGGGTTT'), strand.dna_sequence)

    def test_ligate__small_nicked_design_no_ligation_yet(self) -> None:
        """
        0        8        16
        ACGTACGA AACCGGTA
    0   [------> [------>
        <------] <------]
        TGCATGCT TTGGCCAT

        AAACCCGG TTTGGGCC
    1   [------> [------>
        <------] <------]
        TTTGGGCC AAACCCGG
        """
        design = self.small_nicked_design
        self.assertEqual(8, len(design.strands))
        self.assertIn(sc.Strand([sc.Domain(0, True, 0, 8)]), design.strands)
        self.assertIn(sc.Strand([sc.Domain(0, False, 0, 8)]), design.strands)
        self.assertIn(sc.Strand([sc.Domain(1, True, 0, 8)]), design.strands)
        self.assertIn(sc.Strand([sc.Domain(1, False, 0, 8)]), design.strands)
        self.assertIn(sc.Strand([sc.Domain(0, True, 8, 16)]), design.strands)
        self.assertIn(sc.Strand([sc.Domain(0, False, 8, 16)]), design.strands)
        self.assertIn(sc.Strand([sc.Domain(1, True, 8, 16)]), design.strands)
        self.assertIn(sc.Strand([sc.Domain(1, False, 8, 16)]), design.strands)
        # DNA
        strand = strand_matching(design.strands, 0, True, 0, 8)
        self.assertEqual(remove_whitespace('ACGTACGA'), strand.dna_sequence)
        strand = strand_matching(design.strands, 0, True, 8, 16)
        self.assertEqual(remove_whitespace('AACCGGTA'), strand.dna_sequence)

        strand = strand_matching(design.strands, 0, False, 0, 8)
        self.assertEqual(remove_whitespace('TCGTACGT'), strand.dna_sequence)
        strand = strand_matching(design.strands, 0, False, 8, 16)
        self.assertEqual(remove_whitespace('TACCGGTT'), strand.dna_sequence)

        strand = strand_matching(design.strands, 1, True, 0, 8)
        self.assertEqual(remove_whitespace('AAACCCGG'), strand.dna_sequence)
        strand = strand_matching(design.strands, 1, True, 8, 16)
        self.assertEqual(remove_whitespace('TTTGGGCC'), strand.dna_sequence)

        strand = strand_matching(design.strands, 1, False, 0, 8)
        self.assertEqual(remove_whitespace('CCGGGTTT'), strand.dna_sequence)
        strand = strand_matching(design.strands, 1, False, 8, 16)
        self.assertEqual(remove_whitespace('GGCCCAAA'), strand.dna_sequence)

    def test_add_nick__small_design_H0_forward(self) -> None:
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

        design = self.small_design
        design.add_nick(helix=0, offset=8, forward=True)
        self.assertEqual(5, len(design.strands))
        # two new Strands
        self.assertIn(sc.Strand([sc.Domain(0, True, 0, 8)]), design.strands)
        self.assertIn(sc.Strand([sc.Domain(0, True, 8, 16)]), design.strands)
        # existing Strands
        self.assertIn(sc.Strand([sc.Domain(0, False, 0, 16)]), design.strands)
        self.assertIn(sc.Strand([sc.Domain(1, True, 0, 16)]), design.strands)
        self.assertIn(sc.Strand([sc.Domain(1, False, 0, 16)]), design.strands)
        # DNA
        strand = strand_matching(design.strands, 0, True, 0, 8)
        self.assertEqual(remove_whitespace('ACGTACGA'), strand.dna_sequence)
        strand = strand_matching(design.strands, 0, True, 8, 16)
        self.assertEqual(remove_whitespace('AACCGGTA'), strand.dna_sequence)
        strand = strand_matching(design.strands, 0, False, 0, 16)
        self.assertEqual(remove_whitespace('TACCGGTT TCGTACGT'), strand.dna_sequence)
        strand = strand_matching(design.strands, 1, True, 0, 16)
        self.assertEqual(remove_whitespace('AAACCCGG TTTGGGCC'), strand.dna_sequence)
        strand = strand_matching(design.strands, 1, False, 0, 16)
        self.assertEqual(remove_whitespace('GGCCCAAA CCGGGTTT'), strand.dna_sequence)

    def test_ligate__small_nicked_design_ligate_all(self) -> None:
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
        design = self.small_nicked_design
        design.ligate(0, 8, True)
        design.ligate(0, 8, False)
        design.ligate(1, 8, True)
        design.ligate(1, 8, False)
        self.assertEqual(4, len(design.strands))
        self.assertIn(sc.Strand([sc.Domain(0, True, 0, 16)]), design.strands)
        self.assertIn(sc.Strand([sc.Domain(0, False, 0, 16)]), design.strands)
        self.assertIn(sc.Strand([sc.Domain(1, True, 0, 16)]), design.strands)
        self.assertIn(sc.Strand([sc.Domain(1, False, 0, 16)]), design.strands)
        # DNA
        strand = strand_matching(design.strands, 0, True, 0, 16)
        self.assertEqual(remove_whitespace('ACGTACGA AACCGGTA'), strand.dna_sequence)
        strand = strand_matching(design.strands, 0, False, 0, 16)
        self.assertEqual(remove_whitespace('TACCGGTT TCGTACGT'), strand.dna_sequence)
        strand = strand_matching(design.strands, 1, True, 0, 16)
        self.assertEqual(remove_whitespace('AAACCCGG TTTGGGCC'), strand.dna_sequence)
        strand = strand_matching(design.strands, 1, False, 0, 16)
        self.assertEqual(remove_whitespace('GGCCCAAA CCGGGTTT'), strand.dna_sequence)

    def test_add_nick__small_design_H0_reverse(self) -> None:
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
        design = self.small_design
        design.add_nick(helix=0, offset=8, forward=False)
        self.assertEqual(5, len(design.strands))
        # two new Strands
        self.assertIn(sc.Strand([sc.Domain(0, False, 0, 8)]), design.strands)
        self.assertIn(sc.Strand([sc.Domain(0, False, 8, 16)]), design.strands)
        # existing Strands
        self.assertIn(sc.Strand([sc.Domain(0, True, 0, 16)]), design.strands)
        self.assertIn(sc.Strand([sc.Domain(1, True, 0, 16)]), design.strands)
        self.assertIn(sc.Strand([sc.Domain(1, False, 0, 16)]), design.strands)
        # DNA
        strand = strand_matching(design.strands, 0, True, 0, 16)
        self.assertEqual(remove_whitespace('ACGTACGA AACCGGTA'), strand.dna_sequence)
        strand = strand_matching(design.strands, 0, False, 8, 16)
        self.assertEqual(remove_whitespace('TACCGGTT'), strand.dna_sequence)
        strand = strand_matching(design.strands, 0, False, 0, 8)
        self.assertEqual(remove_whitespace('TCGTACGT'), strand.dna_sequence)
        strand = strand_matching(design.strands, 1, True, 0, 16)
        self.assertEqual(remove_whitespace('AAACCCGG TTTGGGCC'), strand.dna_sequence)
        strand = strand_matching(design.strands, 1, False, 0, 16)
        self.assertEqual(remove_whitespace('GGCCCAAA CCGGGTTT'), strand.dna_sequence)

    def test_add_nick__small_design_H1_forward(self) -> None:
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
        design = self.small_design
        design.add_nick(helix=1, offset=8, forward=True)
        self.assertEqual(5, len(design.strands))
        # two new Strands
        self.assertIn(sc.Strand([sc.Domain(1, True, 0, 8)]), design.strands)
        self.assertIn(sc.Strand([sc.Domain(1, True, 8, 16)]), design.strands)
        # existing Strands
        self.assertIn(sc.Strand([sc.Domain(1, False, 0, 16)]), design.strands)
        self.assertIn(sc.Strand([sc.Domain(0, True, 0, 16)]), design.strands)
        self.assertIn(sc.Strand([sc.Domain(0, False, 0, 16)]), design.strands)
        # DNA
        strand = strand_matching(design.strands, 0, True, 0, 16)
        self.assertEqual(remove_whitespace('ACGTACGA AACCGGTA'), strand.dna_sequence)
        strand = strand_matching(design.strands, 0, False, 0, 16)
        self.assertEqual(remove_whitespace('TACCGGTT TCGTACGT'), strand.dna_sequence)
        strand = strand_matching(design.strands, 1, True, 0, 8)
        self.assertEqual(remove_whitespace('AAACCCGG'), strand.dna_sequence)
        strand = strand_matching(design.strands, 1, True, 8, 16)
        self.assertEqual(remove_whitespace('TTTGGGCC'), strand.dna_sequence)
        strand = strand_matching(design.strands, 1, False, 0, 16)
        self.assertEqual(remove_whitespace('GGCCCAAA CCGGGTTT'), strand.dna_sequence)

    def test_add_nick__small_design_H1_reverse(self) -> None:
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
        design = self.small_design
        design.add_nick(helix=1, offset=8, forward=False)
        self.assertEqual(5, len(design.strands))
        # two new Strands
        self.assertIn(sc.Strand([sc.Domain(1, False, 0, 8)]), design.strands)
        self.assertIn(sc.Strand([sc.Domain(1, False, 8, 16)]), design.strands)
        # existing Strands
        self.assertIn(sc.Strand([sc.Domain(1, True, 0, 16)]), design.strands)
        self.assertIn(sc.Strand([sc.Domain(0, True, 0, 16)]), design.strands)
        self.assertIn(sc.Strand([sc.Domain(0, False, 0, 16)]), design.strands)
        # DNA
        strand = strand_matching(design.strands, 0, True, 0, 16)
        self.assertEqual(remove_whitespace('ACGTACGA AACCGGTA'), strand.dna_sequence)
        strand = strand_matching(design.strands, 0, False, 0, 16)
        self.assertEqual(remove_whitespace('TACCGGTT TCGTACGT'), strand.dna_sequence)
        strand = strand_matching(design.strands, 1, True, 0, 16)
        self.assertEqual(remove_whitespace('AAACCCGG TTTGGGCC'), strand.dna_sequence)
        strand = strand_matching(design.strands, 1, False, 0, 8)
        self.assertEqual(remove_whitespace('CCGGGTTT'), strand.dna_sequence)
        strand = strand_matching(design.strands, 1, False, 8, 16)
        self.assertEqual(remove_whitespace('GGCCCAAA'), strand.dna_sequence)

    def test_add_full_crossover__small_design_H0_forward(self) -> None:
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
        design = self.small_design
        design.add_full_crossover(helix=0, helix2=1, offset=8, forward=True)
        self.assertEqual(4, len(design.strands))
        # two new Strands
        self.assertIn(sc.Strand([
            sc.Domain(0, True, 0, 8),
            sc.Domain(1, False, 0, 8),
        ]), design.strands)
        self.assertIn(sc.Strand([
            sc.Domain(1, False, 8, 16),
            sc.Domain(0, True, 8, 16),
        ]), design.strands)
        # existing Strands
        self.assertIn(sc.Strand([sc.Domain(0, False, 0, 16)]), design.strands)
        self.assertIn(sc.Strand([sc.Domain(1, True, 0, 16)]), design.strands)
        # DNA
        strand = strand_matching(design.strands, 0, False, 0, 16)
        self.assertEqual(remove_whitespace('TACCGGTT TCGTACGT'), strand.dna_sequence)
        strand = strand_matching(design.strands, 1, True, 0, 16)
        self.assertEqual(remove_whitespace('AAACCCGG TTTGGGCC'), strand.dna_sequence)
        strand = strand_matching(design.strands, 0, True, 0, 8)
        self.assertEqual(remove_whitespace('ACGTACGA CCGGGTTT'), strand.dna_sequence)
        strand = strand_matching(design.strands, 1, False, 8, 16)
        self.assertEqual(remove_whitespace('GGCCCAAA AACCGGTA'), strand.dna_sequence)

    def test_add_full_crossover__small_design_H0_reverse(self) -> None:
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
        design = self.small_design
        design.add_full_crossover(helix=0, helix2=1, offset=8, forward=False)
        self.assertEqual(4, len(design.strands))
        # two new Strands
        self.assertIn(sc.Strand([
            sc.Domain(1, True, 0, 8),
            sc.Domain(0, False, 0, 8),
        ]), design.strands)
        self.assertIn(sc.Strand([
            sc.Domain(0, False, 8, 16),
            sc.Domain(1, True, 8, 16),
        ]), design.strands)
        # existing Strands
        self.assertIn(sc.Strand([sc.Domain(0, True, 0, 16)]), design.strands)
        self.assertIn(sc.Strand([sc.Domain(1, False, 0, 16)]), design.strands)
        # DNA
        strand = strand_matching(design.strands, 0, True, 0, 16)
        self.assertEqual(remove_whitespace('ACGTACGA AACCGGTA'), strand.dna_sequence)
        strand = strand_matching(design.strands, 0, False, 8, 16)
        self.assertEqual(remove_whitespace('TACCGGTT TTTGGGCC'), strand.dna_sequence)
        strand = strand_matching(design.strands, 1, True, 0, 8)
        self.assertEqual(remove_whitespace('AAACCCGG TCGTACGT'), strand.dna_sequence)
        strand = strand_matching(design.strands, 1, False, 0, 16)
        self.assertEqual(remove_whitespace('GGCCCAAA CCGGGTTT'), strand.dna_sequence)

    def test_add_half_crossover__small_design_H0_reverse_8(self) -> None:
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
        design = self.small_design
        design.add_nick(helix=0, offset=8, forward=False)
        design.add_nick(helix=1, offset=8, forward=True)
        design.add_half_crossover(helix=0, helix2=1, offset=8, forward=False)
        self.assertEqual(5, len(design.strands))
        # three new Strands
        self.assertIn(sc.Strand([
            sc.Domain(1, True, 0, 8),
        ]), design.strands)
        self.assertIn(sc.Strand([
            sc.Domain(0, False, 0, 8),
        ]), design.strands)
        self.assertIn(sc.Strand([
            sc.Domain(0, False, 8, 16),
            sc.Domain(1, True, 8, 16),
        ]), design.strands)
        # existing Strands
        self.assertIn(sc.Strand([sc.Domain(0, True, 0, 16)]), design.strands)
        self.assertIn(sc.Strand([sc.Domain(1, False, 0, 16)]), design.strands)
        # DNA
        strand = strand_matching(design.strands, 0, True, 0, 16)
        self.assertEqual(remove_whitespace('ACGTACGA AACCGGTA'), strand.dna_sequence)
        strand = strand_matching(design.strands, 0, False, 8, 16)
        self.assertEqual(remove_whitespace('TACCGGTT TTTGGGCC'), strand.dna_sequence)
        strand = strand_matching(design.strands, 0, False, 0, 8)
        self.assertEqual(remove_whitespace('TCGTACGT'), strand.dna_sequence)
        strand = strand_matching(design.strands, 1, True, 0, 8)
        self.assertEqual(remove_whitespace('AAACCCGG'), strand.dna_sequence)
        strand = strand_matching(design.strands, 1, False, 0, 16)
        self.assertEqual(remove_whitespace('GGCCCAAA CCGGGTTT'), strand.dna_sequence)

    def test_add_half_crossover__both_scaffolds_preserved(self) -> None:
        """
        0       8        16

    0   [------- -------> scaffold


    1   <------- -------] scaffold

        """
        design: sc.Design = self.two_strand_design
        design.strands[0].set_scaffold(True)
        design.strands[1].set_scaffold(True)
        design.add_half_crossover(helix=0, helix2=1, offset=15, forward=True)
        self.assertEqual(1, len(design.strands))
        self.assertTrue(design.strands[0].is_scaffold)
        self.assertIsNotNone(design.scaffold)

    def test_add_half_crossover__first_scaffold_preserved(self) -> None:
        """
        0       8        16

    0   [------- -------> scaffold


    1   <------- -------]

        """
        design: sc.Design = self.two_strand_design
        design.strands[0].set_scaffold(True)
        design.add_half_crossover(helix=0, helix2=1, offset=15, forward=True)
        self.assertEqual(1, len(design.strands))
        self.assertTrue(design.strands[0].is_scaffold)
        self.assertIsNotNone(design.scaffold)

    def test_add_half_crossover__second_scaffold_preserved(self) -> None:
        """
        0       8        16

    0   [------- ------->


    1   <------- -------] scaffold

        """
        design: sc.Design = self.two_strand_design
        design.strands[1].set_scaffold(True)
        design.add_half_crossover(helix=0, helix2=1, offset=15, forward=True)
        self.assertEqual(1, len(design.strands))
        self.assertTrue(design.strands[0].is_scaffold)
        self.assertIsNotNone(design.scaffold)

    def test_add_half_crossover__neither_scaffold_preserved(self) -> None:
        """
        0       8        16

    0   [------- ------->


    1   <------- -------]

        """
        design: sc.Design = self.two_strand_design
        design.add_half_crossover(helix=0, helix2=1, offset=15, forward=True)
        self.assertEqual(1, len(design.strands))
        self.assertFalse(design.strands[0].is_scaffold)
        self.assertIsNone(design.scaffold)

    def test_add_half_crossover__small_design_H0_reverse_0(self) -> None:
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
        design = self.small_design
        design.add_half_crossover(helix=0, helix2=1, offset=0, forward=False)
        self.assertEqual(3, len(design.strands))
        # one new Strand
        self.assertIn(sc.Strand([
            sc.Domain(0, False, 0, 16),
            sc.Domain(1, True, 0, 16),
        ]), design.strands)
        # existing Strands
        self.assertIn(sc.Strand([sc.Domain(0, True, 0, 16)]), design.strands)
        self.assertIn(sc.Strand([sc.Domain(1, False, 0, 16)]), design.strands)
        # DNA
        strand = strand_matching(design.strands, 0, True, 0, 16)
        self.assertEqual(remove_whitespace('ACGTACGA AACCGGTA'), strand.dna_sequence)
        strand = strand_matching(design.strands, 0, False, 0, 16)
        self.assertEqual(remove_whitespace('TACCGGTT TCGTACGT AAACCCGG TTTGGGCC'), strand.dna_sequence)
        strand = strand_matching(design.strands, 1, False, 0, 16)
        self.assertEqual(remove_whitespace('GGCCCAAA CCGGGTTT'), strand.dna_sequence)

    def test_add_half_crossover__small_design_H0_reverse_15(self) -> None:
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
        design = self.small_design
        design.add_half_crossover(helix=0, helix2=1, offset=15, forward=False)
        self.assertEqual(3, len(design.strands))
        # one new Strand
        self.assertIn(sc.Strand([
            sc.Domain(1, True, 0, 16),
            sc.Domain(0, False, 0, 16),
        ]), design.strands)
        # existing Strands
        self.assertIn(sc.Strand([sc.Domain(0, True, 0, 16)]), design.strands)
        self.assertIn(sc.Strand([sc.Domain(1, False, 0, 16)]), design.strands)
        # DNA
        strand = strand_matching(design.strands, 0, True, 0, 16)
        self.assertEqual(remove_whitespace('ACGTACGA AACCGGTA'), strand.dna_sequence)
        strand = strand_matching(design.strands, 1, True, 0, 16)
        self.assertEqual(remove_whitespace('AAACCCGG TTTGGGCC TACCGGTT TCGTACGT'), strand.dna_sequence)
        strand = strand_matching(design.strands, 1, False, 0, 16)
        self.assertEqual(remove_whitespace('GGCCCAAA CCGGGTTT'), strand.dna_sequence)

    def test_add_half_crossover__small_design_illegal(self) -> None:
        """
        0        8        16
    0   [------- ------->
        <------- -------- ?
                          |
    1   [------- -------- ?
        <------- -------]
        """
        with self.assertRaises(sc.IllegalDesignError):
            self.small_design.add_half_crossover(helix=0, helix2=1, offset=16, forward=False)

    def test_add_full_crossover__small_design_illegal(self) -> None:
        """
        0        8        16
    0   [------- ------->
        <------- -------+ ?
                        | |
    1   [------- -------+ ?
        <------- -------]
        """
        with self.assertRaises(sc.IllegalDesignError):
            self.small_design.add_full_crossover(helix=0, helix2=1, offset=16, forward=False)

    def test_add_full_crossover__small_design_illegal_only_one_helix_has_domain(self) -> None:
        """
        0        8        16
    0   [------- ------->
        <------+ +------]
               | |
    1   [--->  ? ?
        <---]
        """
        design = sc.Design(strands=[
            sc.Strand([sc.Domain(0, True, 0, 16)]),
            sc.Strand([sc.Domain(0, False, 0, 16)]),
            sc.Strand([sc.Domain(1, True, 0, 5)]),
            sc.Strand([sc.Domain(1, False, 0, 5)]),
        ], grid=sc.square)
        with self.assertRaises(sc.IllegalDesignError):
            design.add_full_crossover(helix=0, helix2=1, offset=10, forward=False)

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

    @staticmethod
    def add_nicks(design: sc.Design):
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

    def test_add_nick__6_helix_rectangle(self) -> None:
        self.add_nicks(self.origami)
        self.assertEqual(25, len(self.origami.strands))
        for helix in range(0, len(self.origami.helices), 2):
            # even helix
            self.assertIn(sc.Strand([sc.Domain(helix, True, 0, 96)]), self.origami.strands)
            self.assertIn(sc.Strand([sc.Domain(helix, False, 0, 40)]), self.origami.strands)
            self.assertIn(sc.Strand([sc.Domain(helix, False, 40, 72)]), self.origami.strands)
            self.assertIn(sc.Strand([sc.Domain(helix, False, 72, 96)]), self.origami.strands)
            # odd helix
            if helix + 1 < len(self.origami.helices) - 1:
                self.assertIn(sc.Strand([sc.Domain(helix + 1, False, 0, 96)]), self.origami.strands)
            else:
                # nick in scaffold on bottom helix
                self.assertIn(sc.Strand([sc.Domain(helix + 1, False, 0, 48)]), self.origami.strands)
                self.assertIn(sc.Strand([sc.Domain(helix + 1, False, 48, 96)]), self.origami.strands)
            self.assertIn(sc.Strand([sc.Domain(helix + 1, True, 0, 24)]), self.origami.strands)
            self.assertIn(sc.Strand([sc.Domain(helix + 1, True, 24, 56)]), self.origami.strands)
            self.assertIn(sc.Strand([sc.Domain(helix + 1, True, 56, 96)]), self.origami.strands)

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

    @staticmethod
    def add_crossovers_after_nicks(design: sc.Design):
        # scaffold seam crossovers
        design.add_full_crossover(helix=1, helix2=2, offset=48, forward=False)
        design.add_full_crossover(helix=3, helix2=4, offset=48, forward=False)

        # staple crossovers
        design.add_full_crossover(helix=0, helix2=1, offset=16, forward=False)
        design.add_full_crossover(helix=0, helix2=1, offset=80, forward=False)
        design.add_full_crossover(helix=1, helix2=2, offset=32, forward=True)
        design.add_full_crossover(helix=1, helix2=2, offset=64, forward=True)
        design.add_full_crossover(helix=2, helix2=3, offset=16, forward=False)
        design.add_full_crossover(helix=2, helix2=3, offset=80, forward=False)
        design.add_full_crossover(helix=3, helix2=4, offset=32, forward=True)
        design.add_full_crossover(helix=3, helix2=4, offset=64, forward=True)
        design.add_full_crossover(helix=4, helix2=5, offset=16, forward=False)
        design.add_full_crossover(helix=4, helix2=5, offset=80, forward=False)

        # The left and right edge crossovers need to be added last to ensure the Strands remain
        # non-circular during all intermediate stages.

        # scaffold left crossovers
        design.add_half_crossover(helix=0, helix2=1, offset=0, forward=True)
        design.add_half_crossover(helix=2, helix2=3, offset=0, forward=True)
        design.add_half_crossover(helix=4, helix2=5, offset=0, forward=True)

        # scaffold right crossovers
        design.add_half_crossover(helix=0, helix2=1, offset=95, forward=True)
        design.add_half_crossover(helix=2, helix2=3, offset=95, forward=True)
        design.add_half_crossover(helix=4, helix2=5, offset=95, forward=True)

    def test_add_nick_then_add_crossovers__6_helix_rectangle(self) -> None:
        self.add_nicks(self.origami)
        self.add_crossovers_after_nicks(self.origami)

        self.assertEqual(19, len(self.origami.strands))

        # staples left edge
        # {"helix": 1, "forward": true, "start": 0, "end": 16},
        # {"helix": 0, "forward": false, "start": 0, "end": 16}
        stap = sc.Strand([
            sc.Domain(1, True, 0, 16),
            sc.Domain(0, False, 0, 16),
        ])
        self.assertIn(stap, self.origami.strands)

        # {"helix": 3, "forward": true, "start": 0, "end": 16},
        # {"helix": 2, "forward": false, "start": 0, "end": 16}
        stap = sc.Strand([
            sc.Domain(3, True, 0, 16),
            sc.Domain(2, False, 0, 16),
        ])
        self.assertIn(stap, self.origami.strands)

        # {"helix": 5, "forward": true, "start": 0, "end": 16},
        # {"helix": 4, "forward": false, "start": 0, "end": 16}
        stap = sc.Strand([
            sc.Domain(5, True, 0, 16),
            sc.Domain(4, False, 0, 16),
        ])
        self.assertIn(stap, self.origami.strands)

        # staples right edge
        # {"helix": 0, "forward": false, "start": 80, "end": 96},
        # {"helix": 1, "forward": true, "start": 80, "end": 96}
        stap = sc.Strand([
            sc.Domain(0, False, 80, 96),
            sc.Domain(1, True, 80, 96),
        ])
        self.assertIn(stap, self.origami.strands)

        # {"helix": 2, "forward": false, "start": 80, "end": 96},
        # {"helix": 3, "forward": true, "start": 80, "end": 96}
        stap = sc.Strand([
            sc.Domain(2, False, 80, 96),
            sc.Domain(3, True, 80, 96),
        ])
        self.assertIn(stap, self.origami.strands)

        # {"helix": 4, "forward": false, "start": 80, "end": 96},
        # {"helix": 5, "forward": true, "start": 80, "end": 96}
        stap = sc.Strand([
            sc.Domain(4, False, 80, 96),
            sc.Domain(5, True, 80, 96),
        ])
        self.assertIn(stap, self.origami.strands)

        # staples remainder
        # {"helix": 0, "forward": false, "start": 40, "end": 72}
        stap = sc.Strand([sc.Domain(0, False, 40, 72)])
        self.assertIn(stap, self.origami.strands)

        # {"helix": 2, "forward": false, "start": 32, "end": 40},
        # {"helix": 1, "forward": true, "start": 32, "end": 56}
        stap = sc.Strand([
            sc.Domain(2, False, 32, 40),
            sc.Domain(1, True, 32, 56),
        ])
        self.assertIn(stap, self.origami.strands)

        # {"helix": 1, "forward": true, "start": 56, "end": 64},
        # {"helix": 2, "forward": false, "start": 40, "end": 64}
        stap = sc.Strand([
            sc.Domain(1, True, 56, 64),
            sc.Domain(2, False, 40, 64),
        ])
        self.assertIn(stap, self.origami.strands)

        # {"helix": 4, "forward": false, "start": 32, "end": 40},
        # {"helix": 3, "forward": true, "start": 32, "end": 56}
        stap = sc.Strand([
            sc.Domain(4, False, 32, 40),
            sc.Domain(3, True, 32, 56),
        ])
        self.assertIn(stap, self.origami.strands)

        # {"helix": 3, "forward": true, "start": 56, "end": 64},
        # {"helix": 4, "forward": false, "start": 40, "end": 64}
        stap = sc.Strand([
            sc.Domain(3, True, 56, 64),
            sc.Domain(4, False, 40, 64),
        ])
        self.assertIn(stap, self.origami.strands)

        # {"helix": 5, "forward": true, "start": 24, "end": 56}
        stap = sc.Strand([sc.Domain(5, True, 24, 56)])
        self.assertIn(stap, self.origami.strands)

        # {"helix": 0, "forward": false, "start": 16, "end": 40},
        # {"helix": 1, "forward": true, "start": 16, "end": 24}
        stap = sc.Strand([
            sc.Domain(0, False, 16, 40),
            sc.Domain(1, True, 16, 24),
        ])
        self.assertIn(stap, self.origami.strands)

        # {"helix": 1, "forward": true, "start": 24, "end": 32},
        # {"helix": 2, "forward": false, "start": 16, "end": 32},
        # {"helix": 3, "forward": true, "start": 16, "end": 24}
        stap = sc.Strand([
            sc.Domain(1, True, 24, 32),
            sc.Domain(2, False, 16, 32),
            sc.Domain(3, True, 16, 24),
        ])
        self.assertIn(stap, self.origami.strands)

        # {"helix": 3, "forward": true, "start": 24, "end": 32},
        # {"helix": 4, "forward": false, "start": 16, "end": 32},
        # {"helix": 5, "forward": true, "start": 16, "end": 24}
        stap = sc.Strand([
            sc.Domain(3, True, 24, 32),
            sc.Domain(4, False, 16, 32),
            sc.Domain(5, True, 16, 24),
        ])
        self.assertIn(stap, self.origami.strands)

        # {"helix": 5, "forward": true, "start": 56, "end": 80},
        # {"helix": 4, "forward": false, "start": 72, "end": 80}
        stap = sc.Strand([
            sc.Domain(5, True, 56, 80),
            sc.Domain(4, False, 72, 80),
        ])
        self.assertIn(stap, self.origami.strands)

        # {"helix": 2, "forward": false, "start": 64, "end": 72},
        # {"helix": 1, "forward": true, "start": 64, "end": 80},
        # {"helix": 0, "forward": false, "start": 72, "end": 80}
        stap = sc.Strand([
            sc.Domain(2, False, 64, 72),
            sc.Domain(1, True, 64, 80),
            sc.Domain(0, False, 72, 80),
        ])
        self.assertIn(stap, self.origami.strands)

        # {"helix": 4, "forward": false, "start": 64, "end": 72},
        # {"helix": 3, "forward": true, "start": 64, "end": 80},
        # {"helix": 2, "forward": false, "start": 72, "end": 80}
        stap = sc.Strand([
            sc.Domain(4, False, 64, 72),
            sc.Domain(3, True, 64, 80),
            sc.Domain(2, False, 72, 80),
        ])
        self.assertIn(stap, self.origami.strands)

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
            sc.Domain(5, False, 0, 48),
            sc.Domain(4, True, 0, 48),
            sc.Domain(3, False, 0, 48),
            sc.Domain(2, True, 0, 48),
            sc.Domain(1, False, 0, 48),
            sc.Domain(0, True, 0, 96),
            sc.Domain(1, False, 48, 96),
            sc.Domain(2, True, 48, 96),
            sc.Domain(3, False, 48, 96),
            sc.Domain(4, True, 48, 96),
            sc.Domain(5, False, 48, 96),
        ])
        self.assertIn(scaf, self.origami.strands)


class TestAutocalculatedData(unittest.TestCase):

    def test_helix_min_max_offsets_illegal_explicitly_specified(self) -> None:
        helices = [sc.Helix(min_offset=5, max_offset=5)]
        with self.assertRaises(sc.IllegalDesignError):
            sc.Design(helices=helices, strands=[], grid=sc.square)

    def test_helix_min_max_offsets_illegal_autocalculated(self) -> None:
        helices = [sc.Helix(min_offset=5)]
        ss = sc.Domain(0, True, 0, 4)
        strand = sc.Strand([ss])
        with self.assertRaises(sc.IllegalDesignError):
            sc.Design(helices=helices, strands=[strand], grid=sc.square)

    def test_helix_min_max_offsets(self) -> None:
        helices = [sc.Helix(), sc.Helix(min_offset=-5), sc.Helix(max_offset=5),
                   sc.Helix(min_offset=5, max_offset=10)]
        ss_0 = sc.Domain(helix=0, forward=True, start=20, end=25)
        ss_1 = sc.Domain(helix=1, forward=False, start=-5, end=30)
        ss_2 = sc.Domain(helix=2, forward=True, start=0, end=5)
        ss_3 = sc.Domain(helix=3, forward=False, start=5, end=10)
        strand = sc.Strand([ss_0, ss_1, ss_2, ss_3])
        design = sc.Design(helices=helices, strands=[strand], grid=sc.square)
        self.assertEqual(0, design.helices[0].min_offset)
        self.assertEqual(25, design.helices[0].max_offset)
        self.assertEqual(-5, design.helices[1].min_offset)
        self.assertEqual(30, design.helices[1].max_offset)
        self.assertEqual(0, design.helices[2].min_offset)
        self.assertEqual(5, design.helices[2].max_offset)
        self.assertEqual(5, design.helices[3].min_offset)
        self.assertEqual(10, design.helices[3].max_offset)

    def test_helix_max_offset(self) -> None:
        helices = [sc.Helix(), sc.Helix(max_offset=8), sc.Helix()]
        ss_0 = sc.Domain(helix=0, forward=True, start=5, end=10)
        ss_1 = sc.Domain(helix=1, forward=False, start=2, end=6)
        ss_2 = sc.Domain(helix=2, forward=True, start=0, end=5)
        strand = sc.Strand([ss_0, ss_1, ss_2])
        design = sc.Design(helices=helices, strands=[strand], grid=sc.square)
        self.assertEqual(10, design.helices[0].max_offset)
        self.assertEqual(8, design.helices[1].max_offset)
        self.assertEqual(5, design.helices[2].max_offset)


class TestSetHelixIdx(unittest.TestCase):
    def test_set_helix_idx(self) -> None:
        helices = [sc.Helix(max_offset=20), sc.Helix(max_offset=20), sc.Helix(max_offset=20)]
        ss_0 = sc.Domain(helix=0, forward=True, start=0, end=6)
        ss_1 = sc.Domain(helix=1, forward=True, start=0, end=7)
        ss_2 = sc.Domain(helix=2, forward=True, start=0, end=8)
        strand = sc.Strand([ss_0, ss_1, ss_2])
        design = sc.Design(helices=helices, strands=[strand], grid=sc.square)
        design.set_helix_idx(2, 5)
        design.set_helix_idx(1, 3)
        design.set_helix_idx(0, 2)
        self.assertEqual(3, len(design.helices))
        self.assertTrue(2 in design.helices)
        self.assertTrue(3 in design.helices)
        self.assertTrue(5 in design.helices)
        h2 = design.helices[2]
        h3 = design.helices[3]
        h5 = design.helices[5]
        self.assertEqual(2, h2.idx)
        self.assertEqual(3, h3.idx)
        self.assertEqual(5, h5.idx)

        ss2 = design.domain_at(2, 0, True)
        ss3 = design.domain_at(3, 0, True)
        ss5 = design.domain_at(5, 0, True)
        self.assertEqual(2, ss2.helix)
        self.assertEqual(3, ss3.helix)
        self.assertEqual(5, ss5.helix)
        self.assertEqual(6, ss2.end)
        self.assertEqual(7, ss3.end)
        self.assertEqual(8, ss5.end)


class TestDesignPitchYawRollOfHelix(unittest.TestCase):
    def setUp(self) -> None:
        n = 'north'
        helix = sc.Helix(max_offset=12, group=n, grid_position=(1, 2), roll=4)

        group_north = sc.HelixGroup(position=sc.Position3D(x=0, y=-200, z=0), grid=sc.square, pitch=12,
                                    yaw=40, roll=32)
        self.design = sc.Design(helices=[helix], groups={n: group_north}, strands=[])
        self.helix = helix

    def test_design_pitch_of_helix(self) -> None:
        self.assertEqual(12, self.design.pitch_of_helix(self.helix))

    def test_design_yaw_of_helix(self) -> None:
        self.assertEqual(40, self.design.yaw_of_helix(self.helix))

    def test_design_roll_of_helix(self) -> None:
        self.assertEqual(36, self.design.roll_of_helix(self.helix))


class TestHelixGroups(unittest.TestCase):
    def setUp(self) -> None:
        n = 'north'
        e = 'east'
        s = 'south'
        w = 'west'
        helices = [
            sc.Helix(max_offset=20, group=n, grid_position=(1, 1)),  # 0
            sc.Helix(max_offset=21, group=n, grid_position=(0, 1)),  # 1
            sc.Helix(max_offset=19, group=n, grid_position=(0, 2)),  # 2
            sc.Helix(max_offset=18, group=n, grid_position=(1, 2)),  # 3
            sc.Helix(max_offset=17, group=n, grid_position=(2, 2)),  # 4
            sc.Helix(max_offset=16, group=n, grid_position=(2, 1)),  # 5
            sc.Helix(max_offset=24, group=s),  # 6
            sc.Helix(max_offset=25, group=s),  # 7
            sc.Helix(max_offset=26, group=w, position=sc.Position3D(x=0, y=0, z=0)),  # 8
            sc.Helix(max_offset=27, group=w, position=sc.Position3D(x=0, y=2.5, z=0)),  # 9
            sc.Helix(idx=13, max_offset=22, group=e),  # 13
            sc.Helix(idx=15, max_offset=23, group=e),  # 15
        ]
        group_north = sc.HelixGroup(position=sc.Position3D(x=0, y=-200, z=0), grid=sc.honeycomb)
        group_south = sc.HelixGroup(position=sc.Position3D(x=0, y=70, z=0), helices_view_order=[7, 6],
                                    grid=sc.square)
        group_east = sc.HelixGroup(position=sc.Position3D(x=0, y=0, z=100), pitch=45, grid=sc.square)
        group_west = sc.HelixGroup()
        groups = {
            n: group_north,
            e: group_east,
            s: group_south,
            w: group_west,
        }
        self.design = sc.Design(helices=helices, groups=groups, strands=[])
        self.n = n
        self.e = e
        self.s = s
        self.w = w

    def test_helix_groups(self) -> None:
        self._asserts_for_fixture(self.design)

    def test_helix_groups_to_from_JSON(self) -> None:
        n = self.n
        e = self.e
        s = self.s
        w = self.w
        design_json_str = self.design.to_json()

        design_json_map = json.loads(design_json_str)
        groups_map = design_json_map[sc.groups_key]
        group_n = groups_map[n]
        group_e = groups_map[e]
        group_s = groups_map[s]
        group_w = groups_map[w]

        pos_n = group_n[sc.position_key]
        self.assertAlmostEqual(0, pos_n['x'])
        self.assertAlmostEqual(-200, pos_n['y'])
        self.assertAlmostEqual(0, pos_n['z'])

        pos_s = group_e[sc.position_key]
        self.assertAlmostEqual(0, pos_s['x'])
        self.assertAlmostEqual(0, pos_s['y'])
        self.assertAlmostEqual(100, pos_s['z'])

        pos_w = group_s[sc.position_key]
        self.assertAlmostEqual(0, pos_w['x'])
        self.assertAlmostEqual(70, pos_w['y'])
        self.assertAlmostEqual(0, pos_w['z'])

        pos_e = group_w[sc.position_key]
        self.assertAlmostEqual(0, pos_e['x'])
        self.assertAlmostEqual(0, pos_e['y'])
        self.assertAlmostEqual(0, pos_e['z'])

        helices_map = design_json_map[sc.helices_key]
        self.assertEqual(12, len(helices_map))
        helix0_map = helices_map[0]
        helix1_map = helices_map[1]
        helix2_map = helices_map[2]
        helix3_map = helices_map[3]
        helix4_map = helices_map[4]
        helix5_map = helices_map[5]
        helix6_map = helices_map[6]
        helix7_map = helices_map[7]
        helix8_map = helices_map[8]
        helix9_map = helices_map[9]
        helix13_map = helices_map[10]
        helix15_map = helices_map[11]

        self.assertEqual(n, helix0_map[sc.group_key])
        self.assertEqual(n, helix1_map[sc.group_key])
        self.assertEqual(n, helix2_map[sc.group_key])
        self.assertEqual(n, helix3_map[sc.group_key])
        self.assertEqual(n, helix4_map[sc.group_key])
        self.assertEqual(n, helix5_map[sc.group_key])
        self.assertEqual(s, helix6_map[sc.group_key])
        self.assertEqual(s, helix7_map[sc.group_key])
        self.assertEqual(w, helix8_map[sc.group_key])
        self.assertEqual(w, helix9_map[sc.group_key])
        self.assertEqual(e, helix13_map[sc.group_key])
        self.assertEqual(e, helix15_map[sc.group_key])

        self.assertEqual(0, helix0_map[sc.idx_on_helix_key])
        self.assertEqual(1, helix1_map[sc.idx_on_helix_key])
        self.assertEqual(2, helix2_map[sc.idx_on_helix_key])
        self.assertEqual(3, helix3_map[sc.idx_on_helix_key])
        self.assertEqual(4, helix4_map[sc.idx_on_helix_key])
        self.assertEqual(5, helix5_map[sc.idx_on_helix_key])
        self.assertEqual(6, helix6_map[sc.idx_on_helix_key])
        self.assertEqual(7, helix7_map[sc.idx_on_helix_key])
        self.assertEqual(8, helix8_map[sc.idx_on_helix_key])
        self.assertEqual(9, helix9_map[sc.idx_on_helix_key])
        self.assertEqual(13, helix13_map[sc.idx_on_helix_key])
        self.assertEqual(15, helix15_map[sc.idx_on_helix_key])

        design_from_json = sc.Design.from_scadnano_json_str(design_json_str)
        self._asserts_for_fixture(design_from_json)

    def test_helix_groups_fail_nonexistent(self) -> None:
        helices = [
            sc.Helix(max_offset=20, group="north"),
            sc.Helix(max_offset=21, group="east"),
        ]
        group_north = sc.HelixGroup(position=sc.Position3D(x=0, y=-200, z=0), grid=sc.honeycomb)
        groups = {self.n: group_north}
        with self.assertRaises(sc.IllegalDesignError):
            sc.Design(helices=helices, groups=groups, strands=[])

    def _asserts_for_fixture(self, design: sc.Design):
        n = self.n
        e = self.e
        s = self.s
        w = self.w
        groups = design.groups
        if groups is None:
            return  # this makes MyPy shut up about how groups might be None

        self.assertEqual(4, len(groups))

        self.assertSequenceEqual([0, 1, 2, 3, 4, 5], groups[n].helices_view_order)
        self.assertSequenceEqual([7, 6], groups[s].helices_view_order)
        self.assertSequenceEqual([8, 9], groups[w].helices_view_order)
        self.assertSequenceEqual([13, 15], groups[e].helices_view_order)

        self.assertEqual(sc.Grid.honeycomb, groups[n].grid)
        self.assertEqual(sc.Grid.square, groups[e].grid)
        self.assertEqual(sc.Grid.square, groups[s].grid)
        self.assertEqual(sc.Grid.none, groups[w].grid)

        self.assertAlmostEqual(0, groups[n].pitch)
        self.assertAlmostEqual(45, groups[e].pitch)
        self.assertAlmostEqual(0, groups[s].pitch)
        self.assertAlmostEqual(0, groups[w].pitch)

        # test auto-assignment of grid_positions based on helices view order
        self.assertSequenceEqual([0, 1], design.helices[6].grid_position)
        self.assertSequenceEqual([0, 0], design.helices[7].grid_position)
        self.assertSequenceEqual([0, 0], design.helices[13].grid_position)
        self.assertSequenceEqual([0, 1], design.helices[15].grid_position)

    def test_JSON_bad_uses_groups_and_top_level_grid(self) -> None:
        json_str = '''
{
  "grid": "none",
  "groups": {
    "north": {
      "position": {"x": 0, "y": -200, "z": 0},
      "grid": "honeycomb"
    },
    "east": {
      "position": {"x": 0, "y": 0, "z": 100},
      "pitch": 45,
      "grid": "square"
    }
  },
  "helices": [
    {"group": "north", "max_offset": 20, "grid_position": [1, 1]},
    {"group": "north", "max_offset": 21, "grid_position": [0, 1]},
    {"group": "east", "max_offset": 22, "grid_position": [0, 13]},
    {"group": "east", "max_offset": 23, "grid_position": [0, 15]}
  ],
  "strands": [
    {
      "color": "#f74308",
      "domains": [
        {"helix": 0, "forward": true, "start": 0, "end": 8},
        {"helix": 1, "forward": false, "start": 0, "end": 8}
      ]
    },
    {
      "color": "#57bb00",
      "domains": [
        {"helix": 2, "forward": true, "start": 0, "end": 8},
        {"helix": 3, "forward": false, "start": 0, "end": 8}
      ]
    }
  ]
}
        '''
        with self.assertRaises(sc.IllegalDesignError):
            sc.Design.from_scadnano_json_str(json_str)

    def test_JSON_bad_uses_groups_and_top_level_helices_view_order(self) -> None:
        json_str = '''
{
  "helices_view_order": [3, 2, 1, 0],
  "groups": {
    "north": {
      "position": {"x": 0, "y": -200, "z": 0},
      "grid": "honeycomb"
    },
    "east": {
      "position": {"x": 0, "y": 0, "z": 100},
      "pitch": 45,
      "grid": "square"
    }
  },
  "helices": [
    {"group": "north", "max_offset": 20, "grid_position": [1, 1]},
    {"group": "north", "max_offset": 21, "grid_position": [0, 1]},
    {"group": "east", "max_offset": 22, "grid_position": [0, 13]},
    {"group": "east", "max_offset": 23, "grid_position": [0, 15]}
  ],
  "strands": [
    {
      "color": "#f74308",
      "domains": [
        {"helix": 0, "forward": true, "start": 0, "end": 8},
        {"helix": 1, "forward": false, "start": 0, "end": 8}
      ]
    },
    {
      "color": "#57bb00",
      "domains": [
        {"helix": 2, "forward": true, "start": 0, "end": 8},
        {"helix": 3, "forward": false, "start": 0, "end": 8}
      ]
    }
  ]
}
        '''
        with self.assertRaises(sc.IllegalDesignError):
            sc.Design.from_scadnano_json_str(json_str)

    def test_JSON_bad_no_groups_but_helices_reference_groups(self) -> None:
        json_str = '''
{
  "grid": "square",
  "helices": [
    {"group": "north", "max_offset": 20, "grid_position": [0, 0]},
    {"group": "north", "max_offset": 21, "grid_position": [0, 1]},
    {"group": "east", "max_offset": 22, "grid_position": [0, 2]},
    {"group": "east", "max_offset": 23, "grid_position": [0, 3]}
  ],
  "strands": [
    {
      "color": "#f74308",
      "domains": [
        {"helix": 0, "forward": true, "start": 0, "end": 8},
        {"helix": 1, "forward": false, "start": 0, "end": 8}
      ]
    },
    {
      "color": "#57bb00",
      "domains": [
        {"helix": 2, "forward": true, "start": 0, "end": 8},
        {"helix": 3, "forward": false, "start": 0, "end": 8}
      ]
    }
  ]
}

'''
        with self.assertRaises(sc.IllegalDesignError):
            sc.Design.from_scadnano_json_str(json_str)


class TestNames(unittest.TestCase):

    def test_strand_domain_names_json(self) -> None:
        helices = [sc.Helix(max_offset=100), sc.Helix(max_offset=100), sc.Helix(max_offset=100)]
        strand0 = sc.Strand(name='strand0', domains=[
            sc.Domain(0, True, 0, 8, name='domain_forward0'),
            sc.Domain(0, True, 8, 16, name='domain_forward1'),
        ])
        strand1 = sc.Strand(name='strand1', domains=[
            sc.Domain(0, False, 0, 8, name='domain_reverse0'),
            sc.Loopout(3, name='loopout'),
            sc.Domain(2, True, 0, 8, name='domain_forward2'),
        ])
        strand2 = sc.Strand(domains=[
            sc.Domain(1, True, 0, 8),
            sc.Domain(1, True, 8, 16, name='domain_forward0'),
        ])
        strands = [strand0, strand1, strand2]
        design = sc.Design(helices=helices, strands=strands, grid=sc.square)

        json_map = design.to_json_serializable(suppress_indent=False)

        self.assertIn(sc.strand_name_key, json_map[sc.strands_key][0].keys())
        self.assertIn(sc.strand_name_key, json_map[sc.strands_key][1].keys())
        self.assertNotIn(sc.strand_name_key, json_map[sc.strands_key][2].keys())
        self.assertEqual('strand0', json_map[sc.strands_key][0][sc.strand_name_key])
        self.assertEqual('strand1', json_map[sc.strands_key][1][sc.strand_name_key])

        self.assertIn(sc.domain_name_key, json_map[sc.strands_key][0]['domains'][0])
        self.assertIn(sc.domain_name_key, json_map[sc.strands_key][0]['domains'][1])
        self.assertEqual('domain_forward0', json_map[sc.strands_key][0]['domains'][0][sc.domain_name_key])
        self.assertEqual('domain_forward1', json_map[sc.strands_key][0]['domains'][1][sc.domain_name_key])

        self.assertIn(sc.domain_name_key, json_map[sc.strands_key][1]['domains'][0])
        self.assertIn(sc.domain_name_key, json_map[sc.strands_key][1]['domains'][1])
        self.assertIn(sc.domain_name_key, json_map[sc.strands_key][1]['domains'][2])
        self.assertEqual('domain_reverse0', json_map[sc.strands_key][1]['domains'][0][sc.domain_name_key])
        self.assertEqual('loopout', json_map[sc.strands_key][1]['domains'][1][sc.domain_name_key])
        self.assertEqual('domain_forward2', json_map[sc.strands_key][1]['domains'][2][sc.domain_name_key])

        self.assertIn(sc.domain_name_key, json_map[sc.strands_key][2]['domains'][1])
        self.assertEqual('domain_forward0', json_map[sc.strands_key][2]['domains'][1][sc.domain_name_key])

        design_from_json = sc.Design.from_scadnano_json_map(json_map)

        self.assertEqual(3, len(design_from_json.strands))
        self.assertEqual('strand0', design_from_json.strands[0].name)
        self.assertEqual('strand1', design_from_json.strands[1].name)
        self.assertEqual(None, design_from_json.strands[2].name)

        self.assertEqual('domain_forward0', design_from_json.strands[0].domains[0].name)
        self.assertEqual('domain_forward1', design_from_json.strands[0].domains[1].name)

        self.assertEqual('domain_reverse0', design_from_json.strands[1].domains[0].name)
        self.assertEqual('loopout', design_from_json.strands[1].domains[1].name)
        self.assertEqual('domain_forward2', design_from_json.strands[1].domains[2].name)

        self.assertEqual(None, design_from_json.strands[2].domains[0].name)
        self.assertEqual('domain_forward0', design_from_json.strands[2].domains[1].name)

    def test_strand_names_can_be_nonunique(self) -> None:
        helices = [sc.Helix(max_offset=100)]
        strand0 = sc.Strand(name='strand0', domains=[
            sc.Domain(0, True, 0, 8, name='domain_forward0'),
            sc.Domain(0, True, 8, 16, name='domain_forward1'),
        ])
        strand1 = sc.Strand(name='strand0', domains=[
            sc.Domain(0, False, 0, 8, name='domain_reverse0'),
            sc.Domain(0, False, 8, 16, name='domain_reverse1'),
        ])
        strands = [strand0, strand1]
        sc.Design(helices=helices, strands=strands, grid=sc.square)


class TestJSON(unittest.TestCase):

    def test_grid_design_level_converted_to_enum_from_string(self) -> None:
        # reproduces an error where the grid was stored as a string instead of the Grid enum type
        json_str = '''
{
  "version": "0.14.0",
  "grid": "square",
  "helices": [
    {"grid_position": [0, 0]},
    {"grid_position": [0, 1]}
  ],
  "strands": []
}
    '''
        design = sc.Design.from_scadnano_json_str(json_str)
        grid = design.grid
        self.assertTrue(type(grid) is sc.Grid)
        self.assertEqual(sc.Grid.square, grid)

    def test_grid_helix_group_level_converted_to_enum_from_string(self) -> None:
        # reproduces an error where the grid was stored as a string instead of the Grid enum type
        json_str = '''
{
  "version": "0.15.0",
  "groups": {
    "north": {
      "position": {"x": 0, "y": -10, "z": 0},
      "grid": "honeycomb"
    },
    "east": {
      "position": {"x": 0, "y": 0, "z": 10},
      "grid": "square"
    },
    "south": {
      "position": {"x": 0, "y": 10, "z": 0},
      "grid": "hex"
    },
    "west": {
      "position": {"x": 0, "y": 0, "z": -10},
      "grid": "none"
    }
  },
  "helices": [
    {"group": "north", "max_offset": 20, "grid_position": [0, 0], "idx": 0},
    {"group": "north", "max_offset": 21, "grid_position": [1, 0], "idx": 1},
    {"group": "north", "max_offset": 19, "grid_position": [1, 1], "idx": 2},
    {"group": "north", "max_offset": 18, "grid_position": [0, 1], "idx": 3},
    {"group": "north", "max_offset": 17, "grid_position": [-1, 1], "idx": 4},
    {"group": "north", "max_offset": 16, "grid_position": [-1, 0], "idx": 5},
    {"group": "south", "max_offset": 24, "grid_position": [0, 1], "idx": 6},
    {"group": "south", "max_offset": 25, "grid_position": [0, 0], "idx": 7},
    {"group": "west", "max_offset": 26, "position": {"x": 0, "y": 0, "z": 0}, "idx": 8},
    {"group": "west", "max_offset": 27, "position": {"x": 0, "y": 3, "z": 0}, "idx": 9},
    {"group": "east", "max_offset": 22, "grid_position": [0, 0], "idx": 13},
    {"group": "east", "max_offset": 23, "grid_position": [0, 1], "idx": 15}
  ],
  "strands": []
}
    '''
        design = sc.Design.from_scadnano_json_str(json_str)

        grid_n = design.groups['north'].grid
        grid_e = design.groups['east'].grid
        grid_s = design.groups['south'].grid
        grid_w = design.groups['west'].grid
        for grid in [grid_n, grid_e, grid_s, grid_w]:
            self.assertTrue(type(grid) is sc.Grid)

        self.assertEqual(sc.Grid.honeycomb, grid_n)
        self.assertEqual(sc.Grid.square, grid_e)
        self.assertEqual(sc.Grid.hex, grid_s)
        self.assertEqual(sc.Grid.none, grid_w)

    def test_legacy_idt_name_import__no_strand_name(self) -> None:
        # tests proper importing of old format when name was a subfield of idt;
        # ensures if that exists and no Strand.name field exists, the idt.name is used as Strand.name
        json_str = '''
    {
      "version": "0.14.0",
      "grid": "square",
      "helices": [
        {"grid_position": [0, 0]},
        {"grid_position": [0, 1]}
      ],
      "strands": [
        {
          "color": "#f74308",
          "sequence": "TATTATAGTCTTACCCTGAC",
          "idt": {"name": "staple1", "scale": "100nm", "purification": "HPLC", "plate": "plate1", "well": "A1"},
          "domains": [
            {"helix": 0, "forward": true, "start": 0, "end": 10},
            {"helix": 1, "forward": false, "start": 0, "end": 10}
          ]
        }
      ]
    }
    '''
        design = sc.Design.from_scadnano_json_str(json_str)
        self.assertEqual(1, len(design.strands))
        strand = design.strands[0]
        self.assertEqual('staple1', strand.name)
        self.assertEqual('100nm', strand.idt.scale)
        self.assertEqual('HPLC', strand.idt.purification)
        self.assertEqual('plate1', strand.idt.plate)
        self.assertEqual('A1', strand.idt.well)

    def test_legacy_idt_name_import__strand_name_exists(self) -> None:
        # tests proper importing of old format when name was a subfield of idt;
        # ensures if both exist, we use Strand.name
        json_str = '''
    {
      "version": "0.14.0",
      "grid": "square",
      "helices": [
        {"grid_position": [0, 0]},
        {"grid_position": [0, 1]}
      ],
      "strands": [
        {
          "name": "staple1 strand level",
          "color": "#f74308",
          "sequence": "TATTATAGTCTTACCCTGAC",
          "idt": {"name": "staple1", "scale": "100nm", "purification": "HPLC", "plate": "plate1", "well": "A1"},
          "domains": [
            {"helix": 0, "forward": true, "start": 0, "end": 10},
            {"helix": 1, "forward": false, "start": 0, "end": 10}
          ]
        }
      ]
    }
    '''
        design = sc.Design.from_scadnano_json_str(json_str)
        self.assertEqual(1, len(design.strands))
        strand = design.strands[0]
        self.assertEqual('staple1 strand level', strand.name)

    def test_Helix_major_tick_start_default_min_offset(self) -> None:
        helices = [
            sc.Helix(min_offset=10, max_offset=100),
            sc.Helix(max_offset=100),
            sc.Helix(major_tick_start=15),
        ]
        design = sc.Design(helices=helices, strands=[], grid=sc.square)
        self.assertEqual(10, design.helices[0].major_tick_start)
        self.assertEqual(0, design.helices[1].major_tick_start)
        self.assertEqual(15, design.helices[2].major_tick_start)

        design_json_map = design.to_json_serializable(suppress_indent=False)
        self.assertNotIn(sc.major_tick_start_key, design_json_map['helices'][0])
        self.assertNotIn(sc.major_tick_start_key, design_json_map['helices'][1])
        self.assertIn(sc.major_tick_start_key, design_json_map['helices'][2])
        self.assertEqual(15, design_json_map['helices'][2][sc.major_tick_start_key])

        # this isn't related to major_tick_start, but it failed for some reason so let's check it
        self.assertIn(sc.grid_position_key, design_json_map['helices'][0])
        self.assertIn(sc.grid_position_key, design_json_map['helices'][1])
        self.assertIn(sc.grid_position_key, design_json_map['helices'][2])
        self.assertSequenceEqual([0, 0], design_json_map['helices'][0][sc.grid_position_key])
        self.assertSequenceEqual([0, 1], design_json_map['helices'][1][sc.grid_position_key])
        self.assertSequenceEqual([0, 2], design_json_map['helices'][2][sc.grid_position_key])

        design_json_str = json.dumps(design_json_map)
        design = sc.Design.from_scadnano_json_str(design_json_str)
        self.assertEqual(10, design.helices[0].major_tick_start)
        self.assertEqual(0, design.helices[1].major_tick_start)
        self.assertEqual(15, design.helices[2].major_tick_start)

    def test_Helix_major_tick_periodic_distances(self) -> None:
        grid = sc.square
        helices = [
            sc.Helix(major_tick_start=10, max_offset=30, major_tick_distance=5),
            sc.Helix(major_tick_start=10, max_offset=30, major_tick_periodic_distances=[2, 3]),
            sc.Helix(major_tick_start=10, max_offset=30, major_ticks=[10, 20, 30]),
            sc.Helix(major_tick_start=10, max_offset=30),
            sc.Helix(max_offset=30),
        ]
        design = sc.Design(helices=helices, strands=[], grid=grid)
        self.assertEqual(10, design.helices[0].major_tick_start)

        self.assertSequenceEqual([10, 15, 20, 25, 30], design.helices[0].calculate_major_ticks(grid))
        self.assertSequenceEqual([10, 12, 15, 17, 20, 22, 25, 27, 30],
                                 design.helices[1].calculate_major_ticks(grid))
        self.assertSequenceEqual([10, 20, 30], design.helices[2].calculate_major_ticks(grid))
        self.assertSequenceEqual([10, 18, 26], design.helices[3].calculate_major_ticks(grid))
        self.assertSequenceEqual([0, 8, 16, 24], design.helices[4].calculate_major_ticks(grid))

        design_json_map = design.to_json_serializable(suppress_indent=False)

        h0 = design_json_map['helices'][0]
        self.assertNotIn(sc.major_ticks_key, h0)
        self.assertNotIn(sc.major_tick_periodic_distances_key, h0)
        self.assertIn(sc.major_tick_distance_key, h0)
        self.assertEqual(5, h0[sc.major_tick_distance_key])

        h1 = design_json_map['helices'][1]
        self.assertNotIn(sc.major_ticks_key, h1)
        self.assertIn(sc.major_tick_periodic_distances_key, h1)
        self.assertNotIn(sc.major_tick_distance_key, h1)
        self.assertSequenceEqual([2, 3], h1[sc.major_tick_periodic_distances_key])

        h2 = design_json_map['helices'][2]
        self.assertIn(sc.major_ticks_key, h2)
        self.assertNotIn(sc.major_tick_distance_key, h2)
        self.assertNotIn(sc.major_tick_periodic_distances_key, h2)
        self.assertSequenceEqual([10, 20, 30], h2[sc.major_ticks_key])

        h3 = design_json_map['helices'][3]
        self.assertNotIn(sc.major_ticks_key, h3)
        self.assertNotIn(sc.major_tick_distance_key, h3)
        self.assertNotIn(sc.major_tick_periodic_distances_key, h3)
        self.assertIn(sc.major_tick_start_key, h3)

        h4 = design_json_map['helices'][4]
        self.assertNotIn(sc.major_ticks_key, h4)
        self.assertNotIn(sc.major_tick_distance_key, h4)
        self.assertNotIn(sc.major_tick_periodic_distances_key, h4)
        self.assertNotIn(sc.major_tick_start_key, h4)

    def test_default_helices_view_order_with_nondefault_helix_idxs_in_default_order(self) -> None:
        helices = [sc.Helix(idx=1, max_offset=100), sc.Helix(idx=3, max_offset=100)]
        design = sc.Design(helices=helices, strands=[])
        self.assertListEqual([1, 3], design.helices_view_order)

        # [1, 3] is default so json should not contain key
        design_json_ser = design.to_json_serializable(suppress_indent=False)
        self.assertFalse(sc.helices_view_order_key in design_json_ser)

    def test_default_helices_view_order_with_nondefault_helix_idxs_in_nondefault_order(self) -> None:
        helices = [sc.Helix(idx=1, max_offset=100), sc.Helix(idx=3, max_offset=100)]
        design = sc.Design(helices=helices, strands=[], helices_view_order=[3, 1])
        self.assertListEqual([3, 1], design.helices_view_order)

        # [1, 3] is default so json should not contain key
        design_json_ser = design.to_json_serializable(suppress_indent=False)
        actual_view_order = design_json_ser[sc.helices_view_order_key]
        self.assertListEqual([3, 1], actual_view_order)

    def test_strand_labels(self) -> None:
        helices = [sc.Helix(max_offset=100), sc.Helix(max_offset=100)]
        strand0_expected = sc.Strand([sc.Domain(0, True, 0, 10)], label={
            'name': 'strand 0',
            'num_domains': 1,
        })
        strand1_expected = sc.Strand([sc.Domain(0, False, 0, 10), sc.Domain(1, True, 0, 10)], label={
            'name': 'strand 1',
            'num_domains': 2,
        })
        strands = [strand0_expected, strand1_expected]
        design = sc.Design(helices=helices, strands=strands, grid=sc.square)
        json_str = design.to_json()
        design_from_json = sc.Design.from_scadnano_json_str(json_str)
        strand0 = design_from_json.strands[0]
        strand1 = design_from_json.strands[1]
        self.assertDictEqual(strand0_expected.label, strand0.label)
        self.assertDictEqual(strand1_expected.label, strand1.label)

    def test_strand_idt(self) -> None:
        helices = [sc.Helix(max_offset=100), sc.Helix(max_offset=100)]
        idt = sc.IDTFields(scale='25nm', purification='HPLC', plate='plate1', well='A2')
        strand0_expected = sc.Strand([sc.Domain(0, True, 0, 10)], name='strand1', idt=idt)
        strands = [strand0_expected]
        design = sc.Design(helices=helices, strands=strands, grid=sc.square)
        json_str = design.to_json()
        design_from_json = sc.Design.from_scadnano_json_str(json_str)
        strand0 = design_from_json.strands[0]
        self.assertEqual(strand0_expected.idt, strand0.idt)
        self.assertEqual('strand1', strand0.idt_export_name())
        self.assertEqual('25nm', strand0.idt.scale)
        self.assertEqual('HPLC', strand0.idt.purification)
        self.assertEqual('plate1', strand0.idt.plate)
        self.assertEqual('A2', strand0.idt.well)

    def test_domain_labels(self) -> None:
        helices = [sc.Helix(max_offset=100), sc.Helix(max_offset=100)]
        dom00_expected = sc.Domain(0, True, 0, 10, label='domain 00')
        dom10_expected = sc.Domain(0, False, 0, 10)
        dom11_expected = sc.Domain(1, True, 0, 10, label='domain 11')
        strand0 = sc.Strand([dom00_expected])
        strand1 = sc.Strand([dom10_expected, dom11_expected])
        strands = [strand0, strand1]
        design = sc.Design(helices=helices, strands=strands, grid=sc.square)
        json_str = design.to_json()
        design_from_json = sc.Design.from_scadnano_json_str(json_str)
        dom00 = design_from_json.strands[0].domains[0]
        dom10 = design_from_json.strands[1].domains[0]
        dom11 = design_from_json.strands[1].domains[1]
        self.assertEqual(dom00_expected.label, dom00.label)
        self.assertIsNone(dom10.label)
        self.assertEqual(dom11_expected.label, dom11.label)

    def test_nondefault_geometry(self) -> None:
        geometry_expected = sc.Geometry(rise_per_base_pair=10.0, helix_radius=4.0, bases_per_turn=11.0,
                                        minor_groove_angle=10.0,
                                        inter_helix_gap=5.0)
        design = sc.Design(helices=[], strands=[], geometry=geometry_expected)
        json_str = design.to_json()
        design_from_json = sc.Design.from_scadnano_json_str(json_str)
        geometry_actual = design_from_json.geometry
        self.assertAlmostEqual(geometry_expected.rise_per_base_pair, geometry_actual.rise_per_base_pair)
        self.assertAlmostEqual(geometry_expected.helix_radius, geometry_actual.helix_radius)
        self.assertAlmostEqual(geometry_expected.bases_per_turn, geometry_actual.bases_per_turn)
        self.assertAlmostEqual(geometry_expected.minor_groove_angle, geometry_actual.minor_groove_angle)
        self.assertAlmostEqual(geometry_expected.inter_helix_gap, geometry_actual.inter_helix_gap)

    def test_nondefault_geometry_some_default(self) -> None:
        geometry_expected = sc.Geometry(rise_per_base_pair=10.0, minor_groove_angle=10.0, inter_helix_gap=5.0)
        design = sc.Design(helices=[], strands=[], geometry=geometry_expected)
        json_str = design.to_json()
        design_from_json = sc.Design.from_scadnano_json_str(json_str)
        geometry_actual = design_from_json.geometry
        self.assertAlmostEqual(geometry_expected.rise_per_base_pair, geometry_actual.rise_per_base_pair)
        self.assertAlmostEqual(geometry_expected.helix_radius, geometry_actual.helix_radius)
        self.assertAlmostEqual(geometry_expected.bases_per_turn, geometry_actual.bases_per_turn)
        self.assertAlmostEqual(geometry_expected.minor_groove_angle, geometry_actual.minor_groove_angle)
        self.assertAlmostEqual(geometry_expected.inter_helix_gap, geometry_actual.inter_helix_gap)

    def test_lack_of_NoIndent_on_helix_if_position_or_major_ticks_present(self) -> None:
        helices = [sc.Helix(position=sc.Position3D(0, 0, 0))]
        strands = []
        design = sc.Design(helices=helices, strands=strands)
        json_map = design.to_json_serializable(suppress_indent=True)
        helix_json = json_map[sc.helices_key][0]
        self.assertTrue(isinstance(helix_json, sc.NoIndent))
        # self.assertTrue(isinstance(helix_json[sc.position_key], sc.NoIndent))

    def test_NoIndent_on_helix_without_position_or_major_ticks_present(self) -> None:
        helices = [sc.Helix()]
        strands = []
        # will only be NoIndent if we don't have a Helix.position field, so make sure one isn't created
        # by choosing a grid
        design = sc.Design(helices=helices, strands=strands, grid=sc.square)
        json_map = design.to_json_serializable(suppress_indent=True)
        helix_json = json_map[sc.helices_key][0]
        self.assertTrue(isinstance(helix_json, sc.NoIndent))

    def test_error_when_grid_missing(self) -> None:
        json_str = """
        { 
          "helices": [{"grid_position": [0,0]}],
          "strands": [ 
            { 
              "color": 26316, 
              "domains": [ {"helix": 0, "forward": true, "start": 0, "end": 32} ]
            } 
          ] 
        }
        """
        with self.assertRaises(sc.IllegalDesignError) as ex:
            sc.Design.from_scadnano_json_str(json_str)
        msg = ex.exception.args[0]
        self.assertTrue('grid' in msg)

    def test_error_when_domain_helix_missing(self) -> None:
        json_str = """
        { 
          "grid": "square",
          "helices": [{"grid_position": [0,0]}],
          "strands": [ 
            { 
              "color": 26316, 
              "domains": [ {"forward": true, "start": 0, "end": 32} ]
            } 
          ] 
        }
        """
        with self.assertRaises(sc.IllegalDesignError) as ex:
            sc.Design.from_scadnano_json_str(json_str)
        msg = ex.exception.args[0]
        self.assertTrue('helix' in msg)

    def test_error_when_domain_forward_and_right_missing(self) -> None:
        json_str = """
        { 
          "grid": "square",
          "helices": [{"grid_position": [0,0]}],
          "strands": [ 
            { 
              "color": 26316, 
              "domains": [ {"helix": 0, "start": 0, "end": 32} ]
            } 
          ] 
        }
        """
        with self.assertRaises(sc.IllegalDesignError) as ex:
            sc.Design.from_scadnano_json_str(json_str)
        msg = ex.exception.args[0]
        self.assertTrue('forward' in msg)
        self.assertTrue('right' in msg)

    def test_error_when_domain_start_missing(self) -> None:
        json_str = """
        { 
          "grid": "square",
          "helices": [{"grid_position": [0,0]}],
          "strands": [ 
            { 
              "color": 26316, 
              "domains": [ {"helix": 0, "forward": true, "end": 32} ]
            } 
          ] 
        }
        """
        with self.assertRaises(sc.IllegalDesignError) as ex:
            sc.Design.from_scadnano_json_str(json_str)
        msg = ex.exception.args[0]
        self.assertTrue('start' in msg)

    def test_error_when_domain_end_missing(self) -> None:
        json_str = """
        { 
          "grid": "square",
          "helices": [{"grid_position": [0,0]}],
          "strands": [ 
            { 
              "color": 26316, 
              "domains": [ {"helix": 0, "forward": true, "start": 0 } ]
            } 
          ] 
        }
        """
        with self.assertRaises(sc.IllegalDesignError) as ex:
            sc.Design.from_scadnano_json_str(json_str)
        msg = ex.exception.args[0]
        self.assertTrue('end' in msg)

    def test_error_when_strands_missing(self) -> None:
        json_str = """
        { 
          "grid": "square",
          "helices": [{"grid_position": [0,0]}]
        }
        """
        with self.assertRaises(sc.IllegalDesignError) as ex:
            sc.Design.from_scadnano_json_str(json_str)
        msg = ex.exception.args[0]
        self.assertTrue('strands' in msg)

    def test_legacy_right_key(self) -> None:
        json_str = """
        { 
          "grid": "square",
          "helices": [{"grid_position": [0,0]}],
          "strands": [ 
            { 
              "color": 26316, 
              "domains": [ {"helix": 0, "right": true, "start": 0, "end": 5 } ]
            } 
          ] 
        }
        """
        d = sc.Design.from_scadnano_json_str(json_str)
        self.assertEqual(True, d.strands[0].domains[0].forward)

    def test_legacy_dna_sequence_key(self) -> None:
        json_str = """
        { 
          "grid": "square",
          "helices": [{"grid_position": [0,0]}],
          "strands": [ 
            { 
              "color": 26316, 
              "dna_sequence": "ACGTA",
              "domains": [ {"helix": 0, "right": true, "start": 0, "end": 5 } ]
            } 
          ] 
        }
        """
        d = sc.Design.from_scadnano_json_str(json_str)
        self.assertEqual("ACGTA", d.strands[0].dna_sequence)

    def test_legacy_substrands_key(self) -> None:
        json_str = """
        { 
          "grid": "square",
          "helices": [{"grid_position": [0,0]}],
          "strands": [ 
            { 
              "color": 26316, 
              "substrands": [ {"helix": 0, "forward": true, "start": 0, "end": 5 } ]
            } 
          ] 
        }
        """
        d = sc.Design.from_scadnano_json_str(json_str)
        self.assertEqual(0, d.strands[0].domains[0].helix)
        self.assertEqual(True, d.strands[0].domains[0].forward)
        self.assertEqual(0, d.strands[0].domains[0].start)
        self.assertEqual(5, d.strands[0].domains[0].end)

    def test_color_specified_with_integer(self) -> None:
        # addresses https://github.com/UC-Davis-molecular-computing/scadnano-python-package/issues/58
        # 0066cc hex is 26316 decimal
        json_str = """
        { 
          "grid": "square",
          "helices": [{"grid_position": [0,0]}],
          "strands": [ 
            { 
              "color": 26316, 
              "domains": [ {"helix": 0, "forward": true, "start": 0, "end": 32} ]
            } 
          ] 
        }
        """
        d = sc.Design.from_scadnano_json_str(json_str)
        expected_color_hex = '#0066cc'
        actual_color_hex = d.strands[0].color.to_json_serializable(False)
        self.assertEqual(expected_color_hex, actual_color_hex)

    def test_single_helix_group_and_helices_specify_pitch_and_yaw(self) -> None:
        json_str = """
        {
          "helices": [
            {
              "group": "north",
              "position": {"x": 1, "y": 2, "z": 3},
              "pitch": 4,
              "roll": 5,
              "yaw": 6
            },
            {
              "group": "north",
              "position": {"x": 3, "y": 2, "z": 3},
              "roll": 15
            }
          ],
          "groups": {
            "north": {
              "position": {"x": 0, "y": -200, "z": 0},
              "pitch": 21,
              "yaw": 13,
              "grid": "none"
            }
          },
          "strands": [ 
            { 
              "color": "#0066cc", 
              "domains": [ {"helix": 0, "forward": true, "start": 0, "end": 32} ]
            } 
          ] 
        }
        """
        d = sc.Design.from_scadnano_json_str(json_str)
        helix0 = d.helices[0]
        helix1 = d.helices[1]
        self.assertEqual(sc.Position3D(1, 2, 3), helix0.position)
        self.assertEqual(5, helix0.roll)
        # Helix 0 should have been moved to a new helix group
        pitch_25_yaw_19_group_name = f'pitch_25_yaw_19'
        pitch_25_yaw_19_group = d.groups[pitch_25_yaw_19_group_name]
        self.assertEqual(25, pitch_25_yaw_19_group.pitch)
        self.assertEqual(19, pitch_25_yaw_19_group.yaw)
        self.assertEqual(pitch_25_yaw_19_group_name, helix0.group)
        self.assertEqual(sc.Position3D(3, 2, 3), helix1.position)
        self.assertEqual(15, helix1.roll)
        self.assertEqual("north", helix1.group)
        self.assertEqual(2, len(d.groups))

    def test_only_individual_helices_specify_pitch_and_yaw(self) -> None:
        json_str = """
        {
          "helices": [
            {
              "group": "north",
              "position": {"x": 1, "y": 2, "z": 3},
              "pitch": 25,
              "yaw": 19,
              "roll": 5
            },
            {
              "group": "north",
              "position": {"x": 3, "y": 2, "z": 3},
              "pitch": 21,
              "yaw": 13,
              "roll": 15
            }
          ],
          "groups": {
            "north": {
              "position": {"x": 0, "y": -200, "z": 0},
              "grid": "none"
            }
          },
          "strands": [ 
            { 
              "color": "#0066cc", 
              "domains": [ {"helix": 0, "forward": true, "start": 0, "end": 32} ]
            } 
          ] 
        }
        """
        d = sc.Design.from_scadnano_json_str(json_str)
        helix0 = d.helices[0]
        helix1 = d.helices[1]

        # Helix 0 should have been moved to a new helix group
        pitch_25_yaw_19_group_name = f'pitch_25.0_yaw_19.0'
        pitch_25_yaw_19_group = d.groups[pitch_25_yaw_19_group_name]
        self.assertEqual(sc.Position3D(1, 2, 3), helix0.position)
        self.assertEqual(25, pitch_25_yaw_19_group.pitch)
        self.assertEqual(19, pitch_25_yaw_19_group.yaw)
        self.assertEqual(5, helix0.roll)
        self.assertEqual(pitch_25_yaw_19_group_name, helix0.group)

        # Helix 1 should have been moved to a new helix group
        pitch_21_yaw_13_group_name = f'pitch_21.0_yaw_13.0'
        pitch_21_yaw_13_group = d.groups[pitch_21_yaw_13_group_name]
        self.assertEqual(sc.Position3D(3, 2, 3), helix1.position)
        self.assertEqual(21, pitch_21_yaw_13_group.pitch)
        self.assertEqual(13, pitch_21_yaw_13_group.yaw)
        self.assertEqual(15, helix1.roll)
        self.assertEqual(pitch_21_yaw_13_group_name, helix1.group)

        self.assertEqual(3, len(d.groups))

    def test_only_helix_groups_specify_pitch_and_yaw(self) -> None:
        json_str = """
        {
          "helices": [
            {
              "group": "north",
              "position": {"x": 1, "y": 2, "z": 3},
              "roll": 5
            },
            {
              "group": "south",
              "position": {"x": 3, "y": 2, "z": 3},
              "roll": 15
            }
          ],
          "groups": {
            "north": {
              "position": {"x": 0, "y": -200, "z": 0},
              "pitch": 21,
              "yaw": 13,
              "grid": "none"
            },
            "south": {
              "position": {"x": 0, "y": -400, "z": 0},
              "pitch": 23,
              "yaw": 98,
              "grid": "none"
            }
          },
          "strands": [
            {
              "color": "#0066cc",
              "domains": [ {"helix": 0, "forward": true, "start": 0, "end": 32} ]
            }
          ]
        }
        """
        d = sc.Design.from_scadnano_json_str(json_str)
        helix0 = d.helices[0]
        helix1 = d.helices[1]

        north_str = 'north'
        south_str = 'south'
        north_group = d.groups[north_str]
        south_group = d.groups[south_str]
        self.assertEqual(2, len(d.groups))

        self.assertEqual(sc.Position3D(1, 2, 3), helix0.position)
        self.assertEqual(5, helix0.roll)
        self.assertEqual(21, north_group.pitch)
        self.assertEqual(13, north_group.yaw)
        self.assertEqual(north_str, helix0.group)

        self.assertEqual(sc.Position3D(3, 2, 3), helix1.position)
        self.assertEqual(15, helix1.roll)
        self.assertEqual(23, south_group.pitch)
        self.assertEqual(98, south_group.yaw)
        self.assertEqual(south_str, helix1.group)

    def test_both_helix_groups_and_helices_do_not_specify_pitch_nor_yaw(self) -> None:
        json_str = """
        {
          "helices": [
            {
              "group": "north",
              "position": {"x": 1, "y": 2, "z": 3},
              "roll": 5
            },
            {
              "group": "south",
              "position": {"x": 3, "y": 2, "z": 3},
              "roll": 15
            }
          ],
          "groups": {
            "north": {
              "position": {"x": 0, "y": -200, "z": 0},
              "grid": "none"
            },
            "south": {
              "position": {"x": 0, "y": -400, "z": 0},
              "grid": "none"
            }
          },
          "strands": [
            {
              "color": "#0066cc",
              "domains": [ {"helix": 0, "forward": true, "start": 0, "end": 32} ]
            }
          ]
        }
        """
        d = sc.Design.from_scadnano_json_str(json_str)
        helix0 = d.helices[0]
        helix1 = d.helices[1]

        north_str = 'north'
        south_str = 'south'
        north_group = d.groups[north_str]
        south_group = d.groups[south_str]
        self.assertEqual(2, len(d.groups))

        self.assertEqual(sc.Position3D(1, 2, 3), helix0.position)
        self.assertEqual(5, helix0.roll)
        self.assertEqual(0, north_group.pitch)
        self.assertEqual(0, north_group.yaw)
        self.assertEqual(north_str, helix0.group)

        self.assertEqual(sc.Position3D(3, 2, 3), helix1.position)
        self.assertEqual(15, helix1.roll)
        self.assertEqual(0, south_group.pitch)
        self.assertEqual(0, south_group.yaw)
        self.assertEqual(south_str, helix1.group)

    def test_multiple_helix_groups_helices_specify_pitch_and_yaw(self) -> None:
        json_str = """
        {
          "helices": [
            {
              "group": "north",
              "position": {"x": 1, "y": 2, "z": 3},
              "pitch": 4,
              "roll": 5,
              "yaw": 6            },
            {
              "group": "south",
              "position": {"x": 3, "y": 2, "z": 3},
              "roll": 15
            }
          ],
          "groups": {
            "north": {
              "position": {"x": 0, "y": -200, "z": 0},
              "pitch": 21,
              "yaw": 13,
              "grid": "none"
            },
            "south": {
              "position": {"x": 0, "y": -400, "z": 0},
              "pitch": 23,
              "yaw": 98,
              "grid": "none"
            }
          },
          "strands": [
            {
              "color": "#0066cc",
              "domains": [ {"helix": 0, "forward": true, "start": 0, "end": 32} ]
            }
          ]
        }
        """
        d = sc.Design.from_scadnano_json_str(json_str)
        helix0 = d.helices[0]
        helix1 = d.helices[1]

        north_str = 'north'
        south_str = 'south'
        north_group = d.groups[north_str]
        south_group = d.groups[south_str]
        self.assertEqual(2, len(d.groups))

        self.assertEqual(sc.Position3D(1, 2, 3), helix0.position)
        self.assertEqual(5, helix0.roll)
        self.assertEqual(25, north_group.pitch)
        self.assertEqual(19, north_group.yaw)
        self.assertEqual(north_str, helix0.group)

        self.assertEqual(sc.Position3D(3, 2, 3), helix1.position)
        self.assertEqual(15, helix1.roll)
        self.assertEqual(23, south_group.pitch)
        self.assertEqual(98, south_group.yaw)
        self.assertEqual(south_str, helix1.group)

    def test_multiple_helix_groups_helices_specify_pitch_and_yaw_invalid(self) -> None:
        json_str = """
        {
          "helices": [
            {
              "group": "north",
              "position": {"x": 1, "y": 2, "z": 3},
              "pitch": 4,
              "roll": 5,
              "yaw": 6
            },
            {
              "group": "north",
              "position": {"x": 3, "y": 2, "z": 3},
              "pitch": 10,
              "roll": 5,
              "yaw": 84
            }
          ],
          "groups": {
            "north": {
              "position": {"x": 0, "y": -200, "z": 0},
              "grid": "none"
            },
            "east": {
              "position": {"x": 0, "y": 0, "z": 100},
              "pitch": 45,
              "grid": "square"
            }
          },
          "strands": [ 
            { 
              "color": "#0066cc", 
              "domains": [ {"helix": 0, "forward": true, "start": 0, "end": 32} ]
            } 
          ] 
        }
        """
        # Should fail because multiple helices in same helix group are non-parallel
        with self.assertRaises(sc.IllegalDesignError):
            sc.Design.from_scadnano_json_str(json_str)

    def test_position_specified_with_origin_keyword(self) -> None:
        # addresses https://github.com/UC-Davis-molecular-computing/scadnano-python-package/issues/59
        json_str = """
        { 
          "grid": "none",
          "helices": [{
            "origin": {"x": 1, "y": 2, "z": 3}, 
            "pitch": 4, 
            "roll": 5, 
            "yaw": 6 
          }],
          "strands": [ 
            { 
              "color": "#0066cc", 
              "domains": [ {"helix": 0, "forward": true, "start": 0, "end": 32} ]
            } 
          ] 
        }
        """
        d = sc.Design.from_scadnano_json_str(json_str)
        expected_position = sc.Position3D(1, 2, 3)
        expected_pitch = 4
        expected_roll = 5
        expected_yaw = 6
        actual_position = d.helices[0].position
        expected_group_name = f'pitch_{expected_pitch}_yaw_{expected_yaw}'
        expected_group = d.groups[expected_group_name]
        actual_pitch = expected_group.pitch
        actual_roll = d.helices[0].roll
        actual_yaw = expected_group.yaw
        self.assertEqual(expected_position, actual_position)
        self.assertEqual(expected_pitch, actual_pitch)
        self.assertEqual(expected_roll, actual_roll)
        self.assertEqual(expected_yaw, actual_yaw)
        self.assertEqual(expected_group_name, d.helices[0].group)

    def test_json_tristan_example_issue_32(self) -> None:
        json_str = """
        { 
          "version": "0.3.0", 
          "grid": "square",
          "helices": [ 
            {"grid_position": [0, 0]}, 
            {"max_offset": 32, "grid_position": [0, 1]} 
          ], 
          "strands": [ 
            { 
              "color": "#0066cc", 
              "domains": [ {"helix": 0, "forward": true, "start": 0, "end": 32} ], 
              "is_scaffold": true 
            } 
          ] 
        }
        """
        sc.Design.from_scadnano_json_str(json_str)

    def test_to_json__hairpin(self) -> None:
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
        ss_f = sc.Domain(helix=0, forward=True, start=0, end=5)
        loop = sc.Loopout(length=5)
        ss_r = sc.Domain(helix=0, forward=False, start=0, end=5)
        strand_forward = sc.Strand([ss_f, loop, ss_r])
        design = sc.Design(strands=[strand_forward], grid=sc.square)
        design.assign_dna(strand_forward, 'AAACC TGCAC')
        design.to_json()
        # should be no error getting here

    def test_to_json__roll(self) -> None:
        helix = sc.Helix(roll=90)
        ss_f = sc.Domain(helix=0, forward=True, start=0, end=5)
        ss_r = sc.Domain(helix=0, forward=False, start=0, end=5)
        strand_f = sc.Strand([ss_f])
        strand_r = sc.Strand([ss_r])
        design = sc.Design(helices=[helix], strands=[strand_f, strand_r], grid=sc.square)
        design.to_json()
        # should be no error getting here


class TestIllegalStructuresPrevented(unittest.TestCase):

    # def test_to_json__error_if_DNAOrigamiDesign_no_scaffold(self) -> None:
    #     # we are allowed to delay assigning a scaffold to a DNAOrigamiDesign,
    #     # but to_json should fail
    #     st_l = sc.Strand([
    #         sc.Substrand(1, True, 0, 8, insertions=[(4, 2)]),
    #         sc.Substrand(0, False, 0, 8, deletions=[3]),
    #     ])
    #     st_r = sc.Strand([
    #         sc.Substrand(0, False, 8, 16),
    #         sc.Substrand(1, True, 8, 16),
    #     ])
    #     scaf = sc.Strand([
    #         sc.Substrand(1, False, 0, 8, insertions=[(4, 2)]),
    #         sc.Substrand(0, True, 0, 16, deletions=[3]),
    #         sc.Loopout(3),
    #         sc.Substrand(1, False, 8, 16, deletions=[]),
    #     ])
    #     design_pre_json = sc.DNAOrigamiDesign(strands=[st_l, st_r, scaf], grid=sc.square)
    #
    #     with self.assertRaises(sc.IllegalDesignError):
    #         design_pre_json.to_json()

    def test_strands_not_specified_in_Design_constructor(self) -> None:
        design = sc.Design(helices=[])
        self.assertEqual(0, len(design.helices))
        self.assertEqual(0, len(design.strands))

    def test_helices_not_specified_in_Design_constructor(self) -> None:
        design = sc.Design(strands=[])
        self.assertEqual(0, len(design.helices))
        self.assertEqual(0, len(design.strands))

    def test_strands_and_helices_not_specified_in_Design_constructor(self) -> None:
        design = sc.Design()
        self.assertEqual(0, len(design.helices))
        self.assertEqual(0, len(design.strands))

    def test_consecutive_domains_loopout(self) -> None:
        helices = [sc.Helix(max_offset=10)]
        ss1 = sc.Domain(0, True, 0, 3)
        ss2 = sc.Loopout(4)
        ss3 = sc.Loopout(4)
        with self.assertRaises(sc.IllegalDesignError):
            sc.Strand([ss1, ss2, ss3])

        # XXX: we used to allow Strands to violate the loopout rules and caught it only at the design level
        # now the Strand constructor checks, so that means we can't set up a bad Strand for the Design check
        # strand = sc.Strand([ss1, ss2])
        # strand.domains.append(ss3)
        # with self.assertRaises(sc.IllegalDesignError):
        #     sc.Design(helices=helices, strands=[strand], grid=sc.square)

    def test_singleton_loopout(self) -> None:
        helices = [sc.Helix(max_offset=10)]
        loopout = sc.Loopout(4)
        with self.assertRaises(sc.StrandError):
            sc.Strand([loopout])

        # XXX: we used to allow Strands to violate the loopout rules and caught it only at the design level
        # now the Strand constructor checks, so that means we can't set up a bad Strand for the Design check
        # strand = sc.Strand([])
        # strand.domains.append(loopout)
        # with self.assertRaises(sc.StrandError):
        #     sc.Design(helices=helices, strands=[strand], grid=sc.square)

    def test_strand_offset_beyond_maxbases(self) -> None:
        helices = [sc.Helix(max_offset=10)]
        ss1 = sc.Domain(0, True, 0, 20)
        strands = [sc.Strand([ss1])]
        with self.assertRaises(sc.StrandError):
            sc.Design(helices=helices, strands=strands)

    def test_assign_dna__conflicting_sequences_directly_assigned(self) -> None:
        ss_right = sc.Domain(0, True, 0, 5)
        ss_left = sc.Domain(0, False, 0, 5)
        strand_right = sc.Strand([ss_right])
        strand_left = sc.Strand([ss_left])
        design = sc.Design(strands=[strand_left, strand_right])
        design.assign_dna(strand_right, 'ACGTT')
        with self.assertRaises(sc.IllegalDesignError):
            design.assign_dna(strand_right, 'TTTTT')

    def test_assign_dna__conflicting_sequences_indirectly_assigned(self) -> None:
        ss_right = sc.Domain(0, True, 0, 5)
        ss_left = sc.Domain(0, False, 0, 5)
        strand_right = sc.Strand([ss_right])
        strand_left = sc.Strand([ss_left])
        design = sc.Design(strands=[strand_left, strand_right])
        design.assign_dna(strand_right, 'ACGTT')
        with self.assertRaises(sc.IllegalDesignError):
            design.assign_dna(strand_left, 'GGGGG')

    def test_overlapping_caught_in_strange_counterexample(self) -> None:
        # found this counterexample as a simplified version of something caught in practice
        s1_left_ss0 = sc.Domain(0, False, 0, 5)
        s1_ss1 = sc.Domain(0, True, 0, 15)
        s1_right_ss0 = sc.Domain(0, False, 5, 15)
        s1 = sc.Strand([s1_left_ss0, s1_ss1, s1_right_ss0])

        s2_ss1 = sc.Domain(0, True, 10, 20)
        s2_ss0 = sc.Domain(0, False, 10, 20)
        s2 = sc.Strand([s2_ss1, s2_ss0])

        strands = [s1, s2]

        with self.assertRaises(sc.IllegalDesignError):
            sc.Design(strands=strands, grid=sc.square)

    def test_major_tick_outside_range(self) -> None:
        with self.assertRaises(sc.IllegalDesignError):
            sc.Helix(max_offset=9, major_ticks=[2, 5, 10])

    def test_major_tick_just_inside_range(self) -> None:
        sc.Helix(max_offset=9, major_ticks=[0, 5, 9])

    def test_two_illegally_overlapping_strands(self) -> None:
        ss_bot = sc.Domain(helix=0, forward=False, start=0, end=9)
        ss_top = sc.Domain(helix=0, forward=False, start=0, end=9)
        strand_bot = sc.Strand(domains=[ss_bot])
        strand_top = sc.Strand(domains=[ss_top])
        strands = [strand_bot, strand_top]
        with self.assertRaises(sc.IllegalDesignError):
            sc.Design(grid=sc.square, strands=strands)

    def test_two_nonconsecutive_illegally_overlapping_strands(self) -> None:
        ss_top1 = sc.Domain(helix=0, forward=False, start=0, end=5)
        ss_bot = sc.Domain(helix=0, forward=True, start=2, end=9)
        ss_top2 = sc.Domain(helix=0, forward=False, start=4, end=8)
        strand_bot = sc.Strand(domains=[ss_bot])
        strand_top1 = sc.Strand(domains=[ss_top1])
        strand_top2 = sc.Strand(domains=[ss_top2])
        strands = [strand_bot, strand_top1, strand_top2]
        with self.assertRaises(sc.IllegalDesignError):
            sc.Design(grid=sc.square, strands=strands)

    def test_four_legally_leapfrogging_strands(self) -> None:
        ss_top1 = sc.Domain(helix=0, forward=False, start=0, end=20)
        ss_bot1 = sc.Domain(helix=0, forward=True, start=10, end=30)
        ss_top2 = sc.Domain(helix=0, forward=False, start=20, end=40)
        ss_bot2 = sc.Domain(helix=0, forward=True, start=30, end=50)
        strand_bot1 = sc.Strand(domains=[ss_bot1])
        strand_bot2 = sc.Strand(domains=[ss_bot2])
        strand_top1 = sc.Strand(domains=[ss_top1])
        strand_top2 = sc.Strand(domains=[ss_top2])
        strands = [strand_bot1, strand_bot2, strand_top1, strand_top2]
        sc.Design(grid=sc.square, strands=strands)

    def test_strand_references_nonexistent_helix(self) -> None:
        h1 = sc.Helix(max_offset=9)
        h2 = sc.Helix(max_offset=9)
        ss_bot = sc.Domain(helix=2, forward=False, start=0, end=9)
        ss_top = sc.Domain(helix=3, forward=False, start=0, end=9)
        strand_bot = sc.Strand(domains=[ss_bot])
        strand_top = sc.Strand(domains=[ss_top])
        strands = [strand_bot, strand_top]
        with self.assertRaises(sc.IllegalDesignError):
            sc.Design(grid=sc.square, helices=[h1, h2], strands=strands)


class TestInsertRemoveDomains(unittest.TestCase):

    def setUp(self) -> None:
        helices = [sc.Helix(max_offset=100) for _ in range(4)]
        self.design = sc.Design(helices=helices, strands=[])
        self.design.strand(0, 0).to(3).cross(1).to(0).cross(2).to(3).with_sequence('ACA TCT GTG')
        self.strand = self.design.strands[0]

    def test_3_helix_before_design(self) -> None:
        expected_strand_before = sc.Strand([
            sc.Domain(0, True, 0, 3),
            sc.Domain(1, False, 0, 3),
            sc.Domain(2, True, 0, 3),
        ], dna_sequence='ACA TCT GTG'.replace(' ', ''))
        self.assertEqual(expected_strand_before, self.strand)

    def test_insert_domain_with_sequence(self) -> None:
        helices = [sc.Helix(max_offset=100) for _ in range(4)]
        design = sc.Design(helices=helices, strands=[])
        design.strand(0, 0).to(3).cross(1).to(0).cross(3).to(3).with_sequence('ACA TCT GTG')
        strand = design.strands[0]

        expected_strand_before = sc.Strand([
            sc.Domain(0, True, 0, 3),
            sc.Domain(1, False, 0, 3),
            sc.Domain(3, True, 0, 3),
        ])  # , dna_sequence='ACA TCT GTG'.replace(' ', ''))
        self.assertEqual(expected_strand_before, design.strands[0])

        domain = sc.Domain(2, True, 0, 3)
        design.insert_domain(strand, 2, domain)
        expected_strand = sc.Strand([
            sc.Domain(0, True, 0, 3),
            sc.Domain(1, False, 0, 3),
            sc.Domain(2, True, 0, 3),
            sc.Domain(3, True, 0, 3),
        ], dna_sequence='ACA TCT ??? GTG'.replace(' ', ''))
        self.assertEqual(expected_strand, design.strands[0])

    def test_append_domain_with_sequence(self) -> None:
        domain = sc.Domain(3, False, 0, 3)
        self.design.append_domain(self.strand, domain)
        expected_strand = sc.Strand([
            sc.Domain(0, True, 0, 3),
            sc.Domain(1, False, 0, 3),
            sc.Domain(2, True, 0, 3),
            sc.Domain(3, False, 0, 3),
        ], dna_sequence='ACA TCT GTG ???'.replace(' ', ''))
        self.assertEqual(expected_strand, self.strand)

    def test_remove_first_domain_with_sequence(self) -> None:
        self.design.remove_domain(self.strand, self.strand.domains[0])
        expected_strand = sc.Strand([
            sc.Domain(1, False, 0, 3),
            sc.Domain(2, True, 0, 3),
        ], dna_sequence='    TCT GTG'.replace(' ', ''))
        self.assertEqual(expected_strand, self.strand)

    def test_remove_middle_domain_with_sequence(self) -> None:
        self.design.remove_domain(self.strand, self.strand.domains[1])
        expected_strand = sc.Strand([
            sc.Domain(0, True, 0, 3),
            sc.Domain(2, True, 0, 3),
        ], dna_sequence='ACA     GTG'.replace(' ', ''))
        self.assertEqual(expected_strand, self.strand)

    def test_remove_last_domain_with_sequence(self) -> None:
        self.design.remove_domain(self.strand, self.strand.domains[2])
        expected_strand = sc.Strand([
            sc.Domain(0, True, 0, 3),
            sc.Domain(1, False, 0, 3),
        ], dna_sequence='ACA TCT GTG'.replace(' ', ''))
        self.assertEqual(expected_strand, self.strand)


class TestLabels(unittest.TestCase):

    def setUp(self) -> None:
        helices = [sc.Helix(max_offset=100) for _ in range(10)]
        self.design = sc.Design(helices=helices, strands=[], grid=sc.square)

    def test_with_label__str(self) -> None:
        label = 'abc'
        self.design.strand(0, 0).to(5).cross(1).to(0).with_label(label)
        actual_strand = self.design.strands[0]
        expected_strand = sc.Strand(domains=[
            sc.Domain(0, True, 0, 5),
            sc.Domain(0, False, 0, 5),
        ], label=label)

        self.assertEqual(expected_strand.label, actual_strand.label)

    def test_with_label__dict(self) -> None:
        label = {'name': 'abc', 'type': 3}
        self.design.strand(0, 0).to(5).cross(1).to(0).with_label(label)
        actual_strand = self.design.strands[0]
        expected_strand = sc.Strand(domains=[
            sc.Domain(0, True, 0, 5),
            sc.Domain(0, False, 0, 5),
        ], label=label)

        self.assertDictEqual(expected_strand.label, actual_strand.label)

    def test_with_domain_label(self) -> None:
        label0: Union[str, Dict[str, Any]] = 'abc'
        label1: Union[str, Dict[str, Any]] = {'name': 'abc', 'type': 3}
        self.design.strand(0, 0).to(5).with_domain_label(label0).cross(1).to(0).with_domain_label(label1)
        actual_strand = self.design.strands[0]
        expected_strand = sc.Strand(domains=[
            sc.Domain(0, True, 0, 5, label=label0),
            sc.Domain(0, False, 0, 5, label=label1),
        ])

        self.assertEqual(expected_strand.domains[0].label, actual_strand.domains[0].label)
        self.assertDictEqual(expected_strand.domains[1].label, actual_strand.domains[1].label)

    def test_with_domain_label__and__with_label(self) -> None:
        strand_label = 'xyz'
        label0: Union[str, Dict[str, Any]] = 'abc'
        label1: Union[str, Dict[str, Any]] = {'name': 'abc', 'type': 3}
        self.design.strand(0, 0).to(5).with_domain_label(label0).cross(1).to(0).with_domain_label(label1) \
            .with_label(strand_label)
        actual_strand = self.design.strands[0]
        expected_strand = sc.Strand(domains=[
            sc.Domain(0, True, 0, 5, label=label0),
            sc.Domain(0, False, 0, 5, label=label1),
        ], label=strand_label)

        self.assertEqual(expected_strand.label, actual_strand.label)
        self.assertEqual(expected_strand.domains[0].label, actual_strand.domains[0].label)
        self.assertDictEqual(expected_strand.domains[1].label, actual_strand.domains[1].label)


def set_colors_black(*strands) -> None:
    for strand in strands:
        strand.set_color(sc.Color(r=0, g=0, b=0))


class TestCircularStrandsLegalMods(unittest.TestCase):
    '''
    Tests that circular strands cannot have 5' or 3' mods.
    '''

    def setUp(self) -> None:
        helices = [sc.Helix(max_offset=10) for _ in range(2)]
        self.design = sc.Design(helices=helices, strands=[])
        self.design.strand(0, 0).move(10).cross(1).move(-10)
        self.strand = self.design.strands[0]
        r'''
        0  [--------\
                    |
        1  <--------/
        '''

    def test_can_add_internal_mod_to_circular_strand(self) -> None:
        self.strand.set_circular()
        self.assertTrue(self.strand.circular)
        self.strand.set_modification_internal(2, mod.biotin_int)
        self.assertEqual(1, len(self.strand.modifications_int))

    def test_cannot_make_strand_circular_if_5p_mod(self) -> None:
        self.strand.set_modification_5p(mod.biotin_5p)
        with self.assertRaises(sc.StrandError):
            self.strand.set_circular(True)

    def test_cannot_make_strand_circular_if_3p_mod(self) -> None:
        self.strand.set_modification_3p(mod.biotin_3p)
        with self.assertRaises(sc.StrandError):
            self.strand.set_circular(True)

    def test_add_5p_mod_to_circular_strand(self) -> None:
        self.strand.set_circular(True)
        with self.assertRaises(sc.StrandError):
            self.strand.set_modification_5p(mod.biotin_5p)

    def test_add_3p_mod_to_circular_strand(self) -> None:
        self.strand.set_circular(True)
        with self.assertRaises(sc.StrandError):
            self.strand.set_modification_3p(mod.biotin_3p)


class TestCircularStrandEdits(unittest.TestCase):
    '''
    Tests that circular strand edits (nicking, ligating, adding crossovers/loopouts)
    works properly with circular strands.
    '''

    def setUp(self) -> None:
        helices = [sc.Helix(max_offset=10) for _ in range(3)]
        self.design = sc.Design(helices=helices, strands=[])
        self.design.strand(0, 0).move(10).cross(1).move(-10)
        self.design.strand(0, 15).move(5).cross(1).move(-10).cross(0).move(5)
        self.design.strand(0, 20).move(10)
        self.design.strand(1, 30).move(-10)
        self.design.strand(0, 30).move(10).loopout(1, 5).move(-10).cross(2).move(10).as_circular()
        self.design.strand(0, 40).move(10).cross(1).move(-10).as_circular()
        self.num_strands = len(self.design.strands)
        r'''
           0          10         20         30         40
            strand 0   strand 1   strand 2   strand 4   strand 5
        0  [--------\ /--->[---\ [--------> [--------\ /--------\
                    | |        |                     ) |        |
        1  <--------/ \--------/ <--------] /--------/ \--------/
                                  strand 3  |
        2                                   \-------->cross to 5' end here
        '''

    def test_add_crossover_from_linear_strand_to_itself_makes_it_circular(self) -> None:
        # add crossover to strand 0
        r'''
           0
            strand 0
        0  /--------\
           |        |
        1  \--------/
        '''
        self.assertEqual(2, len(self.design.strands[0].domains))
        self.assertFalse(self.design.strands[0].circular)
        self.design.add_half_crossover(0, 1, 0, True)
        self.assertTrue(self.design.strands[0].circular)
        self.assertEqual(self.num_strands, len(self.design.strands))
        self.assertEqual(2, len(self.design.strands[0].domains))

    def test_add_nick_to_2_domain_circular_strand_makes_it_linear_nick_first_domain(self) -> None:
        r'''
           40
            strand 5
        0  /--->[---\
           |        |
        1  \--------/
        '''
        self.assertTrue(self.design.strands[5].circular)

        # nick strand 5
        self.design.add_nick(0, 45, True)
        self.assertFalse(self.design.strands[5].circular)
        self.assertEqual(self.num_strands, len(self.design.strands))

        strand = self.design.strands[5]
        self.assertEqual(3, len(strand.domains))
        self.assertEqual(3, len(strand.bound_domains()))

        d0, d1, d2 = strand.bound_domains()
        self.assertEqual(0, d0.helix)
        self.assertEqual(True, d0.forward)
        self.assertEqual(45, d0.start)
        self.assertEqual(50, d0.end)

        self.assertEqual(1, d1.helix)
        self.assertEqual(False, d1.forward)
        self.assertEqual(40, d1.start)
        self.assertEqual(50, d1.end)

        self.assertEqual(0, d2.helix)
        self.assertEqual(True, d2.forward)
        self.assertEqual(40, d2.start)
        self.assertEqual(45, d2.end)

        for domain in strand.domains:
            self.assertIs(strand, domain.strand())

    def test_add_nick_to_2_domain_circular_strand_makes_it_linear_nick_second_domain(self) -> None:
        # nick strand 5
        r'''
           40
            strand 5
        0  /--------\
           |        |
        1  \---]<---/
        '''
        self.design.add_nick(1, 45, False)
        strand = self.design.strands[5]
        self.assertFalse(strand.circular)
        self.assertEqual(self.num_strands, len(self.design.strands))

        self.assertEqual(3, len(strand.domains))
        self.assertEqual(3, len(strand.bound_domains()))

        d0, d1, d2 = strand.bound_domains()
        self.assertEqual(1, d0.helix)
        self.assertEqual(False, d0.forward)
        self.assertEqual(40, d0.start)
        self.assertEqual(45, d0.end)

        self.assertEqual(0, d1.helix)
        self.assertEqual(True, d1.forward)
        self.assertEqual(40, d1.start)
        self.assertEqual(50, d1.end)

        self.assertEqual(1, d2.helix)
        self.assertEqual(False, d2.forward)
        self.assertEqual(45, d2.start)
        self.assertEqual(50, d2.end)

        for domain in strand.domains:
            self.assertIs(strand, domain.strand())

    def test_add_nick_to_3_domain_circular_strand_makes_it_linear_nick_first_domain(self) -> None:
        # nick strand 4
        r'''
           30         40
            strand 4
        0  [--->[---\
                     ) loopout length 5
        1  /--------/
           |
        2  \-------->cross to 5' end here
        '''
        self.design.add_nick(0, 35, True)
        strand = self.design.strands[4]
        self.assertFalse(strand.circular)
        self.assertEqual(self.num_strands, len(self.design.strands))

        self.assertEqual(5, len(strand.domains))
        self.assertEqual(4, len(strand.bound_domains()))

        d0, loopout, d1, d2, d3 = strand.domains
        self.assertIsInstance(loopout, sc.Loopout)

        self.assertEqual(0, d0.helix)
        self.assertEqual(True, d0.forward)
        self.assertEqual(35, d0.start)
        self.assertEqual(40, d0.end)

        self.assertEqual(1, d1.helix)
        self.assertEqual(False, d1.forward)
        self.assertEqual(30, d1.start)
        self.assertEqual(40, d1.end)

        self.assertEqual(2, d2.helix)
        self.assertEqual(True, d2.forward)
        self.assertEqual(30, d2.start)
        self.assertEqual(40, d2.end)

        self.assertEqual(0, d3.helix)
        self.assertEqual(True, d3.forward)
        self.assertEqual(30, d3.start)
        self.assertEqual(35, d3.end)

        for domain in strand.domains:
            self.assertIs(strand, domain.strand())

    def test_add_nick_to_3_domain_circular_strand_makes_it_linear_nick_middle_domain(self) -> None:
        # nick strand 4
        r'''
           30         40
            strand 4
        0  [--------\
                     ) loopout length 5
        1  /---]<---/
           |
        2  \-------->cross to 5' end here
        '''
        self.design.add_nick(1, 35, False)
        strand = self.design.strands[4]
        self.assertFalse(strand.circular)
        self.assertEqual(self.num_strands, len(self.design.strands))

        self.assertEqual(5, len(strand.domains))
        self.assertEqual(4, len(strand.bound_domains()))

        d0, d1, d2, loopout, d3 = strand.domains
        self.assertIsInstance(loopout, sc.Loopout)

        self.assertEqual(1, d0.helix)
        self.assertEqual(False, d0.forward)
        self.assertEqual(30, d0.start)
        self.assertEqual(35, d0.end)

        self.assertEqual(2, d1.helix)
        self.assertEqual(True, d1.forward)
        self.assertEqual(30, d1.start)
        self.assertEqual(40, d1.end)

        self.assertEqual(0, d2.helix)
        self.assertEqual(True, d2.forward)
        self.assertEqual(30, d2.start)
        self.assertEqual(40, d2.end)

        self.assertEqual(1, d3.helix)
        self.assertEqual(False, d3.forward)
        self.assertEqual(35, d3.start)
        self.assertEqual(40, d3.end)

        for domain in strand.domains:
            self.assertIs(strand, domain.strand())

    def test_add_nick_to_3_domain_circular_strand_makes_it_linear_nick_last_domain(self) -> None:
        # nick strand 4
        r'''
           30         40
            strand 4
        0  [--------\
                     ) loopout length 5
        1  /--------/
           |
        2  \--->[--->cross to 5' end here
        '''
        self.design.add_nick(2, 35, True)
        strand = self.design.strands[4]
        self.assertFalse(strand.circular)
        self.assertEqual(self.num_strands, len(self.design.strands))

        self.assertEqual(5, len(strand.domains))
        self.assertEqual(4, len(strand.bound_domains()))

        d0, d1, loopout, d2, d3 = strand.domains
        self.assertIsInstance(loopout, sc.Loopout)

        self.assertEqual(2, d0.helix)
        self.assertEqual(True, d0.forward)
        self.assertEqual(35, d0.start)
        self.assertEqual(40, d0.end)

        self.assertEqual(0, d1.helix)
        self.assertEqual(True, d1.forward)
        self.assertEqual(30, d1.start)
        self.assertEqual(40, d1.end)

        self.assertEqual(1, d2.helix)
        self.assertEqual(False, d2.forward)
        self.assertEqual(30, d2.start)
        self.assertEqual(40, d2.end)

        self.assertEqual(2, d3.helix)
        self.assertEqual(True, d3.forward)
        self.assertEqual(30, d3.start)
        self.assertEqual(35, d3.end)

        for domain in strand.domains:
            self.assertIs(strand, domain.strand())

    def test_ligate_linear_strand_to_itself_makes_it_circular(self) -> None:
        self.assertFalse(self.design.strands[1].circular)
        self.assertEqual(3, len(self.design.strands[1].domains))

        self.design.ligate(0, 15, True)

        strand = self.design.strands[1]
        self.assertEqual(self.num_strands, len(self.design.strands))
        self.assertTrue(strand.circular)
        self.assertEqual(2, len(strand.domains))

        for domain in strand.domains:
            self.assertIs(strand, domain.strand())

    # TODO: add functionality for removing crossovers and loopouts, and test that here


class TestAddStrand(unittest.TestCase):

    def test_add_strand__with_loopout(self) -> None:
        helices = [sc.Helix(max_offset=10), sc.Helix(max_offset=10)]
        design = sc.Design(helices=helices, strands=[])

        ss1 = sc.Domain(0, True, 0, 10)
        loop = sc.Loopout(4)
        ss2 = sc.Domain(1, False, 0, 10)
        strand = sc.Strand([ss1, loop, ss2])

        design.add_strand(strand)

        self.assertEqual(1, len(design.strands))
        self.assertEqual(strand, design.strands[0])
        self.assertEqual(ss1, design.domain_at(0, 0, True))
        self.assertEqual(ss2, design.domain_at(1, 0, False))

    def test_add_strand__illegal_overlapping_domains(self) -> None:
        helices = [sc.Helix(max_offset=50), sc.Helix(max_offset=50)]
        design = sc.Design(helices=helices, strands=[], grid=sc.square)
        with self.assertRaises(sc.StrandError):
            strand = sc.Strand([
                sc.Domain(0, False, 40, 48),
                sc.Domain(0, False, 32, 48, deletions=[44]),
                sc.Domain(1, True, 32, 40),
            ])
            design.add_strand(strand)


class TestAssignDNA(unittest.TestCase):

    def test_assign_dna__hairpin(self) -> None:
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
        ss_f = sc.Domain(helix=0, forward=True, start=0, end=5)
        loop = sc.Loopout(length=5)
        ss_r = sc.Domain(helix=0, forward=False, start=0, end=5)
        strand_forward = sc.Strand([ss_f, loop, ss_r])
        design = sc.Design(strands=[strand_forward], grid=sc.square)
        design.assign_dna(strand_forward, 'AAACC TGCAC')
        self.assertEqual('AAACC TGCAC GGTTT'.replace(' ', ''), strand_forward.dna_sequence)

    def test_assign_dna__from_strand_with_loopout(self) -> None:
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
        ss_f = sc.Domain(helix=0, forward=True, start=0, end=5)
        loop = sc.Loopout(length=5)
        ss_r = sc.Domain(helix=1, forward=False, start=0, end=5)
        strand_multi = sc.Strand([ss_f, loop, ss_r])

        ss_single0 = sc.Domain(helix=0, forward=False, start=0, end=5)
        strand_single0 = sc.Strand([ss_single0])

        ss_single1 = sc.Domain(helix=1, forward=True, start=0, end=5)
        strand_single1 = sc.Strand([ss_single1])

        design = sc.Design(strands=[strand_multi, strand_single0, strand_single1], grid=sc.square)

        design.assign_dna(strand_multi, 'AAACC TGCAC ATTCG')

        self.assertEqual('AAACC TGCAC ATTCG'.replace(' ', ''), strand_multi.dna_sequence)
        self.assertEqual('GGTTT'.replace(' ', ''), strand_single0.dna_sequence)
        self.assertEqual('CGAAT'.replace(' ', ''), strand_single1.dna_sequence)

    def test_assign_dna__to_strand_with_loopout(self) -> None:
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
        ss_f = sc.Domain(helix=0, forward=True, start=0, end=5)
        loop = sc.Loopout(length=5)
        ss_r = sc.Domain(helix=1, forward=False, start=0, end=5)
        strand_multi = sc.Strand([ss_f, loop, ss_r])

        ss_single0 = sc.Domain(helix=0, forward=False, start=0, end=5)
        strand_single0 = sc.Strand([ss_single0])

        ss_single1 = sc.Domain(helix=1, forward=True, start=0, end=5)
        strand_single1 = sc.Strand([ss_single1])

        design = sc.Design(strands=[strand_multi, strand_single0, strand_single1], grid=sc.square)

        design.assign_dna(strand_single0, 'GGTTT')

        self.assertEqual('AAACC ????? ?????'.replace(' ', ''), strand_multi.dna_sequence)
        self.assertEqual('GGTTT'.replace(' ', ''), strand_single0.dna_sequence)

        design.assign_dna(strand_single1, 'CGAAT')

        self.assertEqual('AAACC ????? ATTCG'.replace(' ', ''), strand_multi.dna_sequence)
        self.assertEqual('GGTTT'.replace(' ', ''), strand_single0.dna_sequence)
        self.assertEqual('CGAAT'.replace(' ', ''), strand_single1.dna_sequence)

    def test_assign_dna__assign_from_strand_multi_other_single(self) -> None:
        """
          01234567
        <-TTTG----GACA-]
        +-AAAC->[-CTGT-+   # helix 0
        |              |
        +-GCTT----AGTA-+   # helix 1
        """
        ss_f_left = sc.Domain(helix=0, forward=True, start=0, end=4)
        ss_f_right = sc.Domain(helix=0, forward=True, start=4, end=8)
        ss_h1 = sc.Domain(helix=1, forward=False, start=0, end=8)
        strand_multi = sc.Strand([ss_f_right, ss_h1, ss_f_left])

        ss_single = sc.Domain(helix=0, forward=False, start=0, end=8)
        strand_single = sc.Strand([ss_single])

        design = sc.Design(strands=[strand_multi, strand_single], grid=sc.square)

        design.assign_dna(strand_multi, 'CTGT ATGA TTCG AAAC')

        self.assertEqual('CTGT ATGA TTCG AAAC'.replace(' ', ''), strand_multi.dna_sequence)
        self.assertEqual('ACAG GTTT'.replace(' ', ''), strand_single.dna_sequence)

    def test_assign_dna__assign_to_strand_multi_other_single(self) -> None:
        """
          01234567
        <-TTTG----GACA-]
        +-AAAC->[-CTGT-+   # helix 0
        |              |
        +-????----????-+   # helix 1
        """
        ss_f_left = sc.Domain(helix=0, forward=True, start=0, end=4)
        ss_f_right = sc.Domain(helix=0, forward=True, start=4, end=8)
        ss_h1 = sc.Domain(helix=1, forward=False, start=0, end=8)
        strand_multi = sc.Strand([ss_f_right, ss_h1, ss_f_left])

        ss_single = sc.Domain(helix=0, forward=False, start=0, end=8)
        strand_single = sc.Strand([ss_single])

        design = sc.Design(strands=[strand_multi, strand_single], grid=sc.square)

        design.assign_dna(strand_single, 'ACAG GTTT')

        self.assertEqual('CTGT ???? ???? AAAC'.replace(' ', ''), strand_multi.dna_sequence)
        self.assertEqual('ACAG GTTT'.replace(' ', ''), strand_single.dna_sequence)

    def test_assign_dna__other_strand_fully_defined_already(self) -> None:
        """
        01234567
        [------>
        CAAAGTCG
        GTTT
        <--]
        """
        ss_r = sc.Domain(helix=0, forward=True, start=0, end=8)
        ss_l = sc.Domain(helix=0, forward=False, start=0, end=4)
        strand_r = sc.Strand(domains=[ss_r])
        strand_l = sc.Strand(domains=[ss_l])
        design = sc.Design(grid=sc.square, strands=[strand_r, strand_l])
        design.assign_dna(strand_r, 'CAAAGTCG')
        design.assign_dna(strand_l, 'TTTG')
        # should not have an error by this point

    def test_assign_dna__other_strand_fully_defined_already_and_other_extends_beyond(self) -> None:
        """
        01234567
        [------>
        CAAAGTCG
          TTCA
          <--]
        """
        ss_r = sc.Domain(helix=0, forward=True, start=0, end=8)
        ss_l = sc.Domain(helix=0, forward=False, start=2, end=6)
        strand_r = sc.Strand(domains=[ss_r])
        strand_l = sc.Strand(domains=[ss_l])
        design = sc.Design(grid=sc.square, strands=[strand_r, strand_l])
        design.assign_dna(strand_r, 'CAAAGTCG')
        design.assign_dna(strand_l, 'ACTT')
        # should not have an error by this point

    def test_assign_dna__other_strand_fully_defined_already_and_self_extends_beyond(self) -> None:
        """
        01234567
        [------>
        CAAAGTCG
          TTCA
          <--]
        """
        ss_r = sc.Domain(helix=0, forward=True, start=0, end=8)
        ss_l = sc.Domain(helix=0, forward=False, start=2, end=6)
        strand_r = sc.Strand(domains=[ss_r])
        strand_l = sc.Strand(domains=[ss_l])
        design = sc.Design(grid=sc.square, strands=[strand_r, strand_l])
        design.assign_dna(strand_l, 'ACTT')
        design.assign_dna(strand_r, 'CAAAGTCG')
        # should not have an error by this point

    def test_assign_dna__two_equal_length_strands_on_one_helix(self) -> None:
        """
        01234
        <---]
        CAAAA
        GTTTT
        [--->
        """
        ss_r = sc.Domain(helix=0, forward=True, start=0, end=5)
        ss_l = sc.Domain(helix=0, forward=False, start=0, end=5)
        strand_r = sc.Strand(domains=[ss_r])
        strand_l = sc.Strand(domains=[ss_l])
        design = sc.Design(grid=sc.square, strands=[strand_r, strand_l])
        design.assign_dna(strand_l, 'AAAAC')
        self.assertEqual('GTTTT', strand_r.dna_sequence)

    def test_assign_dna__assign_seq_with_wildcards(self) -> None:
        """
        01234
        <---]
        C??AA
        G??TT
        [--->
        """
        ss_bot = sc.Domain(helix=0, forward=True, start=0, end=5)
        ss_top = sc.Domain(helix=0, forward=False, start=0, end=5)
        strand_bot = sc.Strand(domains=[ss_bot])
        strand_top = sc.Strand(domains=[ss_top])
        strands = [strand_bot, strand_top]
        design = sc.Design(grid=sc.square, strands=strands)
        design.assign_dna(strand_top, 'AA??C')
        self.assertEqual('G??TT', strand_bot.dna_sequence)

    def test_assign_dna__one_strand_assigned_by_complement_from_two_other_strands(self) -> None:
        """
          0123     4567
        <-AAAC-] <-GGGA-]
        [-TTTG-----CCCT->
        """
        ss_top_left = sc.Domain(0, False, 0, 4)
        ss_top_right = sc.Domain(0, False, 4, 8)
        ss_bot = sc.Domain(0, True, 0, 8)
        st_top_left = sc.Strand([ss_top_left])
        st_top_right = sc.Strand([ss_top_right])
        st_bot = sc.Strand([ss_bot])
        design = sc.Design(strands=[st_bot, st_top_left, st_top_right], grid=sc.square)
        design.assign_dna(st_top_left, 'CAAA')
        self.assertEqual('TTTG????', st_bot.dna_sequence)
        design.assign_dna(st_top_right, 'AGGG')
        self.assertEqual('TTTGCCCT', st_bot.dna_sequence)

    def test_assign_dna__adapter_assigned_from_scaffold_and_tiles(self) -> None:
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
        scaf0_ss = sc.Domain(0, False, 0, 6)
        scaf1_ss = sc.Domain(1, True, 0, 6)
        tile1_ss = sc.Domain(1, True, 6, 12)
        tile0_ss = sc.Domain(0, False, 6, 12)
        adap0_ss = sc.Domain(0, True, 2, 10)
        adap1_ss = sc.Domain(1, False, 2, 10)
        scaf = sc.Strand([scaf1_ss, scaf0_ss])
        adap = sc.Strand([adap0_ss, adap1_ss])
        tile0 = sc.Strand([tile0_ss])
        tile1 = sc.Strand([tile1_ss])

        design = sc.Design(strands=[scaf, adap, tile0, tile1])

        design.assign_dna(tile0, 'AA AATG')
        self.assertEqual('???? CATT ???? ????'.replace(' ', ''), adap.dna_sequence)

        design.assign_dna(tile1, 'TGCC GG')
        self.assertEqual('???? CATT GGCA ????'.replace(' ', ''), adap.dna_sequence)

        design.assign_dna(scaf, 'AA TTTG GAAA TG')
        self.assertEqual('TTTC CATT GGCA CAAA'.replace(' ', ''), adap.dna_sequence)

    def test_assign_dna__adapter_assigned_from_scaffold_and_tiles_with_deletions(self) -> None:
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
        scaf0_ss = sc.Domain(0, False, 0, 6)
        scaf1_ss = sc.Domain(1, True, 0, 6)
        tile1_ss = sc.Domain(1, True, 6, 12)
        tile0_ss = sc.Domain(0, False, 6, 12)
        adap0_ss = sc.Domain(0, True, 2, 10)
        adap1_ss = sc.Domain(1, False, 2, 10)
        scaf = sc.Strand([scaf1_ss, scaf0_ss])
        adap = sc.Strand([adap0_ss, adap1_ss])
        tile0 = sc.Strand([tile0_ss])
        tile1 = sc.Strand([tile1_ss])

        design = sc.Design(strands=[scaf, adap, tile0, tile1])
        design.add_deletion(0, 3)
        design.add_deletion(0, 8)
        design.add_deletion(1, 8)

        design.assign_dna(tile0, 'AA ATG')
        self.assertEqual('??? CAT ??? ????'.replace(' ', ''), adap.dna_sequence)

        design.assign_dna(tile1, 'TGC GG')
        self.assertEqual('??? CAT GCA ????'.replace(' ', ''), adap.dna_sequence)

        design.assign_dna(scaf, 'AA TTTG GAA TG')
        self.assertEqual('TTC CAT GCA CAAA'.replace(' ', ''), adap.dna_sequence)

    def test_assign_dna__adapter_assigned_from_scaffold_and_tiles_with_insertions(self) -> None:
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
        scaf0_ss = sc.Domain(0, False, 0, 6)
        scaf1_ss = sc.Domain(1, True, 0, 6)
        tile1_ss = sc.Domain(1, True, 6, 12)
        tile0_ss = sc.Domain(0, False, 6, 12)
        adap0_ss = sc.Domain(0, True, 2, 10)
        adap1_ss = sc.Domain(1, False, 2, 10)
        scaf = sc.Strand([scaf1_ss, scaf0_ss])
        adap = sc.Strand([adap0_ss, adap1_ss])
        tile0 = sc.Strand([tile0_ss])
        tile1 = sc.Strand([tile1_ss])

        design = sc.Design(strands=[scaf, adap, tile0, tile1])
        design.add_insertion(0, 8, 1)
        design.add_insertion(1, 8, 1)

        design.assign_dna(tile0, 'AA AAATG')
        self.assertEqual('???? CATTT ????? ????'.replace(' ', ''), adap.dna_sequence)

        design.assign_dna(tile1, 'TGCCC GG')
        self.assertEqual('???? CATTT GGGCA ????'.replace(' ', ''), adap.dna_sequence)

        design.assign_dna(scaf, 'AA TTTG GAAA TG')
        self.assertEqual('TTTC CATTT GGGCA CAAA'.replace(' ', ''), adap.dna_sequence)

    def test_assign_dna__dna_sequence_shorter_than_complementary_strand_right_strand_longer(self) -> None:
        """
        <---]
        CAAAA
        GTTTT?????
        [-------->
        """
        ss_long = sc.Domain(helix=0, forward=True, start=0, end=10)
        ss_short = sc.Domain(helix=0, forward=False, start=0, end=5)
        strand_long = sc.Strand(domains=[ss_long])
        strand_short = sc.Strand(domains=[ss_short])
        strands = [strand_long, strand_short]
        design = sc.Design(grid=sc.square, strands=strands)
        design.assign_dna(strand_short, 'AAAAC')
        self.assertEqual('GTTTT?????', strand_long.dna_sequence)

    def test_assign_dna__dna_sequence_shorter_than_complementary_strand_left_strand_longer(self) -> None:
        """
        [--->
        AAAAC
        TTTTG?????
        <--------]
        """
        ss_long = sc.Domain(helix=0, forward=False, start=0, end=10)
        ss_short = sc.Domain(helix=0, forward=True, start=0, end=5)
        strand_long = sc.Strand(domains=[ss_long])
        strand_short = sc.Strand(domains=[ss_short])
        strands = [strand_long, strand_short]
        design = sc.Design(grid=sc.square, strands=strands)
        design.assign_dna(strand_short, 'AAAAC')
        self.assertEqual('?????GTTTT', strand_long.dna_sequence)

    def test_assign_dna__dna_sequence_with_uncomplemented_domain_on_different_helix(self) -> None:
        """
        <---]
        CAAAA
        GTTTT?????
        [--------+
                 |
               <-+
               ???
        """
        ss_long = sc.Domain(helix=0, forward=True, start=0, end=10)
        ss_long_h1 = sc.Domain(helix=0, forward=False, start=7, end=10)
        ss_short = sc.Domain(helix=0, forward=False, start=0, end=5)
        strand_long = sc.Strand(domains=[ss_long, ss_long_h1])
        strand_short = sc.Strand(domains=[ss_short])
        strands = [strand_long, strand_short]
        design = sc.Design(grid=sc.square, strands=strands)
        design.assign_dna(strand_short, 'AAAAC')
        self.assertEqual('GTTTT????????', strand_long.dna_sequence)

    def test_assign_dna__dna_sequence_with_uncomplemented_domain_on_different_helix_wildcards_both_ends(
            self) -> None:
        """
             <---]
             CAAAA
        ?????GTTTT
        [--------+
                 |
               <-+
               ???
        """
        ss_long_h0 = sc.Domain(helix=0, forward=True, start=0, end=10)
        ss_long_h1 = sc.Domain(helix=1, forward=False, start=7, end=10)
        ss_short_h0 = sc.Domain(helix=0, forward=False, start=5, end=10)
        strand_long = sc.Strand(domains=[ss_long_h0, ss_long_h1])
        strand_short = sc.Strand(domains=[ss_short_h0])
        strands = [strand_long, strand_short]
        design = sc.Design(grid=sc.square, strands=strands)
        design.assign_dna(strand_short, 'AAAAC')
        self.assertEqual('?????GTTTT???', strand_long.dna_sequence)

    def test_assign_dna__one_helix_with_one_bottom_strand_and_three_top_strands(self) -> None:
        """
         012   345   678
        -TTT> -GGG> -CCC>
        <AAA---CCC---GGG-
         876   543   210
        """
        ss_bot = sc.Domain(helix=0, forward=False, start=0, end=9)
        ss_top1 = sc.Domain(helix=0, forward=True, start=0, end=3)
        ss_top2 = sc.Domain(helix=0, forward=True, start=3, end=6)
        ss_top3 = sc.Domain(helix=0, forward=True, start=6, end=9)
        strand_bot = sc.Strand(domains=[ss_bot])
        strand_top1 = sc.Strand(domains=[ss_top1])
        strand_top2 = sc.Strand(domains=[ss_top2])
        strand_top3 = sc.Strand(domains=[ss_top3])
        strands = [strand_bot, strand_top1, strand_top2, strand_top3]
        design = sc.Design(grid=sc.square, strands=strands)
        design.assign_dna(strand_bot, 'AAACCCGGG')
        self.assertEqual('CCC', strand_top1.dna_sequence)
        self.assertEqual('GGG', strand_top2.dna_sequence)
        self.assertEqual('TTT', strand_top3.dna_sequence)

    def test_assign_dna__two_helices_with_multiple_domain_intersections(self) -> None:
        """
                012    345   678    901
        M13   [-ACC----TAA---GAA----AAC---+
              +-TGG-]<-ATT-+ CTT----TTG-+ |
              |            | |          | |
              +-GAT----TTC-+ ATG->[-AGT-+ |
              <-CTA----AAG---TAC----TCA---+
        """
        scaf0_ss = sc.Domain(helix=0, forward=True, start=0, end=12)
        scaf1_ss = sc.Domain(helix=1, forward=False, start=0, end=12)
        scaf = sc.Strand(domains=[scaf0_ss, scaf1_ss])

        first_stap0_left_ss = sc.Domain(helix=0, forward=False, start=0, end=3)
        first_stap1_ss = sc.Domain(helix=1, forward=True, start=0, end=6)
        first_stap0_right_ss = sc.Domain(helix=0, forward=False, start=3, end=6)
        first_stap = sc.Strand(domains=[first_stap0_left_ss, first_stap1_ss, first_stap0_right_ss])

        second_stap1_right_ss = sc.Domain(helix=1, forward=True, start=9, end=12)
        second_stap0_ss = sc.Domain(helix=0, forward=False, start=6, end=12)
        second_stap1_left_ss = sc.Domain(helix=1, forward=True, start=6, end=9)
        second_stap = sc.Strand(domains=[second_stap1_right_ss, second_stap0_ss, second_stap1_left_ss])

        strands = [scaf, first_stap, second_stap]
        design = sc.Design(grid=sc.square, strands=strands)
        design.assign_dna(scaf, 'ACC TAA GAA AAC ACT CAT GAA ATC'.replace(' ', ''))
        self.assertEqual('GGT GAT TTC TTA'.replace(' ', ''), first_stap.dna_sequence)
        self.assertEqual('AGT GTT TTC ATG'.replace(' ', ''), second_stap.dna_sequence)

    def test_assign_dna__upper_left_edge_staple_of_16H_origami_rectangle(self) -> None:
        """
        staple <ACATAAGAAAACGGAG--+
        M13   +-TGTATTCTTTTGCCTC> |
              |                   |
              +-GATTTTGTGAGTAGAA- |
               -CTAAAACACTCATCTT--+
        """
        scaf0_ss = sc.Domain(helix=0, forward=True, start=0, end=16)
        scaf1_ss = sc.Domain(helix=1, forward=False, start=0, end=16)
        stap0_ss = sc.Domain(helix=0, forward=False, start=0, end=16)
        stap1_ss = sc.Domain(helix=1, forward=True, start=0, end=16)
        scaf = sc.Strand(domains=[scaf1_ss, scaf0_ss])
        stap = sc.Strand(domains=[stap1_ss, stap0_ss])
        strands = [scaf, stap]
        design = sc.Design(grid=sc.square, strands=strands)

        seq_m13_upper_left = 'AAGATGAGTGTTTTAGTGTATTCTTTTGCCTC'
        design.assign_dna(scaf, seq_m13_upper_left)
        expected_seq_stap_upperleft = 'CTAAAACACTCATCTTGAGGCAAAAGAATACA'
        self.assertEqual(expected_seq_stap_upperleft, stap.dna_sequence)

    def test_assign_dna__2helix_with_deletions(self) -> None:
        r"""
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
        stap_left_ss1 = sc.Domain(1, True, 0, width_h)
        stap_left_ss0 = sc.Domain(0, False, 0, width_h)
        stap_right_ss0 = sc.Domain(0, False, width_h, width)
        stap_right_ss1 = sc.Domain(1, True, width_h, width)
        scaf_ss1_left = sc.Domain(1, False, 0, width_h)
        scaf_ss0 = sc.Domain(0, True, 0, width)
        scaf_ss1_right = sc.Domain(1, False, width_h, width)
        stap_left = sc.Strand([stap_left_ss1, stap_left_ss0])
        stap_right = sc.Strand([stap_right_ss0, stap_right_ss1])
        scaf = sc.Strand([scaf_ss1_left, scaf_ss0, scaf_ss1_right], color=sc.default_scaffold_color)
        strands = [stap_left, stap_right, scaf]
        design = sc.Design(helices=helices, strands=strands, grid=sc.square)
        design.add_deletion(helix=0, offset=1)
        design.add_deletion(helix=0, offset=4)
        design.add_deletion(helix=1, offset=1)
        design.add_deletion(helix=1, offset=4)
        design.assign_dna(scaf, 'AACATCGT')
        self.assertEqual("AACATCGT", scaf.dna_sequence)
        self.assertEqual("TTTG", stap_left.dna_sequence)
        self.assertEqual("GAAC", stap_right.dna_sequence)

    def test_assign_dna__wildcards_simple(self) -> None:
        """
         012   345   678
        -TTC> -GGA> -CCT>
        <AAG---CCT---GGA-
         876   543   210
        """
        ss_bot = sc.Domain(helix=0, forward=False, start=0, end=9)
        ss_top1 = sc.Domain(helix=0, forward=True, start=0, end=3)
        ss_top2 = sc.Domain(helix=0, forward=True, start=3, end=6)
        ss_top3 = sc.Domain(helix=0, forward=True, start=6, end=9)
        strand_bot = sc.Strand(domains=[ss_bot])
        strand_top1 = sc.Strand(domains=[ss_top1])
        strand_top2 = sc.Strand(domains=[ss_top2])
        strand_top3 = sc.Strand(domains=[ss_top3])
        strands = [strand_bot, strand_top1, strand_top2, strand_top3]
        design = sc.Design(grid=sc.square, strands=strands)

        design.assign_dna(strand_top1, 'TTC')
        self.assertEqual('??????GAA', strand_bot.dna_sequence)

        design.assign_dna(strand_top3, 'CCT')
        self.assertEqual('AGG???GAA', strand_bot.dna_sequence)

        design.assign_dna(strand_top2, 'GGA')
        self.assertEqual('AGGTCCGAA', strand_bot.dna_sequence)


class TestAssignDNAToDomains(unittest.TestCase):

    def setUp(self) -> None:
        """
         012   345   678   901   234   567   890
              +---------------+
              |               |
              |         +-----|-------+   +-+
              |         |     |       |   | |
        -ACG> +TTC> -GGA+ -AAC+ -TGC> +TTG+ +GCA>
        <TGC---AAG---CCT---TTG---ACG---AAC---CGT-
         098   765   432   109   876   543   210
        """
        self.dom_bot = sc.Domain(helix=0, forward=False, start=0, end=21)
        self.dom_top0 = sc.Domain(helix=0, forward=True, start=0, end=3)
        self.dom_top3 = sc.Domain(helix=0, forward=True, start=3, end=6)
        self.dom_top6 = sc.Domain(helix=0, forward=True, start=6, end=9)
        self.dom_top9 = sc.Domain(helix=0, forward=True, start=9, end=12)
        self.dom_top12 = sc.Domain(helix=0, forward=True, start=12, end=15)
        self.dom_top15 = sc.Domain(helix=0, forward=True, start=15, end=18)
        self.dom_top18 = sc.Domain(helix=0, forward=True, start=18, end=21)
        self.strand_bot = sc.Strand(domains=[self.dom_bot])
        self.strand_top_small0 = sc.Strand(domains=[self.dom_top0])
        self.strand_top_small12 = sc.Strand(domains=[self.dom_top12])
        self.strand_top_big9 = sc.Strand(domains=[self.dom_top9, self.dom_top3])
        self.strand_top_big6 = sc.Strand(domains=[self.dom_top6, self.dom_top15, self.dom_top18])
        strands = [self.strand_bot, self.strand_top_small0, self.strand_top_small12,
                   self.strand_top_big9, self.strand_top_big6]
        self.design = sc.Design(grid=sc.square, strands=strands)

    def test_assign_dna__wildcards_multiple_overlaps(self) -> None:
        """
         012   345   678   901   234   567   890
              +---------------+
              |               |
              |         +-----|-------+   +-+
              |         |     |       |   | |
        -ACG> +TTC> -GGA+ -AAC+ -TGC> +TTG+ +GCA>
        <TGC---AAG---CCT---TTG---ACG---AAC---CGT-
         098   765   432   109   876   543   210
        """
        self.design.assign_dna(self.strand_top_big9, 'AACTTC')
        self.assertEqual('??? ??? ??? GTT ??? GAA ???'.replace(' ', ''), self.strand_bot.dna_sequence)

        self.design.assign_dna(self.strand_top_small12, 'TGC')
        self.assertEqual('??? ??? GCA GTT ??? GAA ???'.replace(' ', ''), self.strand_bot.dna_sequence)

        self.design.assign_dna(self.strand_top_small0, 'ACG')
        self.assertEqual('??? ??? GCA GTT ??? GAA CGT'.replace(' ', ''), self.strand_bot.dna_sequence)

        self.design.assign_dna(self.strand_top_big6, 'GGATTGGCA')
        self.assertEqual('TGC CAA GCA GTT TCC GAA CGT'.replace(' ', ''), self.strand_bot.dna_sequence)

    def test_assign_dna__domain_sequence_too_long_error(self) -> None:
        with self.assertRaises(sc.IllegalDesignError):
            self.design.assign_dna(self.strand_top_big9, 'AACTTC', domain=self.dom_top9)

    def test_assign_dna__to_individual_domains__wildcards_multiple_overlaps(self) -> None:
        """
         012   345   678   901   234   567   890
              +---------------+
              |               |
              |         +-----|-------+   +-+
              |         |     |       |   | |
        -ACG> +TTC> -GGA+ -AAC+ -TGC> +TTG+ +GCA>
        <TGC---AAG---CCT---TTG---ACG---AAC---CGT-
         098   765   432   109   876   543   210
        """

        """
         012   345   678   901   234   567   890
              +---------------+
              |               |
              |         +-----|-------+   +-+
              |         |     |       |   | |
        -???> +???> -???+ -AAC+ -???> +???+ +???>
        <???---???---???---TTG---???---???---???-
         098   765   432   109   876   543   210
        """
        self.design.assign_dna(self.strand_top_big9, 'AAC', domain=self.dom_top9)
        self.assertEqual('AAC ???'.replace(' ', ''), self.strand_top_big9.dna_sequence)
        self.assertEqual('??? ??? ??? GTT ??? ??? ???'.replace(' ', ''), self.strand_bot.dna_sequence)

        """
         012   345   678   901   234   567   890
              +---------------+
              |               |
              |         +-----|-------+   +-+
              |         |     |       |   | |
        -???> +TTC> -???+ -AAC+ -???> +???+ +???>
        <???---AAG---???---TTG---???---???---???-
         098   765   432   109   876   543   210
        """
        self.design.assign_dna(self.strand_top_big9, 'TTC', domain=self.dom_top3)
        self.assertEqual('AAC TTC'.replace(' ', ''), self.strand_top_big9.dna_sequence)
        self.assertEqual('??? ??? ??? GTT ??? GAA ???'.replace(' ', ''), self.strand_bot.dna_sequence)

        """
         012   345   678   901   234   567   890
              +---------------+
              |               |
              |         +-----|-------+   +-+
              |         |     |       |   | |
        -???> +TTC> -???+ -AAC+ -TGC> +???+ +???>
        <???---AAG---???---TTG---ACG---???---???-
         098   765   432   109   876   543   210
        """
        self.design.assign_dna(self.strand_top_small12, 'TGC')
        self.assertEqual('TGC', self.strand_top_small12.dna_sequence)
        self.assertEqual('??? ??? GCA GTT ??? GAA ???'.replace(' ', ''), self.strand_bot.dna_sequence)

        """
         012   345   678   901   234   567   890
              +---------------+
              |               |
              |         +-----|-------+   +-+
              |         |     |       |   | |
        -ACG> +TTC> -???+ -AAC+ -TGC> +???+ +???>
        <TGC---AAG---???---TTG---ACG---???---???-
         098   765   432   109   876   543   210
        """
        self.design.assign_dna(self.strand_top_small0, 'ACG')
        self.assertEqual('ACG', self.strand_top_small0.dna_sequence)
        self.assertEqual('??? ??? GCA GTT ??? GAA CGT'.replace(' ', ''), self.strand_bot.dna_sequence)

        """
         012   345   678   901   234   567   890
              +---------------+
              |               |
              |         +-----|-------+   +-+
              |         |     |       |   | |
        -ACG> +TTC> -???+ -AAC+ -TGC> +TTG+ +???>
        <TGC---AAG---???---TTG---ACG---AAC---???-
         098   765   432   109   876   543   210
        """
        self.design.assign_dna(self.strand_top_big6, 'TTG', domain=self.dom_top15)
        self.assertEqual('??? TTG ???'.replace(' ', ''), self.strand_top_big6.dna_sequence)
        self.assertEqual('??? CAA GCA GTT ??? GAA CGT'.replace(' ', ''), self.strand_bot.dna_sequence)

        """
         012   345   678   901   234   567   890
              +---------------+
              |               |
              |         +-----|-------+   +-+
              |         |     |       |   | |
        -ACG> +TTC> -???+ -AAC+ -TGC> +TTG+ +GCA>
        <TGC---AAG---???---TTG---ACG---AAC---CGT-
         098   765   432   109   876   543   210
        """
        self.design.assign_dna(self.strand_top_big6, 'GCA', domain=self.dom_top18)
        self.assertEqual('??? TTG GCA'.replace(' ', ''), self.strand_top_big6.dna_sequence)
        self.assertEqual('TGC CAA GCA GTT ??? GAA CGT'.replace(' ', ''), self.strand_bot.dna_sequence)

        """
         012   345   678   901   234   567   890
              +---------------+
              |               |
              |         +-----|-------+   +-+
              |         |     |       |   | |
        -ACG> +TTC> -GGA+ -AAC+ -TGC> +TTG+ +GCA>
        <TGC---AAG---CCT---TTG---ACG---AAC---CGT-
         098   765   432   109   876   543   210
        """
        self.design.assign_dna(self.strand_top_big6, 'GGA', domain=self.dom_top6)
        self.assertEqual('GGA TTG GCA'.replace(' ', ''), self.strand_top_big6.dna_sequence)
        self.assertEqual('TGC CAA GCA GTT TCC GAA CGT'.replace(' ', ''), self.strand_bot.dna_sequence)

    def test_method_chaining_with_domain_sequence(self) -> None:
        """
         012   345   678   901   234   567   890
              +---------------+
              |               |
              |         +-----|-------+   +-+
              |         |     |       |   | |
        -ACG> +TTC> -GGA+ -AAC+ -TGC> +TTG+ +GCA>
        <TGC---AAG---CCT---TTG---ACG---AAC---CGT-
         098   765   432   109   876   543   210
        """
        design = sc.Design(grid=sc.square, helices=[sc.Helix(max_offset=100)], strands=[])

        """
         012   345   678   901   234   567   890
        
        <???---???---???---???---???---???---???-
         098   765   432   109   876   543   210
        """
        design.strand(0, 21).to(0)
        self.assertEqual(1, len(design.strands))
        self.assertEqual(21, design.strands[0].domains[0].dna_length())

        """
         012   345   678   901   234   567   890
        
                                -TGC>
        <???---???---???---???---ACG---???---???-
         098   765   432   109   876   543   210
        """
        design.strand(0, 12).to(15).with_sequence('TGC')
        self.assertEqual('TGC'.replace(' ', ''), design.strands[-1].dna_sequence)
        self.assertEqual('??? ??? GCA ??? ??? ??? ???'.replace(' ', ''), design.strands[0].dna_sequence)

        """
         012   345   678   901   234   567   890
        
        -ACG>                   -TGC>
        <???---???---???---???---ACG---???---???-
         098   765   432   109   876   543   210
        """
        design.strand(0, 0).to(3).with_sequence('ACG', assign_complement=False)
        self.assertEqual('ACG'.replace(' ', ''), design.strands[-1].dna_sequence)
        self.assertEqual('??? ??? GCA ??? ??? ??? ???'.replace(' ', ''), design.strands[0].dna_sequence)

        """
         012   345   678   901   234   567   890
              +---------------+
              |               |
              |               |
              |               |
        -ACG> +TTC>       -AAC+ -TGC>
        <???---AAG---???---TTG---ACG---???---???-
         098   765   432   109   876   543   210
        """
        sb = design.strand(0, 9).to(12)
        sb.with_domain_sequence('AAC')
        sb.cross(0, offset=3)
        sb.to(6)
        sb.with_domain_sequence('TTC')
        self.assertEqual('AAC TTC'.replace(' ', ''), design.strands[-1].dna_sequence)
        self.assertEqual('??? ??? GCA GTT ??? GAA ???'.replace(' ', ''), design.strands[0].dna_sequence)

        """
         012   345   678   901   234   567   890
              +---------------+
              |               |
              |         +-----|-------+   +-+
              |         |     |       |   | |
        -ACG> +TTC> -GGA+ -AAC+ -TGC> +TTG+ +GCA>
        <???---AAG---CCT---TTG---ACG---AAC---CGT-
         098   765   432   109   876   543   210
        """
        design.strand(0, 6).to(9) \
            .with_domain_sequence('GGA') \
            .cross(0, offset=15) \
            .to(18) \
            .with_domain_sequence('TTG') \
            .to(21) \
            .with_domain_sequence('GCA')
        self.assertEqual('GGA TTG GCA'.replace(' ', ''), design.strands[-1].dna_sequence)
        self.assertEqual('TGC CAA GCA GTT TCC GAA ???'.replace(' ', ''), design.strands[0].dna_sequence)


TEST_OFFSETS_AT_DELETION_INSERTIONS = False


class TestSubstrandDNASequenceIn(unittest.TestCase):

    def test_dna_sequence_in__right_then_left(self) -> None:
        ss0 = sc.Domain(0, True, 0, 10)
        ss1 = sc.Domain(1, False, 0, 10)
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

    def test_dna_sequence_in__right_then_left_deletions(self) -> None:
        ss0 = sc.Domain(0, True, 0, 10, deletions=[2, 5, 6])
        ss1 = sc.Domain(1, False, 0, 10, deletions=[2, 6, 7])
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

    def test_dna_sequence_in__right_then_left_insertions(self) -> None:
        ss0 = sc.Domain(0, True, 0, 10, insertions=[(2, 1), (6, 2)])
        ss1 = sc.Domain(1, False, 0, 10, insertions=[(2, 1), (6, 2)])
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

    def test_dna_sequence_in__right_then_left_deletions_and_insertions(self) -> None:
        ss0 = sc.Domain(0, True, 0, 10, deletions=[4], insertions=[(2, 1), (6, 2)])
        ss1 = sc.Domain(1, False, 0, 10, deletions=[4], insertions=[(2, 1), (6, 2)])
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


class TestOxdnaExport(unittest.TestCase):
    def setUp(self) -> None:
        self.OX_UNITS_TO_NM = 0.8518
        self.NM_TO_OX_UNITS = 1.0 / self.OX_UNITS_TO_NM
        self.OX_BASE_DIST = 0.6
        self.BASES_PER_TURN = 10.5
        self.HELIX_ANGLE = math.pi * 2 / self.BASES_PER_TURN
        self.RISE_PER_BASE_PAIR = 0.332
        # square of expected distance between adjacent nucleotide centers of mass
        self.EXPECTED_ADJ_NUC_CM_DIST2 = (2 * self.OX_BASE_DIST * math.sin(self.HELIX_ANGLE / 2)) ** 2 + (
                self.RISE_PER_BASE_PAIR * self.NM_TO_OX_UNITS) ** 2

    def test_basic_design(self) -> None:
        """ 2 double strands of length 7 connected across helices.
                  0      7
        helix 0   [------+
                  +------]
                  |      |
        helix 1   <------+
                  +------>
        """
        helices = [sc.Helix(max_offset=7), sc.Helix(max_offset=7)]
        design = sc.Design(helices=helices, grid=sc.square)
        design.strand(0, 0).to(7).cross(1).move(-7)
        design.strand(0, 7).move(-7).cross(1).move(7)

        # expected values for verification
        expected_num_nucleotides = 7 * 4
        expected_strand_length = 7 * 2

        dat, top = design.to_oxdna_format()
        dat_lines = dat.strip().split('\n')
        top_lines = top.strip().split('\n')

        # check length of output files are as expected (matches # of nucleotides plus header size)
        self.assertEqual(expected_num_nucleotides + 3, len(dat_lines))
        self.assertEqual(expected_num_nucleotides + 1, len(top_lines))

        # find relevant values for nucleotides
        cm_poss = []  # center of mass position
        nbrs_3p = []
        nbrs_5p = []

        for line in dat_lines[3:]:
            data = line.strip().split()
            # make sure there are 15 values per line (3 values per vector * 5 vectors per line)
            # order of vectors: center of mass position, backbone base versor, normal versor, velocity, angular velocity (more info on versors: https://eater.net/quaternions)
            self.assertEqual(15, len(data))

            cm_poss.append(tuple([float(x) for x in data[0:3]]))
            bb_vec = tuple([float(x) for x in data[3:6]])  # backbone base vector
            nm_vec = tuple([float(x) for x in data[6:9]])  # normal vector

            # make sure normal vectors and backbone vectors are unit length
            sqr_bb_vec = sum([x ** 2 for x in bb_vec])
            sqr_nm_vec = sum([x ** 2 for x in nm_vec])
            self.assertAlmostEqual(1.0, sqr_bb_vec)
            self.assertAlmostEqual(1.0, sqr_nm_vec)

            for value in data[9:]:  # values for velocity and angular velocity vectors are 0
                self.assertAlmostEqual(0, float(value))

        strand1_idxs = []
        strand2_idxs = []
        for nuc_idx, line in enumerate(top_lines[1:]):
            data = line.strip().split()
            # make sure there are 4 values per line: strand, base, 3' neighbor, 5' neighbor
            self.assertEqual(4, len(data))

            # make sure there are only 2 strands
            strand_num = int(data[0])
            self.assertIn(strand_num, [1, 2])
            # make sure base is valid
            base = data[1]
            self.assertIn(base, ['A', 'C', 'G', 'T'])

            nbrs_3p.append(int(data[2]))
            nbrs_5p.append(int(data[3]))

            # append start of strand (no 5' neighbor) to list of indexes for strand
            neighbor_5 = int(data[3])
            if neighbor_5 == -1:
                if strand_num == 1:
                    strand1_start = nuc_idx
                    strand1_idxs.append(strand1_start)
                else:
                    strand2_start = nuc_idx
                    strand2_idxs.append(strand2_start)

        # reconstruct strands using indices from oxDNA files
        next_idx = nbrs_3p[strand1_start]
        while next_idx >= 0:
            strand1_idxs.append(next_idx)
            next_idx = nbrs_3p[strand1_idxs[-1]]

        next_idx = nbrs_3p[strand2_start]
        while next_idx >= 0:
            strand2_idxs.append(next_idx)
            next_idx = nbrs_3p[strand2_idxs[-1]]

        # assert that strands are the correct length
        self.assertEqual(expected_strand_length, len(strand1_idxs))
        self.assertEqual(expected_strand_length, len(strand2_idxs))

        for i in range(expected_strand_length - 1):
            # ignore nucleotide distance between domains (on crossover)
            if i == 6:
                continue

            strand1_nuc_idx1 = strand1_idxs[i]
            strand1_nuc_idx2 = strand1_idxs[i + 1]
            strand2_nuc_idx1 = strand2_idxs[i]
            strand2_nuc_idx2 = strand2_idxs[i + 1]

            # find the center of mass for adjacent nucleotides
            s1_cmp1 = cm_poss[strand1_nuc_idx1]
            s1_cmp2 = cm_poss[strand1_nuc_idx2]
            s2_cmp1 = cm_poss[strand2_nuc_idx1]
            s2_cmp2 = cm_poss[strand2_nuc_idx2]

            # calculate and verify squared distance between adjacent nucleotides in a domain
            diff1 = tuple([s1_cmp1[j] - s1_cmp2[j] for j in range(3)])
            diff2 = tuple([s2_cmp1[j] - s2_cmp2[j] for j in range(3)])
            sqr_dist1 = sum([x ** 2 for x in diff1])
            sqr_dist2 = sum([x ** 2 for x in diff2])

            self.assertAlmostEqual(self.EXPECTED_ADJ_NUC_CM_DIST2, sqr_dist1)
            self.assertAlmostEqual(self.EXPECTED_ADJ_NUC_CM_DIST2, sqr_dist2)

    def test_honeycomb_design(self) -> None:
        """ A single strand on a honeycomb grid.
                  0       8
        helix 0   [-------+
                          |
        helix 1   +--------
                  |
        helix 2   +------->
        """
        helices = [sc.Helix(grid_position=(1, 1), max_offset=8), sc.Helix(grid_position=(0, 1), max_offset=8),
                   sc.Helix(grid_position=(0, 2), max_offset=8)]
        design = sc.Design(helices=helices, grid=sc.honeycomb)
        design.strand(0, 0).to(8).cross(1).move(-8).cross(2).to(8)

        # expected values for verification
        expected_num_nucleotides = 8 * 3
        expected_strand_length = 8 * 3

        dat, top = design.to_oxdna_format()
        dat = dat.strip().split('\n')
        top = top.strip().split('\n')

        # check length of output files are as expected (matches # of nucleotides plus header size)
        self.assertEqual(expected_num_nucleotides + 3, len(dat))
        self.assertEqual(expected_num_nucleotides + 1, len(top))

        # find relevant values for nucleotides
        cm_poss = []  # center of mass position
        nbrs_3p = []
        nbrs_5p = []

        for line in dat[3:]:
            data = line.strip().split()
            # make sure there are 15 values per line (3 values per vector * 5 vectors per line)
            # order of vectors: center of mass position, backbone base, normal, velocity, angular velocity
            self.assertEqual(15, len(data))

            cm_poss.append(tuple([float(x) for x in data[0:3]]))
            bb_vec = tuple([float(x) for x in data[3:6]])  # backbone base vector
            nm_vec = tuple([float(x) for x in data[6:9]])  # normal vector

            # make sure normal vectors and backbone vectors are unit length
            sqr_bb_vec = sum([x ** 2 for x in bb_vec])
            sqr_nm_vec = sum([x ** 2 for x in nm_vec])
            self.assertAlmostEqual(1.0, sqr_bb_vec)
            self.assertAlmostEqual(1.0, sqr_nm_vec)

            for value in data[9:]:  # values for velocity and angular velocity vectors are 0
                self.assertAlmostEqual(0, float(value))

        strand1_idxs = []
        for nuc_idx, line in enumerate(top[1:]):
            data = line.strip().split()
            # make sure there are 4 values per line: strand, base, 3' neighbor, 5' neighbor
            self.assertEqual(4, len(data))

            # make sure there's only 1 strand
            strand_num = int(data[0])
            self.assertEqual(1, strand_num)
            # make sure base is valid
            base = data[1]
            self.assertIn(base, ['A', 'C', 'G', 'T'])

            nbrs_3p.append(int(data[2]))
            nbrs_5p.append(int(data[3]))

            # append start of strand (no 5' neighbor) to list of indexes for strand
            neighbor_5 = int(data[3])
            if neighbor_5 == -1:
                strand1_start = nuc_idx
                strand1_idxs.append(strand1_start)

        # reconstruct strand using indices from oxDNA files
        next_idx = nbrs_3p[strand1_start]
        while next_idx >= 0:
            strand1_idxs.append(next_idx)
            next_idx = nbrs_3p[strand1_idxs[-1]]

        # assert that strand is correct length
        self.assertEqual(expected_strand_length, len(strand1_idxs))

        for i in range(expected_strand_length - 1):
            # ignore nucleotide distances between domains (on crossovers)
            if i == 7:
                continue

            if i == 15:
                continue

            strand1_nuc_idx1 = strand1_idxs[i]
            strand1_nuc_idx2 = strand1_idxs[i + 1]

            # find the center of mass for adjacent nucleotides
            s1_cmp1 = cm_poss[strand1_nuc_idx1]
            s1_cmp2 = cm_poss[strand1_nuc_idx2]

            # calculate and verify squared distance between adjacent nucleotides in a domain
            diff1 = tuple([s1_cmp1[j] - s1_cmp2[j] for j in range(3)])
            sqr_dist1 = sum([x ** 2 for x in diff1])

            self.assertAlmostEqual(self.EXPECTED_ADJ_NUC_CM_DIST2, sqr_dist1)

    def test_deletion_design(self) -> None:
        """
        0      7
        [---X-->
        """
        helix = [sc.Helix(max_offset=7)]
        design = sc.Design(helices=helix, grid=sc.square)
        design.strand(0, 0).to(7)
        design.add_deletion(helix=0, offset=4)

        # expected values for verification
        expected_num_nucleotides = 6
        expected_strand_length = 6

        dat, top = design.to_oxdna_format()
        dat = dat.strip().split('\n')
        top = top.strip().split('\n')

        # check length of output files are as expected (matches # of nucleotides plus header size)
        self.assertEqual(expected_num_nucleotides + 3, len(dat))
        self.assertEqual(expected_num_nucleotides + 1, len(top))

        # find relevant values for nucleotides
        cm_poss = []  # center of mass position
        nbrs_3p = []
        nbrs_5p = []

        for line in dat[3:]:
            data = line.strip().split()
            # make sure there are 15 values per line (3 values per vector * 5 vectors per line)
            # order of vectors: center of mass position, backbone base, normal, velocity, angular velocity
            self.assertEqual(15, len(data))

            cm_poss.append(tuple([float(x) for x in data[0:3]]))
            bb_vec = tuple([float(x) for x in data[3:6]])  # backbone base vector
            nm_vec = tuple([float(x) for x in data[6:9]])  # normal vector

            # make sure normal vectors and backbone vectors are unit length
            sqr_bb_vec = sum([x ** 2 for x in bb_vec])
            sqr_nm_vec = sum([x ** 2 for x in nm_vec])
            self.assertAlmostEqual(1.0, sqr_bb_vec)
            self.assertAlmostEqual(1.0, sqr_nm_vec)

            for value in data[9:]:  # values for velocity and angular velocity vectors are 0
                self.assertAlmostEqual(0, float(value))

        strand1_idxs = []
        for nuc_idx, line in enumerate(top[1:]):
            data = line.strip().split()
            # make sure there are 4 values per line: strand, base, 3' neighbor, 5' neighbor
            self.assertEqual(4, len(data))

            # make sure there's only 1 strand
            strand_num = int(data[0])
            self.assertEqual(1, strand_num)
            # make sure base is valid
            base = data[1]
            self.assertIn(base, ['A', 'C', 'G', 'T'])

            nbrs_3p.append(int(data[2]))
            nbrs_5p.append(int(data[3]))

            # append start of strand (no 5' neighbor) to list of indexes for strand
            neighbor_5 = int(data[3])
            if neighbor_5 == -1:
                strand1_start = nuc_idx
                strand1_idxs.append(strand1_start)

        # reconstruct strand using indices from oxDNA files
        next_idx = nbrs_3p[strand1_start]
        while next_idx >= 0:
            strand1_idxs.append(next_idx)
            next_idx = nbrs_3p[strand1_idxs[-1]]

        # assert that strand is correct length
        self.assertEqual(expected_strand_length, len(strand1_idxs))

        for i in range(expected_strand_length - 1):
            strand1_nuc_idx1 = strand1_idxs[i]
            strand1_nuc_idx2 = strand1_idxs[i + 1]

            # find the center of mass for adjacent nucleotides
            s1_cmp1 = cm_poss[strand1_nuc_idx1]
            s1_cmp2 = cm_poss[strand1_nuc_idx2]

            # calculate and verify squared distance between adjacent nucleotides in a domain
            diff1 = tuple([s1_cmp1[j] - s1_cmp2[j] for j in range(3)])
            sqr_dist1 = sum([x ** 2 for x in diff1])

            self.assertAlmostEqual(self.EXPECTED_ADJ_NUC_CM_DIST2, sqr_dist1)

    def test_insertion_design(self) -> None:
        """
        0      7
            ^ insertion of length 1
        [------>
        """
        helix = [sc.Helix(max_offset=10)]
        design = sc.Design(helices=helix, grid=sc.square)
        design.strand(0, 0).to(7)
        design.add_insertion(helix=0, offset=4, length=1)

        # expected values for verification
        expected_num_nucleotides = 8
        expected_strand_length = 8

        dat, top = design.to_oxdna_format()
        dat = dat.strip().split('\n')
        top = top.strip().split('\n')

        # check length of output files are as expected (matches # of nucleotides plus header size)
        self.assertEqual(expected_num_nucleotides + 3, len(dat))
        self.assertEqual(expected_num_nucleotides + 1, len(top))

        # find relevant values for nucleotides
        cm_poss = []  # center of mass position
        nbrs_3p = []
        nbrs_5p = []

        for line in dat[3:]:
            data = line.strip().split()
            # make sure there are 15 values per line (3 values per vector * 5 vectors per line)
            # order of vectors: center of mass position, backbone base, normal, velocity, angular velocity
            self.assertEqual(15, len(data))

            cm_poss.append(tuple([float(x) for x in data[0:3]]))
            bb_vec = tuple([float(x) for x in data[3:6]])  # backbone base vector
            nm_vec = tuple([float(x) for x in data[6:9]])  # normal vector

            # make sure normal vectors and backbone vectors are unit length
            sqr_bb_vec = sum([x ** 2 for x in bb_vec])
            sqr_nm_vec = sum([x ** 2 for x in nm_vec])
            self.assertAlmostEqual(1.0, sqr_bb_vec)
            self.assertAlmostEqual(1.0, sqr_nm_vec)

            for value in data[9:]:  # values for velocity and angular velocity vectors are 0
                self.assertAlmostEqual(0, float(value))

        strand1_idxs = []
        for nuc_idx, line in enumerate(top[1:]):
            data = line.strip().split()
            # make sure there are 4 values per line: strand, base, 3' neighbor, 5' neighbor
            self.assertEqual(4, len(data))

            # make sure there's only 1 strand
            strand_num = int(data[0])
            self.assertEqual(1, strand_num)
            # make sure base is valid
            base = data[1]
            self.assertIn(base, ['A', 'C', 'G', 'T'])

            nbrs_3p.append(int(data[2]))
            nbrs_5p.append(int(data[3]))

            # append start of strand (no 5' neighbor) to list of indexes for strand
            neighbor_5 = int(data[3])
            if neighbor_5 == -1:
                strand1_start = nuc_idx
                strand1_idxs.append(strand1_start)

        # reconstruct strand using indices from oxDNA files
        next_idx = nbrs_3p[strand1_start]
        while next_idx >= 0:
            strand1_idxs.append(next_idx)
            next_idx = nbrs_3p[strand1_idxs[-1]]

        # assert that strand is correct length
        self.assertEqual(expected_strand_length, len(strand1_idxs))

        for i in range(expected_strand_length - 1):
            strand1_nuc_idx1 = strand1_idxs[i]
            strand1_nuc_idx2 = strand1_idxs[i + 1]

            # find the center of mass for adjacent nucleotides
            s1_cmp1 = cm_poss[strand1_nuc_idx1]
            s1_cmp2 = cm_poss[strand1_nuc_idx2]

            # calculate and verify squared distance between adjacent nucleotides in a domain
            diff1 = tuple([s1_cmp1[j] - s1_cmp2[j] for j in range(3)])
            sqr_dist1 = sum([x ** 2 for x in diff1])

            self.assertAlmostEqual(self.EXPECTED_ADJ_NUC_CM_DIST2, sqr_dist1)

    def test_loopout_design(self) -> None:
        """ 2 strands, one with a loopout
            0      7
                ^ loopout at 4 of length = 4 bases
            [------>
            <------]
        """
        helix = [sc.Helix(max_offset=14)]
        design = sc.Design(helices=helix, grid=sc.square)
        design.strand(0, 0).to(4).loopout(0, 4).to(7)
        design.strand(0, 7).to(0)

        # expected values for verification
        expected_num_nucleotides = 7 * 2 + 4
        expected_strand_1_length = 7 + 4  # strand 1 has loopout of 4
        expected_strand_2_length = 7

        dat, top = design.to_oxdna_format()
        dat = dat.strip().split('\n')
        top = top.strip().split('\n')

        # check length of output files are as expected (matches # of nucleotides plus header size)
        self.assertEqual(expected_num_nucleotides + 3, len(dat))
        self.assertEqual(expected_num_nucleotides + 1, len(top))

        # find relevant values for nucleotides
        cm_poss = []  # center of mass position
        nbrs_3p = []
        nbrs_5p = []

        for line in dat[3:]:
            data = line.strip().split()
            # make sure there are 15 values per line (3 values per vector * 5 vectors per line)
            # order of vectors: center of mass position, backbone base, normal, velocity, angular velocity
            self.assertEqual(15, len(data))

            cm_poss.append(tuple([float(x) for x in data[0:3]]))
            bb_vec = tuple([float(x) for x in data[3:6]])  # backbone base vector
            nm_vec = tuple([float(x) for x in data[6:9]])  # normal vector

            # make sure normal vectors and backbone vectors are unit length
            sqr_bb_vec = sum([x ** 2 for x in bb_vec])
            sqr_nm_vec = sum([x ** 2 for x in nm_vec])
            self.assertAlmostEqual(1.0, sqr_bb_vec)
            self.assertAlmostEqual(1.0, sqr_nm_vec)

            for value in data[9:]:  # values for velocity and angular velocity vectors are 0
                self.assertAlmostEqual(0, float(value))

        strand1_idxs = []
        strand2_idxs = []
        for nuc_idx, line in enumerate(top[1:]):
            data = line.strip().split()
            # make sure there are 4 values per line: strand, base, 3' neighbor, 5' neighbor
            self.assertEqual(4, len(data))

            # make sure there are only 2 strands
            strand_num = int(data[0])
            self.assertIn(strand_num, [1, 2])
            # make sure base is valid
            base = data[1]
            self.assertIn(base, ['A', 'C', 'G', 'T'])

            nbrs_3p.append(int(data[2]))
            nbrs_5p.append(int(data[3]))

            # append start of strand (no 5' neighbor) to list of indexes for strand
            neighbor_5 = int(data[3])
            if neighbor_5 == -1:
                if strand_num == 1:
                    strand1_start = nuc_idx
                    strand1_idxs.append(strand1_start)
                else:
                    strand2_start = nuc_idx
                    strand2_idxs.append(strand2_start)

        # reconstruct strands using indices from oxDNA files
        next_idx = nbrs_3p[strand1_start]
        while next_idx >= 0:
            strand1_idxs.append(next_idx)
            next_idx = nbrs_3p[strand1_idxs[-1]]

        next_idx = nbrs_3p[strand2_start]
        while next_idx >= 0:
            strand2_idxs.append(next_idx)
            next_idx = nbrs_3p[strand2_idxs[-1]]

        # assert that strands are the correct length
        self.assertEqual(expected_strand_1_length, len(strand1_idxs))
        self.assertEqual(expected_strand_2_length, len(strand2_idxs))

        # calculate distance between nucleotides preceding and following loopout and confirm that it matches expected distance between nucleotides.
        # loopout here is nucleotides 4 to 7, so check distance between nucleotides 3 to 8
        cm_poss_pre_loopout = cm_poss[strand1_idxs[3]]
        cm_poss_post_loopout = cm_poss[strand1_idxs[8]]
        diff = tuple([cm_poss_pre_loopout[j] - cm_poss_post_loopout[j] for j in range(3)])
        sqr_dist = sum([x ** 2 for x in diff])
        self.assertAlmostEqual(self.EXPECTED_ADJ_NUC_CM_DIST2, sqr_dist)

        for i in range(expected_strand_1_length - 1):

            if i in [3, 4, 5, 6,
                     7]:  # skip nucleotide distances having to do with loopout, as these won't have regular distance between nucleotides (i = 3 denotes distance from 3 to 4, which includes loopout)
                continue

            strand1_nuc_idx1 = strand1_idxs[i]
            strand1_nuc_idx2 = strand1_idxs[i + 1]

            # find the center of mass for adjacent nucleotides
            s1_cmp1 = cm_poss[strand1_nuc_idx1]
            s1_cmp2 = cm_poss[strand1_nuc_idx2]

            # calculate and verify squared distance between adjacent nucleotides in a domain
            diff1 = tuple([s1_cmp1[j] - s1_cmp2[j] for j in range(3)])
            sqr_dist1 = sum([x ** 2 for x in diff1])

            self.assertAlmostEqual(self.EXPECTED_ADJ_NUC_CM_DIST2, sqr_dist1)

        for i in range(expected_strand_2_length - 1):  # check adjacent nucleotide distances for second strand
            strand2_nuc_idx1 = strand2_idxs[i]
            strand2_nuc_idx2 = strand2_idxs[i + 1]

            # find the center of mass for adjacent nucleotides
            s2_cmp1 = cm_poss[strand2_nuc_idx1]
            s2_cmp2 = cm_poss[strand2_nuc_idx2]

            # calculate and verify squared distance between adjacent nucleotides in a domain
            diff2 = tuple([s2_cmp1[j] - s2_cmp2[j] for j in range(3)])
            sqr_dist2 = sum([x ** 2 for x in diff2])

            self.assertAlmostEqual(self.EXPECTED_ADJ_NUC_CM_DIST2, sqr_dist2)
