# import os
# import sys
# sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import scadnano as sc

import unittest

#TODO: add tests for mutation methods on DNADesign

#TODO: mutator methods let me create this strand (which should be illegal); add a test for it
# {
#   "color": {"r": 51, "g": 51, "b": 51},
#   "substrands": [
#     {"helix_idx": 2, "right": false, "start": 40, "end": 48},
#     {"helix_idx": 2, "right": false, "start": 32, "end": 48, "deletions": [44]},
#     {"helix_idx": 3, "right": true, "start": 32, "end": 40}
#   ]
# }

class TestIllegalStructuresPrevented(unittest.TestCase):

    def test_overlapping_caught_in_strange_counterexample(self):
        # found this counterexample as a simplified version of something caught in practice
        helices = [sc.Helix(0, 25), sc.Helix(1, 25)]

        s1_left_ss0 = sc.Substrand(0, sc.left, 0, 5)
        s1_ss1 = sc.Substrand(0, sc.right, 0, 15)
        s1_right_ss0 = sc.Substrand(0, sc.left, 5, 15)
        s1 = sc.Strand([s1_left_ss0, s1_ss1, s1_right_ss0])

        s2_ss1 = sc.Substrand(0, sc.right, 10, 20)
        s2_ss0 = sc.Substrand(0, sc.left, 10, 20)
        s2 = sc.Strand([s2_ss1, s2_ss0])

        strands = [s1, s2]

        with self.assertRaises(sc.IllegalDNADesignError):
            design = sc.DNADesign(helices=helices, strands=strands, grid=sc.square)

    def test_major_tick_outside_range(self):
        with self.assertRaises(sc.IllegalDNADesignError):
            helix = sc.Helix(idx=0, max_bases=9, major_ticks=[2, 5, 10])

    def test_major_tick_just_inside_range(self):
        helix = sc.Helix(idx=0, max_bases=9, major_ticks=[0, 5, 9])

    def test_two_illegally_overlapping_strands(self):
        helix = sc.Helix(idx=0, max_bases=9)
        ss_bot = sc.Substrand(helix_idx=0, right=sc.left, start=0, end=9)
        ss_top = sc.Substrand(helix_idx=0, right=sc.left, start=0, end=9)
        strand_bot = sc.Strand(substrands=[ss_bot])
        strand_top = sc.Strand(substrands=[ss_top])
        strands = [strand_bot, strand_top]
        with self.assertRaises(sc.IllegalDNADesignError):
            sc.DNADesign(grid=sc.square, helices=[helix], strands=strands)

    def test_two_nonconsecutive_illegally_overlapping_strands(self):
        helix = sc.Helix(idx=0, max_bases=9)
        ss_top1 = sc.Substrand(helix_idx=0, right=sc.left, start=0, end=5)
        ss_bot = sc.Substrand(helix_idx=0, right=sc.right, start=2, end=9)
        ss_top2 = sc.Substrand(helix_idx=0, right=sc.left, start=4, end=8)
        strand_bot = sc.Strand(substrands=[ss_bot])
        strand_top1 = sc.Strand(substrands=[ss_top1])
        strand_top2 = sc.Strand(substrands=[ss_top2])
        strands = [strand_bot, strand_top1, strand_top2]
        with self.assertRaises(sc.IllegalDNADesignError):
            sc.DNADesign(grid=sc.square, helices=[helix], strands=strands)

    def test_four_legally_leapfrogging_strands(self):
        helix = sc.Helix(idx=0, max_bases=9)
        ss_top1 = sc.Substrand(helix_idx=0, right=sc.left, start=0, end=20)
        ss_bot1 = sc.Substrand(helix_idx=0, right=sc.right, start=10, end=30)
        ss_top2 = sc.Substrand(helix_idx=0, right=sc.left, start=20, end=40)
        ss_bot2 = sc.Substrand(helix_idx=0, right=sc.right, start=30, end=50)
        strand_bot1 = sc.Strand(substrands=[ss_bot1])
        strand_bot2 = sc.Strand(substrands=[ss_bot2])
        strand_top1 = sc.Strand(substrands=[ss_top1])
        strand_top2 = sc.Strand(substrands=[ss_top2])
        strands = [strand_bot1, strand_bot2, strand_top1, strand_top2]
        sc.DNADesign(grid=sc.square, helices=[helix], strands=strands)

    def test_helices_skip_index(self):
        h1 = sc.Helix(idx=0, max_bases=9)
        h2 = sc.Helix(idx=2, max_bases=9)
        with self.assertRaises(sc.IllegalDNADesignError):
            sc.DNADesign(grid=sc.square, helices=[h1, h2], strands=[])

    def test_helices_repeat_index(self):
        h1 = sc.Helix(idx=0, max_bases=9)
        h2 = sc.Helix(idx=1, max_bases=9)
        h3 = sc.Helix(idx=0, max_bases=9)
        with self.assertRaises(sc.IllegalDNADesignError):
            sc.DNADesign(grid=sc.square, helices=[h1, h2, h3], strands=[])

    def test_strand_references_nonexistent_helix(self):
        h1 = sc.Helix(idx=0, max_bases=9)
        h2 = sc.Helix(idx=1, max_bases=9)
        ss_bot = sc.Substrand(helix_idx=2, right=sc.left, start=0, end=9)
        ss_top = sc.Substrand(helix_idx=3, right=sc.left, start=0, end=9)
        strand_bot = sc.Strand(substrands=[ss_bot])
        strand_top = sc.Strand(substrands=[ss_top])
        strands = [strand_bot, strand_top]
        with self.assertRaises(sc.IllegalDNADesignError):
            sc.DNADesign(grid=sc.square, helices=[h1, h2], strands=strands)


class TestAssignDNA(unittest.TestCase):

    def test_assign_dna__two_equal_length_strands_on_one_helix(self):
        # <---]
        # CAAAA
        # GTTTT
        # [--->
        helix = sc.Helix(idx=0, max_bases=5)
        ss_bot = sc.Substrand(helix_idx=0, right=sc.right, start=0, end=5)
        ss_top = sc.Substrand(helix_idx=0, right=sc.left, start=0, end=5)
        strand_bot = sc.Strand(substrands=[ss_bot])
        strand_top = sc.Strand(substrands=[ss_top])
        strands = [strand_bot, strand_top]
        design = sc.DNADesign(grid=sc.square, helices=[helix], strands=strands)
        design.assign_dna(strand_top, 'AAAAC')
        self.assertEqual('GTTTT', strand_bot.dna_sequence)

    def test_assign_dna__one_strand_assigned_by_complement_from_two_other_strands(self):
        #   0123     4567
        # <-AAAC-] <-GGGA-]
        # [-TTTG-----CCCT->
        helix=sc.Helix(idx=0, max_bases=8)
        ss_top_left = sc.Substrand(0, sc.left, 0, 4)
        ss_top_right = sc.Substrand(0, sc.left, 4, 8)
        ss_bot = sc.Substrand(0,sc.right, 0, 8)
        st_top_left = sc.Strand([ss_top_left])
        st_top_right = sc.Strand([ss_top_right])
        st_bot = sc.Strand([ss_bot])
        design = sc.DNADesign(helices=[helix], strands=[st_bot, st_top_left, st_top_right])
        design.assign_dna(st_top_left, 'CAAA')
        self.assertEqual('TTTG????', st_bot.dna_sequence)
        design.assign_dna(st_top_right, 'AGGG')
        self.assertEqual('TTTGCCCT', st_bot.dna_sequence)

    def test_assign_dna__adapter_assigned_from_scaffold_and_tiles(self):
        #XXX: it appears the behavior this tests (which the other tests miss) is assigning DNA to
        # tile0 first, then to tile1, and adap is connected to each of them on different helices.
        # This means that when tile1 is assigned, we need to ensure when assigning to adap that we
        # keep the old information and don't discard it by simply padding the shorter portion of it on
        # helix 1 with ?'s, but remember the old DNA sequence.

        #        01 2345     6789  01
        # adap    [-TTTC-----CATT-------+
        # scaf <-GT-AAAG-+ <-GTAA--AA-] |
        #                |              |
        #      [-AA-TTTG-+ [-TGCC--GG-> |
        #         <-AAAC-----ACGG-------+
        h0=sc.Helix(idx=0, max_bases=12)
        h1=sc.Helix(idx=1, max_bases=12)
        scaf0_ss = sc.Substrand(0, sc.left, 0, 6)
        scaf1_ss = sc.Substrand(1, sc.right, 0, 6)
        tile1_ss = sc.Substrand(1, sc.right, 6, 12)
        tile0_ss = sc.Substrand(0, sc.left, 6, 12)
        adap0_ss = sc.Substrand(0, sc.right, 2, 10)
        adap1_ss = sc.Substrand(1, sc.left, 2, 10)
        scaf = sc.Strand([scaf1_ss, scaf0_ss])
        adap = sc.Strand([adap0_ss, adap1_ss])
        tile0 = sc.Strand([tile0_ss])
        tile1 = sc.Strand([tile1_ss])

        design = sc.DNADesign(helices=[h0, h1], strands=[scaf, adap, tile0, tile1])

        design.assign_dna(tile0, 'AA AATG')
        self.assertEqual('???? CATT ???? ????'.replace(' ',''), adap.dna_sequence)

        design.assign_dna(tile1, 'TGCC GG')
        self.assertEqual('???? CATT GGCA ????'.replace(' ',''), adap.dna_sequence)

        design.assign_dna(scaf, 'AA TTTG GAAA TG')
        self.assertEqual('TTTC CATT GGCA CAAA'.replace(' ',''), adap.dna_sequence)

    def test_assign_dna__adapter_assigned_from_scaffold_and_tiles_with_deletions(self):
        #XXX: it appears the behavior this tests (which the other tests miss) is assigning DNA to
        # tile0 first, then to tile1, and adap is connected to each of them on different helices.
        # This means that when tile1 is assigned, we need to ensure when assigning to adap that we
        # keep the old information and don't discard it by simply padding the shorter portion of it on
        # helix 1 with ?'s, but remember the old DNA sequence.

        #        01 2345     6789  01
        #            X         X           deletions
        # adap    [-T TC-----CA T-------+
        # scaf <-GT-A AG-+ <-GT A--AA-] |
        #                |              |
        #      [-AA-TTTG-+ [-TG C--GG-> |
        #         <-AAAC-----AC G-------+
        #                      X           deletions
        h0=sc.Helix(idx=0, max_bases=12)
        h1=sc.Helix(idx=1, max_bases=12)
        scaf0_ss = sc.Substrand(0, sc.left, 0, 6)
        scaf1_ss = sc.Substrand(1, sc.right, 0, 6)
        tile1_ss = sc.Substrand(1, sc.right, 6, 12)
        tile0_ss = sc.Substrand(0, sc.left, 6, 12)
        adap0_ss = sc.Substrand(0, sc.right, 2, 10)
        adap1_ss = sc.Substrand(1, sc.left, 2, 10)
        scaf = sc.Strand([scaf1_ss, scaf0_ss])
        adap = sc.Strand([adap0_ss, adap1_ss])
        tile0 = sc.Strand([tile0_ss])
        tile1 = sc.Strand([tile1_ss])

        design = sc.DNADesign(helices=[h0, h1], strands=[scaf, adap, tile0, tile1])
        design.add_deletion(0, 3)
        design.add_deletion(0, 8)
        design.add_deletion(1, 8)

        design.assign_dna(tile0, 'AA ATG')
        self.assertEqual('??? CAT ??? ????'.replace(' ',''), adap.dna_sequence)

        design.assign_dna(tile1, 'TGC GG')
        self.assertEqual('??? CAT GCA ????'.replace(' ',''), adap.dna_sequence)

        design.assign_dna(scaf, 'AA TTTG GAA TG')
        self.assertEqual('TTC CAT GCA CAAA'.replace(' ',''), adap.dna_sequence)

    def test_assign_dna__adapter_assigned_from_scaffold_and_tiles_with_insertions(self):
        #XXX: it appears the behavior this tests (which the other tests miss) is assigning DNA to
        # tile0 first, then to tile1, and adap is connected to each of them on different helices.
        # This means that when tile1 is assigned, we need to ensure when assigning to adap that we
        # keep the old information and don't discard it by simply padding the shorter portion of it on
        # helix 1 with ?'s, but remember the old DNA sequence.

        #        01 2345     678I9  01
        #                       I          insertions
        # adap    [-TTTC-----CATTT-------+
        # scaf <-GT-AAAG-+ <-GTAAA--AA-] |
        #                |              |
        #      [-AA-TTTG-+ [-TGCCC--GG-> |
        #         <-AAAC-----ACGGG-------+
        #                       I          insertions
        h0=sc.Helix(idx=0, max_bases=12)
        h1=sc.Helix(idx=1, max_bases=12)
        scaf0_ss = sc.Substrand(0, sc.left, 0, 6)
        scaf1_ss = sc.Substrand(1, sc.right, 0, 6)
        tile1_ss = sc.Substrand(1, sc.right, 6, 12)
        tile0_ss = sc.Substrand(0, sc.left, 6, 12)
        adap0_ss = sc.Substrand(0, sc.right, 2, 10)
        adap1_ss = sc.Substrand(1, sc.left, 2, 10)
        scaf = sc.Strand([scaf1_ss, scaf0_ss])
        adap = sc.Strand([adap0_ss, adap1_ss])
        tile0 = sc.Strand([tile0_ss])
        tile1 = sc.Strand([tile1_ss])

        design = sc.DNADesign(helices=[h0, h1], strands=[scaf, adap, tile0, tile1])
        design.add_insertion(0, 8, 1)
        design.add_insertion(1, 8, 1)

        design.assign_dna(tile0, 'AA AAATG')
        self.assertEqual('???? CATTT ????? ????'.replace(' ', ''), adap.dna_sequence)

        design.assign_dna(tile1, 'TGCCC GG')
        self.assertEqual('???? CATTT GGGCA ????'.replace(' ', ''), adap.dna_sequence)

        design.assign_dna(scaf, 'AA TTTG GAAA TG')
        self.assertEqual('TTTC CATTT GGGCA CAAA'.replace(' ', ''), adap.dna_sequence)

    def test_assign_dna__dna_sequence_shorter_than_complementary_strand_right_strand_longer(self):
        # <---]
        # CAAAA
        # GTTTT?????
        # [-------->
        helix = sc.Helix(idx=0, max_bases=10)
        ss_long = sc.Substrand(helix_idx=0, right=sc.right, start=0, end=10)
        ss_short = sc.Substrand(helix_idx=0, right=sc.left, start=0, end=5)
        strand_long = sc.Strand(substrands=[ss_long])
        strand_short = sc.Strand(substrands=[ss_short])
        strands = [strand_long, strand_short]
        design = sc.DNADesign(grid=sc.square, helices=[helix], strands=strands)
        design.assign_dna(strand_short, 'AAAAC')
        self.assertEqual('GTTTT?????', strand_long.dna_sequence)

    def test_assign_dna__dna_sequence_shorter_than_complementary_strand_left_strand_longer(self):
        # [--->
        # AAAAC
        # TTTTG?????
        # <--------]
        helix = sc.Helix(idx=0, max_bases=10)
        ss_long = sc.Substrand(helix_idx=0, right=sc.left, start=0, end=10)
        ss_short = sc.Substrand(helix_idx=0, right=sc.right, start=0, end=5)
        strand_long = sc.Strand(substrands=[ss_long])
        strand_short = sc.Strand(substrands=[ss_short])
        strands = [strand_long, strand_short]
        design = sc.DNADesign(grid=sc.square, helices=[helix], strands=strands)
        design.assign_dna(strand_short, 'AAAAC')
        self.assertEqual('?????GTTTT', strand_long.dna_sequence)

    def test_assign_dna__dna_sequence_with_uncomplemented_substrand_on_different_helix(self):
        # <---]
        # CAAAA
        # GTTTT?????
        # [--------+
        #          |
        #        <-+
        #        ???
        helix0 = sc.Helix(idx=0, max_bases=10)
        helix1 = sc.Helix(idx=1, max_bases=10)
        ss_long = sc.Substrand(helix_idx=0, right=sc.right, start=0, end=10)
        ss_long_h1 = sc.Substrand(helix_idx=0, right=sc.left, start=7, end=10)
        ss_short = sc.Substrand(helix_idx=0, right=sc.left, start=0, end=5)
        strand_long = sc.Strand(substrands=[ss_long, ss_long_h1])
        strand_short = sc.Strand(substrands=[ss_short])
        strands = [strand_long, strand_short]
        design = sc.DNADesign(grid=sc.square, helices=[helix0,helix1], strands=strands)
        design.assign_dna(strand_short, 'AAAAC')
        self.assertEqual('GTTTT????????', strand_long.dna_sequence)

    def test_assign_dna__dna_sequence_with_uncomplemented_substrand_on_different_helix_wildcards_both_ends(self):
        #      <---]
        #      CAAAA
        # ?????GTTTT
        # [--------+
        #          |
        #        <-+
        #        ???
        helix0 = sc.Helix(idx=0, max_bases=10)
        helix1 = sc.Helix(idx=1, max_bases=10)
        ss_long_h0 = sc.Substrand(helix_idx=0, right=sc.right, start=0, end=10)
        ss_long_h1 = sc.Substrand(helix_idx=1, right=sc.left, start=7, end=10)
        ss_short_h0 = sc.Substrand(helix_idx=0, right=sc.left, start=5, end=10)
        strand_long = sc.Strand(substrands=[ss_long_h0, ss_long_h1])
        strand_short = sc.Strand(substrands=[ss_short_h0])
        strands = [strand_long, strand_short]
        design = sc.DNADesign(grid=sc.square, helices=[helix0,helix1], strands=strands)
        design.assign_dna(strand_short, 'AAAAC')
        self.assertEqual('?????GTTTT???', strand_long.dna_sequence)

    def test_assign_dna__one_helix_with_one_bottom_strand_and_three_top_strands(self):
        #  012   345   678
        # -TTT> -GGG> -CCC>
        # <AAA---CCC---GGG-
        #  876   543   210
        helix = sc.Helix(idx=0, max_bases=9)
        ss_bot = sc.Substrand(helix_idx=0, right=sc.left, start=0, end=9)
        ss_top1 = sc.Substrand(helix_idx=0, right=sc.right, start=0, end=3)
        ss_top2 = sc.Substrand(helix_idx=0, right=sc.right, start=3, end=6)
        ss_top3 = sc.Substrand(helix_idx=0, right=sc.right, start=6, end=9)
        strand_bot = sc.Strand(substrands=[ss_bot])
        strand_top1 = sc.Strand(substrands=[ss_top1])
        strand_top2 = sc.Strand(substrands=[ss_top2])
        strand_top3 = sc.Strand(substrands=[ss_top3])
        strands = [strand_bot, strand_top1, strand_top2, strand_top3]
        design = sc.DNADesign(grid=sc.square, helices=[helix], strands=strands)
        design.assign_dna(strand_bot, 'AAACCCGGG')
        self.assertEqual('CCC', strand_top1.dna_sequence)
        self.assertEqual('GGG', strand_top2.dna_sequence)
        self.assertEqual('TTT', strand_top3.dna_sequence)

    def test_assign_dna__two_helices_with_multiple_substrand_intersections(self):
        #         012    345   678    901
        # M13   [-ACC----TAA---GAA----AAC---+
        #       +-TGG-]<-ATT-+ CTT----TTG-+ |
        #       |            | |          | |
        #       +-GAT----TTC-+ ATG->[-AGT-+ |
        #       <-CTA----AAG---TAC----TCA---+
        h0 = sc.Helix(idx=0, max_bases=12)
        h1 = sc.Helix(idx=1, max_bases=12)

        scaf0_ss = sc.Substrand(helix_idx=0, right=sc.right, start=0, end=12)
        scaf1_ss = sc.Substrand(helix_idx=1, right=sc.left, start=0, end=12)
        scaf = sc.Strand(substrands=[scaf0_ss, scaf1_ss])

        first_stap0_left_ss = sc.Substrand(helix_idx=0, right=sc.left, start=0, end=3)
        first_stap1_ss = sc.Substrand(helix_idx=1, right=sc.right, start=0, end=6)
        first_stap0_right_ss = sc.Substrand(helix_idx=0, right=sc.left, start=3, end=6)
        first_stap = sc.Strand(substrands=[first_stap0_left_ss, first_stap1_ss, first_stap0_right_ss])

        second_stap1_right_ss = sc.Substrand(helix_idx=1, right=sc.right, start=9, end=12)
        second_stap0_ss = sc.Substrand(helix_idx=0, right=sc.left, start=6, end=12)
        second_stap1_left_ss = sc.Substrand(helix_idx=1, right=sc.right, start=6, end=9)
        second_stap = sc.Strand(substrands=[second_stap1_right_ss, second_stap0_ss, second_stap1_left_ss])

        strands = [scaf, first_stap, second_stap]
        design = sc.DNADesign(grid=sc.square, helices=[h0, h1], strands=strands)
        design.assign_dna(scaf, 'ACC TAA GAA AAC ACT CAT GAA ATC'.replace(' ', ''))
        self.assertEqual('GGT GAT TTC TTA'.replace(' ', ''), first_stap.dna_sequence)
        self.assertEqual('AGT GTT TTC ATG'.replace(' ', ''), second_stap.dna_sequence)

    def test_assign_dna__upper_left_edge_staple_of_16H_origami_rectangle(self):
        # staple <ACATAAGAAAACGGAG--+
        # M13   +-TGTATTCTTTTGCCTC> |
        #       |                   |
        #       +-GATTTTGTGAGTAGAA- |
        #        -CTAAAACACTCATCTT--+
        h0 = sc.Helix(idx=0, max_bases=16)
        h1 = sc.Helix(idx=1, max_bases=16)
        scaf0_ss = sc.Substrand(helix_idx=0, right=sc.right, start=0, end=16)
        scaf1_ss = sc.Substrand(helix_idx=1, right=sc.left, start=0, end=16)
        stap0_ss = sc.Substrand(helix_idx=0, right=sc.left, start=0, end=16)
        stap1_ss = sc.Substrand(helix_idx=1, right=sc.right, start=0, end=16)
        scaf = sc.Strand(substrands=[scaf1_ss, scaf0_ss])
        stap = sc.Strand(substrands=[stap1_ss, stap0_ss])
        strands = [scaf, stap]
        design = sc.DNADesign(grid=sc.square, helices=[h0, h1], strands=strands)

        seq_m13_upper_left = 'AAGATGAGTGTTTTAGTGTATTCTTTTGCCTC'
        design.assign_dna(scaf, seq_m13_upper_left)
        expected_seq_stap_upperleft = 'CTAAAACACTCATCTTGAGGCAAAAGAATACA'
        self.assertEqual(expected_seq_stap_upperleft, stap.dna_sequence)

    def test_assign_dna__2helix_with_deletions(self):
        # scaf index: 2     3  4     5
        # offset:     0 D1  2  3 D4  5
        #             +     -  -     +
        #            /C     A  T     C\
        #           | G     T  A     G |
        # helix 0   | <     +  +     ] |
        #           |       |  |       |
        # helix 1   | [     +  +     > |
        #           | T     T  A     C |
        #            \A     A  T     G/
        #             +     ]  <     +
        # offset:     0 D1  2  3 D4  5
        # scaf index: 1     0  7     6
        width = 6
        width_h = width // 2
        helices = [sc.Helix(0, width), sc.Helix(1, width)]
        stap_left_ss1 = sc.Substrand(1, sc.right, 0, width_h)
        stap_left_ss0 = sc.Substrand(0, sc.left, 0, width_h)
        stap_right_ss0 = sc.Substrand(0, sc.left, width_h, width)
        stap_right_ss1 = sc.Substrand(1, sc.right, width_h, width)
        scaf_ss1_left = sc.Substrand(1, sc.left, 0, width_h)
        scaf_ss0 = sc.Substrand(0, sc.right, 0, width)
        scaf_ss1_right = sc.Substrand(1, sc.left, width_h, width)
        stap_left = sc.Strand([stap_left_ss1, stap_left_ss0])
        stap_right = sc.Strand([stap_right_ss0, stap_right_ss1])
        scaf = sc.Strand([scaf_ss1_left, scaf_ss0, scaf_ss1_right], color=sc.default_scaffold_color)
        strands = [stap_left, stap_right, scaf]
        design = sc.DNADesign(helices=helices, strands=strands, grid=sc.square)
        design.add_deletion(helix_idx=0, offset=1)
        design.add_deletion(helix_idx=0, offset=4)
        design.add_deletion(helix_idx=1, offset=1)
        design.add_deletion(helix_idx=1, offset=4)
        design.assign_dna(scaf, 'AACATCGT')
        self.assertEqual("AACATCGT", scaf.dna_sequence)
        self.assertEqual("TTTG", stap_left.dna_sequence)
        self.assertEqual("GAAC", stap_right.dna_sequence)

    def test_assign_dna__wildcards_simple(self):
        #  012   345   678
        # -TTC> -GGA> -CCT>
        # <AAG---CCT---GGA-
        #  876   543   210
        helix = sc.Helix(idx=0, max_bases=9)
        ss_bot = sc.Substrand(helix_idx=0, right=sc.left, start=0, end=9)
        ss_top1 = sc.Substrand(helix_idx=0, right=sc.right, start=0, end=3)
        ss_top2 = sc.Substrand(helix_idx=0, right=sc.right, start=3, end=6)
        ss_top3 = sc.Substrand(helix_idx=0, right=sc.right, start=6, end=9)
        strand_bot = sc.Strand(substrands=[ss_bot])
        strand_top1 = sc.Strand(substrands=[ss_top1])
        strand_top2 = sc.Strand(substrands=[ss_top2])
        strand_top3 = sc.Strand(substrands=[ss_top3])
        strands = [strand_bot, strand_top1, strand_top2, strand_top3]
        design = sc.DNADesign(grid=sc.square, helices=[helix], strands=strands)

        design.assign_dna(strand_top1, 'TTC')
        self.assertEqual('??????GAA', strand_bot.dna_sequence)

        design.assign_dna(strand_top3, 'CCT')
        self.assertEqual('AGG???GAA', strand_bot.dna_sequence)

        design.assign_dna(strand_top2, 'GGA')
        self.assertEqual('AGGTCCGAA', strand_bot.dna_sequence)


    def test_assign_dna__wildcards_multiple_overlaps(self):
        #  012   345   678   901   234   567
        #       +---------------+
        #       |               |
        #       |         +-----|-------+
        #       |         |     |       |
        # -ACG> +TTC> -GGA+ -AAC+ -TGC> +TTG>
        # <TGC---AAG---CCT---TTG---ACG---AAC---???-
        #  098   765   432   109   876   543   210
        helix = sc.Helix(idx=0, max_bases=9)
        ss_bot = sc.Substrand(helix_idx=0, right=sc.left, start=0, end=21)
        ss_top0 = sc.Substrand(helix_idx=0, right=sc.right, start=0, end=3)
        ss_top3 = sc.Substrand(helix_idx=0, right=sc.right, start=3, end=6)
        ss_top6 = sc.Substrand(helix_idx=0, right=sc.right, start=6, end=9)
        ss_top9 = sc.Substrand(helix_idx=0, right=sc.right, start=9, end=12)
        ss_top12 = sc.Substrand(helix_idx=0, right=sc.right, start=12, end=15)
        ss_top15 = sc.Substrand(helix_idx=0, right=sc.right, start=15, end=18)
        strand_bot = sc.Strand(substrands=[ss_bot])
        strand_top_small0 = sc.Strand(substrands=[ss_top0])
        strand_top_small12 = sc.Strand(substrands=[ss_top12])
        strand_top_big9 = sc.Strand(substrands=[ss_top9,ss_top3])
        strand_top_big6 = sc.Strand(substrands=[ss_top6,ss_top15])
        strands = [strand_bot, strand_top_small0, strand_top_small12, strand_top_big9, strand_top_big6]
        design = sc.DNADesign(grid=sc.square, helices=[helix], strands=strands)

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
        ss0 = sc.Substrand(0, sc.right, 0, 10)
        ss1 = sc.Substrand(1, sc.left, 0, 10)
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
        ss0 = sc.Substrand(0, sc.right, 0, 10, deletions=[2, 5, 6])
        ss1 = sc.Substrand(1, sc.left, 0, 10, deletions=[2, 6, 7])
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
        ss0 = sc.Substrand(0, sc.right, 0, 10, insertions=[(2, 1), (6, 2)])
        ss1 = sc.Substrand(1, sc.left, 0, 10, insertions=[(2, 1), (6, 2)])
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
        ss0 = sc.Substrand(0, sc.right, 0, 10, deletions=[4], insertions=[(2, 1), (6, 2)])
        ss1 = sc.Substrand(1, sc.left, 0, 10, deletions=[4], insertions=[(2, 1), (6, 2)])
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
        self.assertEqual("A",            ss0.dna_sequence_in(0, 0))
        self.assertEqual("AA",           ss0.dna_sequence_in(0, 1))
        self.assertEqual("AAAA",         ss0.dna_sequence_in(0, 2))
        self.assertEqual("AAAAC",        ss0.dna_sequence_in(0, 3))
        self.assertEqual("AAAAC",        ss0.dna_sequence_in(0, 4))
        self.assertEqual("AAAACC",       ss0.dna_sequence_in(0, 5))
        self.assertEqual("AAAACCCCG",    ss0.dna_sequence_in(0, 6))
        self.assertEqual("AAAACCCCGG",   ss0.dna_sequence_in(0, 7))
        self.assertEqual("AAAACCCCGGG",  ss0.dna_sequence_in(0, 8))
        self.assertEqual("AAAACCCCGGGG", ss0.dna_sequence_in(0, 9))
        self.assertEqual( "AAACCCCGGGG", ss0.dna_sequence_in(1, 9))
        self.assertEqual(  "AACCCCGGGG", ss0.dna_sequence_in(2, 9))
        self.assertEqual(    "CCCCGGGG", ss0.dna_sequence_in(3, 9))
        self.assertEqual(     "CCCGGGG", ss0.dna_sequence_in(4, 9))
        self.assertEqual(     "CCCGGGG", ss0.dna_sequence_in(5, 9))
        self.assertEqual(      "CCGGGG", ss0.dna_sequence_in(6, 9))
        self.assertEqual(         "GGG", ss0.dna_sequence_in(7, 9))
        self.assertEqual(          "GG", ss0.dna_sequence_in(8, 9))
        self.assertEqual(           "G", ss0.dna_sequence_in(9, 9))
        #
        self.assertEqual("T",            ss1.dna_sequence_in(9, 9))
        self.assertEqual("TT",           ss1.dna_sequence_in(8, 9))
        self.assertEqual("TTT",          ss1.dna_sequence_in(7, 9))
        self.assertEqual("TTTTAC",       ss1.dna_sequence_in(6, 9))
        self.assertEqual("TTTTACG",      ss1.dna_sequence_in(5, 9))
        self.assertEqual("TTTTACG",      ss1.dna_sequence_in(4, 9))
        self.assertEqual("TTTTACGT",     ss1.dna_sequence_in(3, 9))
        self.assertEqual("TTTTACGTAC",   ss1.dna_sequence_in(2, 9))
        self.assertEqual("TTTTACGTACG",  ss1.dna_sequence_in(1, 9))
        self.assertEqual("TTTTACGTACGT", ss1.dna_sequence_in(0, 9))
        self.assertEqual( "TTTACGTACGT", ss1.dna_sequence_in(0, 8))
        self.assertEqual(  "TTACGTACGT", ss1.dna_sequence_in(0, 7))
        self.assertEqual(   "TACGTACGT", ss1.dna_sequence_in(0, 6))
        self.assertEqual(      "GTACGT", ss1.dna_sequence_in(0, 5))
        self.assertEqual(       "TACGT", ss1.dna_sequence_in(0, 4))
        self.assertEqual(       "TACGT", ss1.dna_sequence_in(0, 3))
        self.assertEqual(        "ACGT", ss1.dna_sequence_in(0, 2))
        self.assertEqual(          "GT", ss1.dna_sequence_in(0, 1))
        self.assertEqual(           "T", ss1.dna_sequence_in(0, 0))

        # if TEST_OFFSETS_AT_DELETION_INSERTIONS:
        #     self.assertEqual("AAAA", ss0.dna_sequence_in(0, 3))
        #     self.assertEqual("AAAAC", ss0.dna_sequence_in(0, 5))
        #     self.assertEqual("AAAACCCCGG", ss0.dna_sequence_in(0, 7))
        #     self.assertEqual("TTTACG", ss1.dna_sequence_in(6, 10))
        #     self.assertEqual("TTTTACG", ss1.dna_sequence_in(4, 10))
        #     self.assertEqual("TTTACGTACGT", ss1.dna_sequence_in(2, 10))
