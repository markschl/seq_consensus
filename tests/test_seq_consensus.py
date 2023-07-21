import unittest
import numpy as np

from seq_consensus.consensus import AlignmentFrequencies, _convert_endgaps, \
    consensus


class TestAlignmentFuncs(unittest.TestCase):
    def test_endgaps(self):
        data = [
            ('---A', '...A'),
            ('A.--', 'A...'),
            ('A', 'A'),
            ('-ATGC--', '.ATGC..')
        ]
        for seq, converted in data:
            seq = np.array(bytearray(seq, 'utf-8'))
            _convert_endgaps(seq, ord('-'), ord('.'))
            seq = bytearray(seq).decode('utf-8')
            self.assertEqual(seq, converted)

    def test_letterfreqs(self):
        seqs = [
            '-ATTGC',
            '-AT-CC',
            '-RT-C-'
        ]
        # no special treatment for end gaps
        exp_freqs = [
            {'-': 3},
            {'A': 2, 'R': 1},
            {'T': 3},
            {'-': 2, 'T': 1},
            {'C': 2, 'G': 1},
            {'-': 1, 'C': 2}
        ]
        af = AlignmentFrequencies(sequences=seqs, free_endgaps=False)
        freqs = [dict(f) for f in af.column_freqs()]
        self.assertEqual(freqs, exp_freqs)
        # ignore end gaps
        exp_freqs = [
            {},
            {'A': 2, 'R': 1},
            {'T': 3},
            {'-': 2, 'T': 1},
            {'C': 2, 'G': 1},
            {'C': 2}
        ]
        af = AlignmentFrequencies(sequences=seqs, free_endgaps=True)
        freqs = [dict(f) for f in af.column_freqs()]
        self.assertEqual(freqs, exp_freqs)

    # def test_replace_ambigs(self):
    #     freqs = [
    #         ({'A': 3}, {'A': 3}),
    #         ({'A': 2, 'R': 1}, {'A': 2.5, 'G': 0.5}),
    #         ({'N': 1}, {'A': 0.25, 'C': 0.25, 'G': 0.25, 'T': 0.25})
    #     ]
    #     for f, exp in freqs:
    #         self.assertEqual(replace_ambig(f, IUPAC_map), exp)

    # def test_consensus_letter(self):
    #     thresholds = [0, 0.5, 0.75, 1]
    #     data = [
    #         ({'A': 1}, ('A', 'A', 'A', 'A')),
    #         ({'G': 1, 'R': 1}, ('G', 'G', 'G', 'R')),
    #         ({'G': 1, 'R': 8}, ('G', 'G', 'R', 'R')),
    #         ({'N': 1, 'A': 1}, ('A', 'A', 'N', 'N')),
    #         ({'A': 3, '-': 1}, ('A', 'A', 'A', '?'))
    #     ]
    #     amap = invert_ambig_map(IUPAC_map)
    #     for f, exp in data:
    #         for t, exp_char in zip(thresholds, exp):
    #             f = replace_ambig(f, IUPAC_map, other_valid=['-'])
    #             self.assertEqual(consensus_letter(f, amap, t), exp_char)

    def test_consensus(self):
        seqs = [
            'GTN-CTA--',
            'GC--TCAA-',
            '-CT-CTAA-',
            'ACCACYAA-'
        ]
        exp_cons = [
            (0,    'GCY-CTAA.', 'GCY-CTAA-', 'GCYCTAA'),
            (0.5,  'GCY-CTAA.', 'GCY-CTAA-', 'GCYCTAA'),
            (0.75, 'RC?-CYAA.', '?C?-CYAA-', 'RC?CYAA'),
            (1,    'RY??YYAA.', '?Y??YYA?-', 'RY??YYAA'),
        ]
        for threshold, free_cons, cons, strip_cons in exp_cons:
            calc_cons = consensus(seqs, threshold,
                                  free_endgaps=True, end_gap_char_out='.')
            self.assertEqual(calc_cons, free_cons)
            calc_cons2 = consensus(seqs, threshold,
                                   free_endgaps=False, end_gap_char_out='.')
            self.assertEqual(calc_cons2, cons)
            calc_cons3 = consensus(seqs, threshold,
                                   free_endgaps=True, strip_gaps=True)
            self.assertEqual(calc_cons3, strip_cons)
