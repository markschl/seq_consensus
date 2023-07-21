#!/usr/bin/env python3

import os
from sys import stderr, stdin
from os import linesep
import re
import argparse
import numpy as np

from seq_consensus import AlignmentFrequencies, IUPAC_map,\
    AlphabetLookupError

try:
    from Bio import SeqIO
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq
    from Bio.SeqIO.FastaIO import FastaWriter
except ImportError:
    print("Biopython should be installed in order to use 'cons_tool.py'",
          file=stderr)
    exit(1)


def calc_consensus(aln_files,
                   output,
                   out_format='fasta',
                   alphabet_map=IUPAC_map,
                   threshold=0.5,
                   min_coverage=None,
                   free_endgaps=True, gap_char='-', end_gap_char='-',
                   maybe_gap_char='?',
                   key=None,
                   header='{filename} consensus',
                   unmatched_header='(other)',
                   wrap=None):
    if out_format == 'fasta':
        out = FastaWriter(output, wrap=wrap)
        def writer(head, seq): out.write_record(
            SeqRecord(Seq(seq), id=head, description=''))
    elif out_format == 'line':
        def writer(_, seq): output.write(seq + linesep)
    else:
        assert out_format == 'raw'
        def writer(_, seq): output.write(seq)

    for aln_file in aln_files:
        sequences = SeqIO.parse(aln_file, 'fasta')
        name = os.path.basename(aln_file.name)
        header_vars = {
            'filename': name,
            'filestem': os.path.splitext(name)[0],
            'path': aln_file.name
        }
        if key is None:
            freqs = AlignmentFrequencies(
                (str(r.seq) for r in sequences),
                alphabet_map=alphabet_map,
                free_endgaps=free_endgaps,
                gap_char=gap_char,
                end_gap_char=end_gap_char
            )
            cons = freqs.consensus(threshold=threshold,
                                   gap_char_out=gap_char,
                                   end_gap_char_out=end_gap_char,
                                   maybe_gap_char=maybe_gap_char)
            _header = header.format(n=freqs.num_seqs(), **header_vars)
            if min_coverage is not None:
                keep = freqs.coverage() >= min_coverage
                assert len(cons) == len(keep)
                cons = ''.join(lt for lt, include in zip(cons, keep)
                               if include)
            writer(_header, cons)
        else:
            pattern = re.compile(key)
            group_freqs = {}
            for rec in sequences:
                # obtain the group
                match = pattern.search(rec.description)
                if match is None:
                    group = unmatched_header
                else:
                    group = match.expand(header)
                # add
                if group not in group_freqs:
                    group_freqs[group] = AlignmentFrequencies(
                        alphabet_map=alphabet_map,
                        free_endgaps=free_endgaps,
                        gap_char=gap_char,
                        end_gap_char=end_gap_char
                    )
                group_freqs[group].add(str(rec.seq))
            if min_coverage is not None:
                cov = np.stack([freqs.coverage()
                                for freqs in group_freqs.values()])
                cov = np.average(cov, axis=0)
                keep = cov >= min_coverage
            else:
                keep = None
            for group, freqs in group_freqs.items():
                cons = freqs.consensus(threshold=threshold,
                                       gap_char_out=gap_char,
                                       end_gap_char_out=end_gap_char,
                                       maybe_gap_char=maybe_gap_char)
                if keep is not None:
                    assert len(cons) == len(keep)
                    cons = ''.join(lt for lt, include in zip(cons, keep)
                                   if include)
                group = group.format(n=freqs.num_seqs(), **header_vars)
                writer(group, cons)


parser = argparse.ArgumentParser(
    description="""
    This tool calculates consensus sequences.
    """)

parser.add_argument('alignment', nargs='*', type=argparse.FileType('r'),
                    default=[stdin],
                    help='''
                    FASTA sequences. Multiple files can be supplied; the
                    consensus is calculated for every of them.
                    Default: standard input
                    ''')
parser.add_argument('-o', '--output', type=argparse.FileType('w'),
                    default='-',
                    help='''
                    Consensus sequence output file.
                    Default: standard output (-)
                    ''')
parser.add_argument('-f', '--out-format',
                    choices={'line', 'raw',
                             'line/fasta', 'raw/fasta',
                             'fasta'},
                    default='line/fasta',
                    help='''
                    Output format. For single consensus sequences, the default
                    is to return only the sequence with a newline. To suppress
                    the newline, use -f raw or -f raw/fasta.
                    With multiple consensus sequences (mutliple files,
                    (--key option), FASTA is returned.
                    Use -f fasta to force FASTA output in all cases.
                    ''')
parser.add_argument('-t', '--threshold', type=float, default=0.5,
                    help='''
                    Frequency threshold indicating the minimum proportion of
                    sequences that need to have the same letter(s) in order
                    to be accepted as consensus (default: 0.5).
                    ''')
parser.add_argument('-c', '--min-coverage', type=float,
                    help='''
                    Minimum coverage, which is the fraction of non-gap bases
                    at a given alignment column. The default is not to filter
                    at all.
                    The filter is applied per file. If there are multiple
                    groups (--key), the average of the coverages from all
                    groups is calculated.
                    ''')
parser.add_argument('-T', '--type', choices={'dna', 'rna'}, default='dna',
                    help='''
                    Sequence type. DNA and RNA are currently supported
                    (default: 'dna')
                    ''')
parser.add_argument('-n', '--no-free-endgaps', action='store_true',
                    help='Do not ignore end gaps when calculating the'
                         'consensus.')
parser.add_argument('-g', '--gap-char', default='-',
                    help='Gap character in alignment (default: "-").')
parser.add_argument('-e', '--end-gap-char', default='-',
                    help='''
                    End gap character in alignment. The default assumption is
                    that end gaps are the same as gaps ('-'). The tool will
                    therefore automatically recognize end gaps (unless
                    --no-free-endgaps is specified). If terminal gaps are
                    different from internal ones in the input,
                    -e/--end-gap-char needs to be correctly specified
                    there will be an error. Also note that no automatic
                    validation of end gaps is done, they are assumed to be
                    correct in the input.
                    ''')
parser.add_argument('-k', '--key',
                    help='''
                    Optional regular expression matching part of the sequence
                    headers. The matched part is used for grouping the
                    sequences. One consensus per group is returned in FASTA
                    format.
                    Regex groups can also be used in more complex situations
                    (see -H/--header).
                    ''')
parser.add_argument('-H', '--out-header',
                    help=r'''
                    Header line of the output FASTA. With a single consensus
                    (no --key), the default is 'consensus (n={n})', where
                    {n} expands to the number of sequences in the alignment.
                    With multiple input files, the headers will be:
                    '{filename} consensus (n={n})'. Apart from 'filename',
                    the following variables are supported: {filestem}
                    (file name without extension), {path} (the full file path).
                    With the (--key option), the matched part of the sequence
                    headers is additionally included in the header. With
                    multiple files *and* --key, the header will be:
                    '{filename}: \g<0> consensus (n={n})'. Note that the
                    matched string is referenced using '\g<0>'. Regex groups
                    in --key are further referenced with '\g<1>' (1st group),
                    '\g<2>' (2nd group), etc. Note that simple backreferences
                    ('\1') do not work).
                    ''')
parser.add_argument('--unmatched-header', default='(other)',
                    help='''
                    FASTA header to use if the --key regex does not
                    match anything.
                    ''')

parser.add_argument('-w', '--wrap', type=int,
                    help='Wrap FASTA lines in the output'
                         '(default: no wrapping)')
args = parser.parse_args()

outfmt = args.out_format
if (args.key is not None or len(args.alignment) > 1):
    if 'fasta' not in outfmt:
        print('Warning: -f line or -f raw are incompatible with --key '
              'or multiple alignment files, '
              'use -f line/fasta or -f raw/fasta to silence this warning.')
    outfmt = 'fasta'
outfmt = outfmt.replace('/fasta', '')

# determine header
header = args.out_header
if header is None:
    if args.key is None:
        if len(args.alignment) > 1:
            header = '{filename} consensus (n={n})'
        else:
            header = 'consensus (n={n})'
    else:
        if len(args.alignment) > 1:
            header = r'{filename}: \g<0> consensus (n={n})'
        else:
            header = r'\g<0> consensus (n={n})'

# determine alphabet
alphabet_map = IUPAC_map
if args.type == 'rna':
    alphabet_map['U'] = 'U'
    del alphabet_map['T']
    for k, v in alphabet_map.items():
        alphabet_map[k] = v.replace('T', 'U')

try:
    calc_consensus(
        args.alignment,
        output=args.output,
        out_format=outfmt,
        alphabet_map=alphabet_map,
        threshold=args.threshold,
        min_coverage=args.min_coverage,
        free_endgaps=not args.no_free_endgaps,
        gap_char=args.gap_char,
        end_gap_char=args.end_gap_char,
        key=args.key,
        header=header,
        unmatched_header=args.unmatched_header,
        wrap=args.wrap
    )
except AlphabetLookupError as e:
    if args.gap_char == '-' and e.letters == ['.']:
        print('cons_tool.py: "." encountered in the input. If end gaps are'
              'different, specify with -e/--end-gap-char')
    else:
        print('''cons_tool.py: Invalid character(s) in the input: "{}".
              Are "--gap-char" and "--end-gap-char" correctly
              specified?'''.format(e.letters))
