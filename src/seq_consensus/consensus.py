from typing import Container, Iterable, Mapping, Optional, Tuple, Union
import numpy as np


IUPAC_map = {
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "R": "AG",
    "Y": "CT",
    "S": "CG",
    "W": "AT",
    "K": "GT",
    "M": "AC",
    "B": "CGT",
    "D": "AGT",
    "H": "ACT",
    "V": "ACG",
    "N": "ACGT",
}


def consensus(
    sequences: Iterable[str],
    threshold: float = 0,
    alphabet_map: Mapping[str, str] = IUPAC_map,
    free_endgaps: bool = True,
    strip_gaps: bool = False,
    gap_char: str = "-",
    end_gap_char: str = "-",
    gap_char_out: str = "-",
    end_gap_char_out: str = "-",
    maybe_gap_char: str = "?",
):
    """
    Calculates the consensus sequence from a multiple alignment, represented by
    an iterable of same-length strings. By default, DNA is assumed (sequences
    may contain IUPAC ambiguities).
    Based on the letter frequencies in a column and `threshold`, the consensus
    letter will be decided.
    Ambiguous letters (such as IUPAC degeneracy codes) are split into the
    corresponding letters, whereby each partially contributes
    (with frequency 1/N) to the letter frequency.
    For example, the DNA ambiguity "Y" contributes half to "A" and half to "T".

    The behaviour of this function is the same as in the Geneious software
    and very similar to the DECIPHER R package.
    Details at
    https://assets.geneious.com/manual/2022.0/static/GeneiousManualse45.html
    and
    http://www2.decipher.codes/index.html

    Args:
        sequences (optional): iterable str (expects equal sequence lengths,
            will fail otherwise)
        threshold: Number between 0 and 1 indicating the proportion of all
            sequences that need the given letter in a column in order to be
            accepted as consensus letter. If the frequency is below the
            threshold, the consensus will be an ambiguity code representing a
            combination of letters whose cumulative frequency is above the
            threshold. With a threshold of 0, the most frequent base will
            always be chosen as consensus, while with a threshold of 1,
            100% of the sequences need the same letter in order to obtain an
            unambiguous consensus.
        alphabet_map: Alphabet map with keys being letters of the expected
            alphabet and values being the letters they translate to.
            This can either be the same letter (unambiguous) or a sorted
            string of ambiguous letters.
        free_endgaps: If True, end gaps will not be taken into account when
            forming the consensus. If False, there will be no distinction
            between internal and end gaps, even if `end_gap_char` is different
            from `gap_char`. Also, `end_gap_char_out` will never appear in the
            output.
        strip_gaps: If True, gaps will not be included in the consensus
        gap_char: Gap character (default: '-').
        end_gap_char: End gap character to be expected in the input
            (default: '-').
            If equal to `gap_char` (which is the default), end gaps will be
            automatically recognized when adding a sequence with `add`.
            If different from `gap_char`, end gaps will be assumed to be
            present in the input and no further end gap parsing is done.
            It is the users's responsibility to make sure the terminal gaps
            are correctly annotated.
        gap_char_out: Gap character to use in the consensus (default: '-')
        end_gap_char_out: End gap character to use in the consensus
            (default: '-'). This character will only be returned if
            `free_endgaps` is True, otherwise terminal gaps are treated like
            internal gaps.
        maybe_gap_char: Character that will be set for columns that have mixed
            gaps and letters, where gaps are frequent enough that the
            situation is ambiguous at the given threshold. Only if the most
            frequent characters are all valid letters or all gaps, an
            unambiguous consensus call is possible.

    Returns (str):
        The consensus sequence
    """
    freqs = AlignmentFrequencies(
        sequences,
        alphabet_map=alphabet_map,
        free_endgaps=free_endgaps,
        gap_char=gap_char,
        end_gap_char=end_gap_char,
    )
    return freqs.consensus(
        threshold=threshold,
        strip_gaps=strip_gaps,
        gap_char_out=gap_char_out,
        end_gap_char_out=end_gap_char_out,
        maybe_gap_char=maybe_gap_char,
    )


class AlignmentFrequencies(object):
    """
    This object holds column-wise letter frequencies, which are obtained from
    a multiple sequence alignment. The object can be constructed given a
    list or other iterator of aligned sequences, or sequences can be gradually
    added using `add`.

    Args:
        sequences: iterable str
            expects equal sequence lengths, will fail otherwise
        free_endgaps: If True, end gaps will not be taken into account.
            Note that if `end_gap_char` is different
            from `gap_char` and `free_endgaps` is not True, end gaps are
            converted to internal gaps.
            Either way, `end_gap_char` will never occur in the output.
        alphabet_map: Alphabet map with keys being letters of the expected
            alphabet and values being the letters they translate to.
            This can either be the same letter (unambiguous) or a sorted
            string of ambiguous letters.
        gap_char: Gap character
            (default: '-')
        end_gap_char: End gap character in the sequences.
            (default: '-').
            If different from `gap_char`, it needs to be correctly specified.
            Otherwise, there will be an error.
    """

    def __init__(
        self,
        sequences: Optional[Iterable[str]] = None,
        alphabet_map: Mapping[str, str] = IUPAC_map,
        free_endgaps: bool = True,
        gap_char: str = "-",
        end_gap_char: str = "-",
    ):
        self.alphabet_map = alphabet_map
        self.free_endgaps = free_endgaps
        self.gap_char = ord(gap_char)
        self._end_gap_char = self.end_gap_char = ord(end_gap_char)
        self._convert = free_endgaps and gap_char == end_gap_char
        if self._convert:
            # needs to be different from gap_char, but will not appear
            # in the output
            self._end_gap_char = 0
        # frequencies
        self._freqs = None
        self._n = 0
        # initialize alphabet
        gap_chars = (
            [self.gap_char]
            if self.gap_char == self._end_gap_char
            else [self.gap_char, self._end_gap_char]
        )
        _alphabet = sorted(alphabet_map.items(), key=lambda item: len(item[1]))
        self._alphabet_chars = np.array(
            bytearray([ord(k) for k, _ in _alphabet] + gap_chars)
        )
        self._non_ambig_i = [i for i, item in enumerate(_alphabet) if len(item[1]) == 1]
        # TODO: item[0] not strictly required
        self._ambigs = [
            (i, item[0], item[1])
            for i, item in enumerate(_alphabet)
            if len(item[1]) > 1
        ]
        self._gap_i = len(_alphabet)
        # add sequences if supplied
        if sequences is not None:
            for seq in sequences:
                self.add(seq)

    def _init_matrix(self, n: int):
        m = len(self._alphabet_chars)
        self._freqs = np.zeros(m * n, dtype=np.uint64).reshape(m, n)

    def add(self, seq: str):
        """Adds a sequence to the internal consensus matrix."""
        # convert to bytearray
        if isinstance(seq, str):
            seq = bytearray(seq, "utf-8")
        else:
            seq = bytearray(seq)
        # to uppercase
        seq = seq.upper()
        # from now on, we need a numpy array
        # seq = np.frombuffer(seq, dtype=np.uint8)
        seq = np.array(seq)

        # initialize consensus matrix if necessary
        if self._freqs is None:
            self._init_matrix(len(seq))
        assert (
            len(seq) == self._freqs.shape[1]
        ), "The input sequences do not all have a consistent length"

        if self._convert:
            _convert_endgaps(
                seq, gap_char=self.gap_char, end_gap_char=self._end_gap_char
            )

        # add to consensus matrix
        idx = np.array(seq) == self._alphabet_chars[:, None]
        if not np.all(np.any(idx, axis=0)):
            not_found = np.where(~np.any(0))[0]
            letters = np.unique(np.array(seq)[not_found])
            raise AlphabetLookupError([chr(letter) for letter in letters])
        self._freqs[idx] += 1
        self._n += 1

    def coverage(self) -> np.array:
        """
        Returns a numpy.array(dtype=numpy.float64) with the fraction of
        non-gap sequences at any position. Note that if the object was
        constructed with `free_endgaps`, terminal gaps will not count as gaps
        at all, potentially resulting in high coverage values.
        """
        matrix = self._normalize_endgaps(self._freqs)
        nongaps = np.sum(
            matrix[
                np.arange(self._gap_i),
            ],
            axis=0,
        )
        gaps = matrix[
            self._gap_i,
        ]
        all = nongaps + gaps
        cov = np.divide(
            nongaps, all, out=np.zeros(len(nongaps), dtype=np.float64), where=all != 0
        )
        return cov

    def matrix(self) -> Tuple[str, np.array]:
        """
        Returns a tuple of the letters ("row names" in the matrix)
        and the corresponding letter frequencies at each position.
        this includes all ambiguous letters as well, even if their frequency
        is zero.
        """
        letters = bytearray(self._alphabet_chars[: self._gap_i + 1]).decode("utf-8")
        matrix = self._normalize_endgaps(self._freqs)
        return letters, matrix

    def normalized_matrix(self) -> Tuple[str, np.array]:
        """
        Like `matrix`, but returns a reduced matrix with frequences from
        ambiguous letters added to the non-ambiguous letters and possible end
        gaps counts removed. This is the matrix ultimately used for calling the
        consensus.
        """
        matrix = self._normalize()
        chars = self._alphabet_chars[self._non_ambig_i + [self._gap_i]]
        return bytearray(chars).decode("utf-8"), matrix

    def num_seqs(self) -> int:
        """
        Returns the number of added sequences
        """
        return self._n

    def _normalize(self):
        """
        "Normalizes" the consensus matrix by removing "free" end gaps
        (or merging counts from two different gap characters), and then
        distributing the frequencies of ambiguous letters to their equivalents.
        Returns the resulting matrix with irrelevant rows removed.
        """
        freqs = self._normalize_endgaps(self._freqs)
        return self._replace_ambigs(freqs)

    def _normalize_endgaps(self, matrix):
        """
        If `free_endgaps` was True, remove end gaps from the frequency matrix.
        Otherwise, convert them to internal gaps if they are different (since
        they are treated like internal gaps, and we don't want to have two
        different characters for them).
        """
        # TODO: may be called too often, maybe modify self._freqs instead
        if self.free_endgaps:
            # ignore end gaps completely
            matrix = np.delete(matrix, self._gap_i + 1, axis=0)
        elif self.end_gap_char != self.gap_char:
            # ...or convert end gaps to internal gaps
            # (in case they were present in the input)
            matrix[self._gap_i, :] += matrix[self._gap_i + 1, :]
            matrix = np.delete(matrix, self._gap_i + 1, axis=0)
        return matrix

    def _replace_ambigs(self, matrix):
        """
        Replace ambiguities by corresponding (fractional) letter frequencies.
        If ambiguities are present, this will convert the
        consensus matrix from numpy.uint64 to to numpy.float64.
        """
        if len(self._ambigs) > 0:
            # found ambiguities -> convert to float
            matrix = matrix.astype(np.float64)
            for i, _, equivalents in self._ambigs:
                weight = 1 / len(equivalents)
                for letter2 in equivalents:
                    i2 = np.where(self._alphabet_chars == ord(letter2))[0]
                    if i2.size == 0:
                        raise AlphabetLookupError(chr(letter2))
                    i2 = i2[0]
                    matrix[i2, :] += matrix[i, :] * weight
            # keep only non-ambiguous columns
            matrix = np.delete(matrix, [i for i, _, _ in self._ambigs], axis=0)
        return matrix

    def column_freqs(self) -> Iterable[Iterable[Tuple[str, int]]]:
        """
        Returns an iterator over all columns, whereby each element is again
        an iterator over letters and frequencies.

        Example:

            >>> from seq_consensus import AlignmentFrequencies
            >>> freqs = AlignmentFrequencies(['AG', 'AR'])
            >>> for colfreq in freqs.column_freqs():
            >>>     print(dict(colfreq))
            {'A': 2}
            {'G': 1, 'R': 1}
        """
        return _column_freqs(*self.matrix())

    def normalized_column_freqs(
        self,
    ) -> Iterable[Iterable[Tuple[str, Union[int, float]]]]:
        """
        Works the same as `column_freqs`, but letter frequencies for
        ambiguities are split to the corresponding letters.

        Example:

            >>> from seq_consensus import AlignmentFrequencies
            >>> freqs = AlignmentFrequencies(['AG', 'AR'])
            >>> for colfreq in freqs.normalized_column_freqs():
            >>>     print(dict(colfreq))
            {'A': 2.0}
            {'A': 0.5, 'G': 1.5}
        """
        return _column_freqs(*self.normalized_matrix())

    def consensus(
        self,
        threshold: float = 0,
        strip_gaps: bool = False,
        gap_char_out: str = "-",
        end_gap_char_out: str = "-",
        maybe_gap_char: str = "?",
    ):
        """
        Calls the consensus sequence given the internal consensus matrix.
        """
        alphabet_inv = _invert_alphabet(self.alphabet_map)
        _gap_chars = [chr(self.gap_char)]
        column_freqs = (list(f) for f in self.normalized_column_freqs())
        cons_letters = (
            _consensus_letter(
                freqs,
                inverse_alphabet=alphabet_inv,
                threshold=threshold,
                gap_chars=_gap_chars,
                maybe_gap_char=maybe_gap_char,
            )
            if len(freqs) > 0
            else end_gap_char_out
            for freqs in column_freqs
        )
        if strip_gaps:
            cons_letters = (
                lt
                for lt in cons_letters
                if lt != self.gap_char and lt != end_gap_char_out
            )
        elif self.gap_char != gap_char_out:
            # convert gaps
            cons_letters = (
                gap_char_out if letter == self.gap_char else letter
                for letter in cons_letters
            )
        # done
        return "".join(cons_letters)


def _column_freqs(letters, matrix) -> Iterable[Iterable[Tuple[str, Union[int, float]]]]:
    for freqs in matrix.T:
        yield ((letter, freq) for letter, freq in zip(letters, freqs) if freq > 0)


def _convert_endgaps(seq: np.array, gap_char, end_gap_char):
    """
    Replaces terminal gaps by another character.
    The sequence is provided as numpy array.
    """
    i0 = np.argmax(seq != gap_char)
    i1 = len(seq) - np.argmax(seq[::-1] != gap_char) - 1
    seq[:i0] = end_gap_char
    seq[i1 + 1 :] = end_gap_char


def _invert_alphabet(alphabet_map):
    return dict(zip(alphabet_map.values(), alphabet_map.keys()))


def _consensus_letter(
    letter_freqs: Iterable[Tuple[str, Union[float, int]]],
    inverse_alphabet: Mapping[str, str],
    threshold: float = 0.5,
    gap_chars: Container[str] = ("-", "."),
    maybe_gap_char: str = "?",
) -> Optional[str]:
    letter_freqs = ((lt, f) for lt, f in letter_freqs if f > 0)
    sorted_freqs = sorted(letter_freqs, key=lambda x: (-x[1], x[0]))
    freqsum = sum(freq for _, freq in sorted_freqs)
    freq_threshold = threshold * freqsum
    letters = []
    cum_freq = 0
    prev_freq = 0
    n_gapchars = 0
    for letter, freq in sorted_freqs:
        # note: freq < prev_freq makes sure that ties are not broken,
        # following the "either all or none" rule of Geneious
        if cum_freq >= freq_threshold and freq < prev_freq:
            break
        cum_freq += freq
        if letter in gap_chars:
            n_gapchars += 1
        letters.append(letter)
        prev_freq = freq
    if len(letters) == 0:
        return None
    if 0 < n_gapchars < len(letters):
        return maybe_gap_char
    letters = "".join(sorted(letters))
    try:
        return inverse_alphabet[letters]
    except KeyError:
        if letters in gap_chars:
            return letters
        raise AlphabetLookupError(letters)


class AlphabetLookupError(Exception):
    def __init__(self, letters):
        self.letters = letters
        if len(letters) == 1:
            m = "The letter '{}' is not in the ambiguity lookup dictionary.".format(
                letters
            )
        else:
            m = "The letters '{}' are not in the ambiguity lookup dictionary.".format(
                letters
            )
        super().__init__(m)
