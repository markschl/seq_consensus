# Consensus sequences from multiple alignments

This is a simple Python 3 library focused on calculating consensus sequences. Ambiguous letters in the input are handled as well. *Numpy* is used under the hood. Currently, DNA/RNA sequences are supported.

The package additionally offers a small utility (`cons_tool`), which allows calculating consensus sequences on the commandline (see below).

## How is the consensus calculated?

The method is identical with the approach by [Geneious](https://assets.geneious.com/manual/2022.0/static/GeneiousManualse45.html) and very similar to the function [*ConsensusSequence*](https://rdrr.io/bioc/DECIPHER/man/ConsensusSequence.html) from the [`DECIPHER` R package](http://www2.decipher.codes) (options a little different).


## Example

```python
from seq_consensus import consensus

seqs = [
    'ATTGC',
    'AT-CC',
    'RT-C-'
]

consensus(seqs, threshold=0.6)
```

This returns:

```
'AT-CC'
```

## Commandline tool

The script `cons_tool` allows using the same functionality from the commandline.
An especially useful feature is the possibility to group sequences
by arbitrary regular expression pattern matched in the sequence headers:

```sh
cons_tool -k 'p:\w+' input.fasta
```

Example output (if taxonomic annotations present in header):

```
>p:Evosea consensus (n=124)
TACKATTTA--RTATTGAC-?TWA?-GKTACTAAAGCATGGGKA-T?AAA?AGGATTAGAGACCCTYGTA
>p:Chordata consensus (n=7065)
TWAYTTTA?--WAW-YWAY-YTGAA-YCCACGAAAGCTAAGAMA-CAAACTGGGATTAGATACCCCACTA
>p:Mollusca consensus (n=843)
TWAWTWTAW--WAW?WWAY-TTGAA-KYYAYGAAAKCTWRGRWA-YAAACTAGGATTAGATACCCTAYTA
>p:Chordata consensus (n=8509)
TWAYTTTA?--WAW-YMAC-TTGAA-CCCACGAAAGCTARGAMA-CAAACTGGGATTAGATACCCCACTA
>p:Platyhelminthes_ consensus (n=130)
TWAWTWTAA--WDW?TKWY-YTGAA-KYYACGAAAGYTAKGWTA-YAAACTGGGATTAGATACCCCATTA
>p:Ascomycotaconsensus (n=280)
TTAWTWTAA--WAA?TDAC-TTGAR-K??ACGAAAGCTWRGRWA-CAAACTAGGATTAGATACCCYABTA
>p:Streptophyta consensus (n=269)
TWAWTWTAW--WAW?TRAY-TTGAR-KY?ACGAAAGCTTRGRKA-CAAACTAGGATTAGATACCCTAKTA
(...)
```
