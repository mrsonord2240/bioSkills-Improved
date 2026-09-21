# pysam Consensus, Comparison and Header Dict

Moved verbatim from `SKILL.md` (2026-09-21). Read for the pysam pileup consensus (teaching only; use `samtools consensus` for real work), the per-position comparison to the reference, and the header dict for writing a BAM.

### Generate Simple Consensus
```python
import pysam
from collections import Counter

def consensus_at_position(bam, chrom, pos):
    bases = Counter()
    for pileup in bam.pileup(chrom, pos, pos + 1, truncate=True):
        if pileup.pos == pos:
            for read in pileup.pileups:
                if not read.is_del and not read.is_refskip:
                    bases[read.alignment.query_sequence[read.query_position]] += 1
    if bases:
        return bases.most_common(1)[0][0]
    return 'N'

with pysam.AlignmentFile('input.bam', 'rb') as bam:
    consensus = consensus_at_position(bam, 'chr1', 1000000)   # 0-based position
    print(f'Consensus at chr1:{1000000 + 1} = {consensus}')
```

### Build Consensus Sequence (Pedagogical Only)

The Python majority-vote consensus below is illustrative, NOT production. `samtools consensus` is Bayesian, quality-aware, and platform-aware; majority vote weights every base equally and produces wrong calls on low-coverage / low-quality regions. Use for teaching pileup iteration mechanics; use `samtools consensus` for any real consensus.

`pileup()` filters before the vote (pysam 0.24.1 defaults): bases with quality < 13, unmapped / secondary / QC-fail / duplicate reads, and orphan reads are dropped, and overlapping mates are counted once. `max_depth` defaults to 8000 (deeper columns are subsampled), so raise it.

`build_consensus` returns exactly `end - start` characters, so index `i` is reference position `start + i`. `pileup()` skips uncovered columns; the function starts from all-`N` and fills only the columns it sees. Building the string by appending per pileup column shifts everything after the first coverage gap (581 false differences vs 1 true on the real chr22 slice). A column deleted in every read counts no base, so it is `N`; insertions are ignored.

```python
import pysam
from collections import Counter

def build_consensus(bam_path, chrom, start, end, min_depth=3):
    """Majority vote over [start, end), 0-based half-open; N where depth < min_depth."""
    consensus = ['N'] * (end - start)

    with pysam.AlignmentFile(bam_path, 'rb') as bam:
        for pileup in bam.pileup(chrom, start, end, truncate=True, max_depth=1_000_000):
            bases = Counter()
            for read in pileup.pileups:
                if not read.is_del and not read.is_refskip:
                    base = read.alignment.query_sequence[read.query_position]
                    bases[base.upper()] += 1

            if sum(bases.values()) >= min_depth:
                consensus[pileup.reference_pos - start] = bases.most_common(1)[0][0]

    return ''.join(consensus)
```

### Compare Consensus to Reference (Python)
```python
def compare_to_ref(bam_path, ref_path, chrom, start, end, min_depth=3):
    """[(1-based position, ref base, consensus base)] for called bases that differ from the reference."""
    consensus = build_consensus(bam_path, chrom, start, end, min_depth)
    with pysam.FastaFile(ref_path) as ref:
        reference = ref.fetch(chrom, start, end).upper()   # soft-masked FASTA is lowercase
    return [(start + i + 1, r, c)
            for i, (c, r) in enumerate(zip(consensus, reference))
            if c != 'N' and c != r]
```
A reference `N` or IUPAC code never equals a called base, so those positions are listed as differences (a 30-base `N` run gave 30 entries). Ties (50/50 columns) go to the first base counted; `samtools consensus` calls them `N` (or an IUPAC code with `--ambig`) and weights bases by quality, so expect a few different calls at het columns and at shallow, low-quality columns (chr22 slice 1952-4617 at `-d 3`: 1 difference by majority vote and by `-m simple --call-fract 0.5 --min-BQ 13`, 2 by the default Bayesian mode; at the default `-d 1`: 5 and 4).

### Header Dict for Writing a BAM (not a .dict file)
`pysam.AlignmentFile(..., 'wb', header=header)` takes this dict. It has no `M5`, so it is not a sequence dictionary: use `samtools dict` for that.
```python
import pysam

def create_dict_header(fasta_path):
    header = {'HD': {'VN': '1.6', 'SO': 'unsorted'}, 'SQ': []}

    with pysam.FastaFile(fasta_path) as ref:
        for name in ref.references:
            length = ref.get_reference_length(name)
            header['SQ'].append({'SN': name, 'LN': length})

    return header

header = create_dict_header('reference.fa')
for sq in header['SQ'][:5]:
    print(f'{sq["SN"]}: {sq["LN"]:,} bp')
```
