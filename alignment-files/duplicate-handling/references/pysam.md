# Duplicate Handling: pysam

## pysam Python Alternative

### Full Pipeline
```python
import pysam

# Sort by name
pysam.sort('-n', '-o', 'namesort.bam', 'input.bam')

# Fixmate
pysam.fixmate('-m', 'namesort.bam', 'fixmate.bam')

# Sort by coordinate
pysam.sort('-o', 'coordsort.bam', 'fixmate.bam')

# Mark duplicates
pysam.markdup('coordsort.bam', 'marked.bam')

# Index
pysam.index('marked.bam')
```

### Check Duplicate Flag
```python
import pysam

with pysam.AlignmentFile('marked.bam', 'rb') as bam:
    total = 0
    duplicates = 0
    for read in bam:
        if read.is_secondary or read.is_supplementary:
            continue
        total += 1
        if read.is_duplicate:
            duplicates += 1

    print(f'Total: {total}')
    print(f'Duplicates: {duplicates}')
    print(f'Rate: {duplicates/total*100:.2f}%')
```

### Filter Out Duplicates
```python
import pysam

with pysam.AlignmentFile('marked.bam', 'rb') as infile:
    with pysam.AlignmentFile('nodup.bam', 'wb', header=infile.header) as outfile:
        for read in infile:
            if not read.is_duplicate:
                outfile.write(read)
```
