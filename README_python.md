# mktest.py - Python3 McDonald-Kreitman Test

This is a Python3 conversion of the original Ruby `mkTest.rb` script by Andrew Kern. It implements the McDonald-Kreitman test for detecting positive selection in protein-coding genes.

## Features

- **McDonald-Kreitman Test**: Standard test comparing fixed differences vs polymorphisms for synonymous and non-synonymous changes
- **Polarized MK Test**: Extended analysis using an outgroup for directional inference
- **Population Genetics Toolkit**: Core functionality for sequence analysis (simplified from original)
- **FASTA Input**: Reads aligned coding sequences in FASTA format

## Requirements

- Python 3.6+
- **scipy** (required for accurate Fisher's exact test)

**Important**: This script requires scipy for optimal statistical accuracy. Install it with:
```bash
pip install scipy
```

The script will run without scipy but will use a less precise manual implementation and display warnings.

## Usage

### Basic McDonald-Kreitman Test
```bash
python3 mktest.py ingroup.fa outgroup.fa
```

### Polarized McDonald-Kreitman Test
```bash
python3 mktest.py ingroup.fa outgroup.fa -p outgroup2.fa
```

### Verbose Mode (shows implementation details)
```bash
python3 mktest.py ingroup.fa outgroup.fa -v
# or
python3 mktest.py ingroup.fa outgroup.fa --verbose
```

This will show which Fisher's exact test implementation is being used and provide additional analysis details.

## Input Format

- **Coding sequences**: Must be in-frame (first base = first codon position)
- **FASTA format**: Standard FASTA files with sequence headers
- **Alignment**: Sequences should be aligned

## Example

Using the provided test data:
```bash
python3 mktest.py kreitmanAdh.fa mauritianaAdh.fa
```

Expected output:
```
aaFix   aaPoly  silFix  silPoly FET_p_val
6       1       8       8       1.76e-01
```

Where:
- **aaFix**: Amino acid (non-synonymous) fixations between species
- **aaPoly**: Amino acid (non-synonymous) polymorphisms within species
- **silFix**: Silent (synonymous) fixations between species  
- **silPoly**: Silent (synonymous) polymorphisms within species
- **FET_p_val**: Fisher's exact test p-value

## Interpretation

The McDonald-Kreitman test evaluates whether the ratio of non-synonymous to synonymous changes differs between fixed differences (between species) and polymorphisms (within species). 

- **Neutral evolution**: Similar ratios for fixations and polymorphisms
- **Positive selection**: Excess of non-synonymous fixations relative to polymorphisms
- **Purifying selection**: Deficit of non-synonymous fixations relative to polymorphisms

## Differences from Original Ruby Version

This Python conversion includes the core McDonald-Kreitman functionality but is streamlined compared to the original 4,200+ line Ruby script. Key differences:

1. **Simplified codon analysis**: Uses direct codon comparison instead of full path algorithms
2. **Basic statistics**: Core population genetics measures without full toolkit
3. **Dependency management**: Optional scipy integration for statistical tests
4. **Code structure**: More modular, object-oriented Python design

## Files

- `mktest.py`: Main script
- `kreitmanAdh.fa`: Test data (D. melanogaster Adh)
- `mauritianaAdh.fa`: Test data (D. mauritiana Adh)
- `codonMatrix.txt`: Codon distance matrix (used by original Ruby version)

## Notes

- The script expects the `codonMatrix.txt` file for full functionality (from original Ruby version)
- Input sequences must be properly aligned coding sequences
- The statistical tests provide approximate results; scipy improves accuracy
- For complex evolutionary scenarios, consider using more sophisticated tools

## Citation

If you use this tool in research, please cite the original Ruby implementation by Andrew Kern and mention this Python conversion.

## License

This code maintains the spirit of the original Ruby implementation - free to use, modify, and distribute with appropriate attribution.