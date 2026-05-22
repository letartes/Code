# SnapGene DNA File Parser

A Python application to read SnapGene format (.dna) files and extract amino acid sequences from CDS (Coding Sequence) features.

## Overview

SnapGene files are binary files commonly used in molecular biology for plasmid design and visualization. They contain:
- DNA sequences
- Plasmid features (promoters, CDS, resistance markers, origins of replication)
- Protein translations
- Restriction enzyme sites
- Metadata (creation date, authors, etc.)

This parser extracts the structured information, particularly useful for protein engineering workflows where you need to quickly access translated sequences.

## Installation

No external dependencies required beyond Python 3.6+. Uses only standard library:
- `xml.etree.ElementTree` - XML parsing
- `re` - Regular expressions for binary data extraction
- `pathlib` - File path handling

## Usage

### Basic Usage

Extract all CDS translations from a SnapGene file:

```bash
python snapgene_parser.py plasmid.dna
```

This will:
1. Display a summary of all features found
2. Print all protein sequences
3. Generate a FASTA file with the same name as the input

### Export to Custom FASTA File

```bash
python snapgene_parser.py plasmid.dna --output my_sequences.fasta
```

### Extract Specific Feature

```bash
python snapgene_parser.py plasmid.dna --feature VEGF165
```

## Examples

### Example 1: Parse Your VEGF165 Plasmid

```bash
python snapgene_parser.py _PNY3247__pUC-GW-Kan_CMV_2PA-VEGF165_8A.dna
```

Output:
```
======================================================================
SnapGene File: _PNY3247__pUC-GW-Kan_CMV_2PA-VEGF165_8A.dna
======================================================================

Plasmid Name: _PNY3247__pUC-GW-Kan_CMV_2PA-VEGF165_8A

Features Found: 8

  • VEGF165
    Type: CDS
    Location(s): 1025-1600
    Protein: 165 aa
    Sequence: MNFLLSWVHWSLALLLYLHHAKWSQAAPMAEGGGQNHHEVVKFMDVYQR...
    
  • CMV promoter
    Type: promoter
    Location(s): 786-989
    ...
```

### Example 2: Batch Processing Multiple Plasmids

```bash
for file in *.dna; do
    python snapgene_parser.py "$file" --output "sequences/${file%.dna}.fasta"
done
```

## API Usage

Use the parser programmatically in your own Python code:

```python
from snapgene_parser import SnapGeneParser

# Parse file
parser = SnapGeneParser('plasmid.dna')

# Get all CDS translations
translations = parser.get_cds_translations()
for feature_name, sequence in translations.items():
    print(f"{feature_name}: {sequence}")

# Export to FASTA
output_file = parser.export_fasta('output.fasta')

# Get all features
features = parser.get_all_features()
for feature in features:
    if feature['type'] == 'CDS':
        print(f"CDS: {feature['name']} at {feature['segments']}")
```

## Output Format

### Console Output

The parser displays:
- File name and plasmid name
- Number of features found
- For each feature:
  - Name
  - Type (promoter, CDS, marker, polyA_signal, etc.)
  - Genomic location(s)
  - For CDS: protein length and sequence preview
  - Additional notes/descriptions

### FASTA Output

Standard FASTA format suitable for alignment, mass spec analysis, or other bioinformatics tools:

```
>plasmid_name|feature_name|165aa
MNFLLSWVHWSLALLLYLHHAKWSQAAPMAEGGGQNHHEVVKFMDVYQRSYCHPIETLVD
IFQEYPDEIEYIFKPSCVPLMRCGGCCNDEGLECVPTEESNITMQIMRIKPHQGQHIGEM
SFLQHNKCECRPKKDRARQENPCGPCSERRKHLFVQDPQTCKCSCKNTDSRCKARQLELN
ERTCRCDKPRR
```

## Features Extracted

The parser identifies and catalogues:

| Feature Type | Description | Example |
|---|---|---|
| `CDS` | Coding sequence | VEGF165 |
| `promoter` | Gene promoter | CMV promoter |
| `enhancer` | Enhancer element | CMV enhancer |
| `polyA_signal` | Polyadenylation signal | hGH poly(A) |
| `rep_origin` | Origin of replication | pBR322 ori |
| `marker_gene` | Selection marker | KanR, AmpR |

## Common Workflows

### Protein Engineering

Extract your engineered protein variants for mass spec analysis:

```bash
python snapgene_parser.py VEGF121_variant.dna --output vegf_variants.fasta
python snapgene_parser.py VEGF165_variant.dna --output vegf_variants.fasta --feature VEGF165
```

### CMC/Analytical Characterization

Get sequences for:
- Intact mass calculation
- Peptide mapping method design
- Reference standard sequences

### Quality Control

Verify that your expression plasmids contain the expected sequences.

## Troubleshooting

### "File not found" error
Ensure the path to your .dna file is correct and the file exists.

### No features found
Some SnapGene files may have minimal metadata. The parser will still attempt to extract DNA sequences from the binary data.

### "Could not parse features XML" warning
Some SnapGene file versions may use slightly different XML structures. The parser tries to be flexible but may miss some edge cases. Contact support if you encounter this.

## Limitations

- Binary format varies slightly between SnapGene versions; parser is optimized for recent versions
- Some metadata may be lost if the file uses non-standard encoding
- Does not parse graphical annotations (colors, shapes, etc.)
- Does not support writing back to SnapGene format (read-only)

## Future Enhancements

Potential additions:
- DNA sequence extraction and export
- Codon optimization analysis
- Primer design for cloning
- SnapGene file writing/modification
- Support for other plasmid formats (GenBank, etc.)

## Author

Simon (Letarte Scientific) - For protein engineering and CMC applications

## License

Use freely for research and development purposes.
