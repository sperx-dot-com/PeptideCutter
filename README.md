# PeptideCutter
A Python tool for calculating monoisotopic masses of protein sequences from FASTA files, with support for peptide fragmentation, methionine oxidation variants, and multiple charge states.

## Features

- **Mass Calculation**: Accurate monoisotopic mass calculation using pyOpenMS
- **Fragmentation Analysis**: 
  - Cut after specific residues (e.g., tryptic digestion simulation)
  - C-terminal truncation (e.g., LysN digestion)
  - All possible cut combinations
- **Methionine Oxidation**: Optional calculation of oxidized methionine variants
- **Charge States**: Automatic m/z calculation for charge states z=1 to z=7
- **Export**: Results exported to Excel format with organized data

## Requirements

```bash
pip install pandas pyopenms biopython tabulate openpyxl
```

## Usage

### Basic Usage

```python
from protein_mass_calculator import calculate_masses

# Simple analysis
calculate_masses("peptides.fasta")
```

### Advanced Usage

```python
# With fragmentation and oxidation analysis
calculate_masses(
    fasta_file="peptides.fasta",
    output_file="results.xlsx",
    cut_after_residues=['R', 'K'],  # Tryptic cuts
    c_terminal_cuts=['K'],      # Remove basic residues from C-terminus
    calculate_oxidized=True           # Include Met oxidation variants
)
```

### Configuration Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `fasta_file` | str | Required | Path to input FASTA file |
| `output_file` | str | `"protein_masses.xlsx"` | Path to output Excel file |
| `cut_after_residues` | list | `[]` | Amino acids to cut after (e.g., `['R', 'K']` for trypsin) |
| `c_terminal_cuts` | list | `[]` | Progressive C-terminal truncation residues |
| `calculate_oxidized` | bool | `True` | Calculate methionine oxidation variants |

## Output

The tool generates an Excel file containing:

- **ID**: Protein/fragment identifier
- **Fragment_Type**: Type of fragment (Original, Cut_after, C-terminal cut)
- **Sequence**: Amino acid sequence
- **Length**: Sequence length
- **#Met**: Number of methionines
- **Monoisotopic Mass [Da]**: Calculated mass
- **Oxidized Mass [Da]**: Mass with all methionines oxidized (if enabled)
- **m/z [z=1-7]**: Mass-to-charge ratios for different charge states
- **Oxidized m/z [z=1-5]**: m/z values for oxidized variants (if enabled)

## Example

```python
# Analyze tryptic peptides with oxidation
df = calculate_masses(
    fasta_file="myproteins.fasta",
    output_file="tryptic_analysis.xlsx",
    cut_after_residues=['R', 'K'],
    calculate_oxidized=True
)

print(f"Analyzed {len(df)} fragments")
```

## How It Works

1. **Sequence Processing**: Reads FASTA file and validates amino acid sequences
2. **Fragment Generation**:
   - Generates base fragments by cutting after specified residues
   - Applies progressive C-terminal truncation to all fragments
   - Creates all possible cut combinations
3. **Mass Calculation**: Uses pyOpenMS for accurate monoisotopic mass determination
4. **Charge State Analysis**: Calculates m/z for multiple protonation states
5. **Oxidation Variants**: Adds +15.9949 Da per methionine when oxidized

## Fragmentation Logic

### Cut After Residues
Generates all possible combinations of cuts after specified amino acids. For a sequence with multiple cut sites, it produces 2^n different fragment sets.

### C-terminal Truncation
Progressively removes C-terminal residues that match the specified list until a non-matching residue is encountered.

**Example:**
- Sequence: `AFLASTEED`
- C-terminal cuts: `['D', 'E']`
- Generates: `AFLASTEE`, `AFLASTE`, `AFLAST`

## Notes

- Invalid amino acid characters are automatically filtered out
- Mass calculations use the proton mass: 1.007276 Da
- Fallback calculation method available if pyOpenMS fails
- Empty or invalid sequences are skipped with warnings

## License

MIT License

## Contributing

Contributions welcome! Please open an issue or submit a pull request.

## Citation

If you use this tool in your research, please cite the underlying libraries:
- pyOpenMS: Röst et al., Nature Methods (2016)
- BioPython: Cock et al., Bioinformatics (2009)
