import pandas as pd
from pyopenms import AASequence, EmpiricalFormula
from Bio import SeqIO
from tabulate import tabulate
import re

def calculate_masses(fasta_file, output_file="protein_masses.xlsx", 
                    cut_after_residues=None, c_terminal_cuts=None, calculate_oxidized=True):
    """
    Calculate monoisotopic masses of protein sequences in a FASTA file using pyOpenMS.
    Includes calculations for methionine oxidized variants, different charge states,
    and fragmentation analysis.
    
    Parameters:
    fasta_file (str): Path to the FASTA file
    output_file (str): Path to the output Excel file
    cut_after_residues (list): List of amino acids to cut after (e.g., ['R'])
    c_terminal_cuts (list): List of C-terminal amino acids to remove progressively (e.g., ['K'])
    calculate_oxidized (bool): Whether to calculate oxidized methionine variants (default: True)
    """
    results = []
    
    # Set default cutting parameters if not provided
    if cut_after_residues is None:
        cut_after_residues = []
    if c_terminal_cuts is None:
        c_terminal_cuts = []
    
    try:
        for record in SeqIO.parse(fasta_file, "fasta"):
            sequence_raw = str(record.seq)
            
            # Clean the sequence to contain only valid amino acid characters
            sequence = re.sub(r'[^ACDEFGHIKLMNPQRSTVWY]', '', sequence_raw.upper())
            
            if not sequence:
                print(f"Warning: Sequence {record.id} contains no valid amino acids after cleaning.")
                continue
            
            # Generate all fragments for analysis - STEP BY STEP APPROACH
            
            # Step 1: Start with original sequence and cut_after fragments
            base_fragments = []
            
            # Add original sequence
            base_fragments.append({
                'sequence': sequence,
                'fragment_type': 'Original',
                'fragment_id': record.id
            })
            
            # Add fragments by cutting after specified residues
            if cut_after_residues:
                cut_fragments = generate_cut_after_fragments(sequence, cut_after_residues)
                for i, frag in enumerate(cut_fragments):
                    base_fragments.append({
                        'sequence': frag,
                        'fragment_type': f'Cut_after_{"|".join(cut_after_residues)}',
                        'fragment_id': f'{record.id}_cut_{i+1}'
                    })
            
            # Step 2: Apply C-terminal cuts to ALL base fragments
            fragments_to_analyze = []
            
            for base_frag in base_fragments:
                # Add the base fragment itself
                fragments_to_analyze.append(base_frag)
                
                # Apply C-terminal cuts to this base fragment
                if c_terminal_cuts:
                    c_term_fragments = generate_c_terminal_cuts(base_frag['sequence'], c_terminal_cuts)
                    print(f"DEBUG: For fragment {base_frag['sequence']} (type: {base_frag['fragment_type']}), generated C-terminal cuts: {c_term_fragments}")
                    
                    for i, c_frag in enumerate(c_term_fragments):
                        fragments_to_analyze.append({
                            'sequence': c_frag,
                            'fragment_type': f"{base_frag['fragment_type']}_Cterm_cut",
                            'fragment_id': f"{base_frag['fragment_id']}_Cterm_{i+1}"
                        })
            
            # Analyze all fragments
            for fragment_info in fragments_to_analyze:
                frag_sequence = fragment_info['sequence']
                
                if not frag_sequence:
                    continue
                    
                try:
                    result = analyze_single_peptide(
                        frag_sequence, 
                        fragment_info['fragment_id'], 
                        fragment_info['fragment_type'],
                        calculate_oxidized
                    )
                    if result:
                        results.append(result)
                        
                except Exception as e:
                    print(f"Warning: Could not process fragment {fragment_info['fragment_id']}: {e}")
                    continue
        
    except FileNotFoundError:
        print(f"Error: File '{fasta_file}' not found.")
        return
    except Exception as e:
        print(f"Error processing file: {e}")
        return
    
    if not results:
        print("No sequences could be processed in the FASTA file.")
        return
    
    # Create DataFrame
    df = pd.DataFrame(results)
    
    # Sort by original ID and fragment type for better organization
    df = df.sort_values(['ID', 'Fragment_Type'])
    
    # Display table
    print("\nProtein Mass Analysis Results:")
    print(tabulate(df.head(500), headers='keys', tablefmt='grid', showindex=False))  # Show first 500 rows
    if len(df) > 500:
        print(f"\n... and {len(df) - 500} more rows")
    
    # Save to Excel
    try:
        df.to_excel(output_file, index=False)
        print(f"\nResults saved to {output_file}")
    except Exception as e:
        print(f"Error saving Excel file: {e}")
        
    return df

def analyze_single_peptide(sequence, peptide_id, fragment_type, calculate_oxidized=True):
    """
    Analyze a single peptide sequence and calculate all mass properties.
    """
    try:
        # Create peptide object from sequence
        peptide = AASequence.fromString(sequence)
        
        # Calculate monoisotopic mass
        mono_mass = peptide.getMonoWeight()
        
        # Count methionines
        met_count = sequence.count('M')
        
        # Calculate m/z values for different charge states (z=1 to z=5)
        charge_states = {}
        for z in range(1, 8):
            mz = (mono_mass + z * 1.007276) / z
            charge_states[f'm/z [z={z}]'] = round(mz, 4)
        
        # Calculate oxidized variants
        oxidized_mass = mono_mass
        oxidized_charge_states = {}
        
        if calculate_oxidized and met_count > 0:
            try:
                # Create oxidized sequence
                oxidized_sequence = ""
                for aa in sequence:
                    if aa == 'M':
                        oxidized_sequence += "M(Oxidation)"
                    else:
                        oxidized_sequence += aa
                        
                oxidized_peptide = AASequence.fromString(oxidized_sequence)
                oxidized_mass = oxidized_peptide.getMonoWeight()
                
            except Exception:
                # Fallback: Calculate oxidized mass by adding 15.9949 Da per methionine
                oxidized_mass = mono_mass + (met_count * 15.9949)
        elif not calculate_oxidized:
            # If not calculating oxidized, keep same as original
            oxidized_mass = mono_mass
        
        # Calculate oxidized m/z for all charge states (only if different from original)
        if calculate_oxidized:
            for z in range(1, 6):
                mz = (oxidized_mass + z * 1.007276) / z
                oxidized_charge_states[f'Oxidized m/z [z={z}]'] = round(mz, 4)
        
        # Prepare result dictionary
        result = {
            'ID': peptide_id,
            'Fragment_Type': fragment_type,
            'Sequence': sequence,
            'Length': len(sequence),
            '#Met': met_count,
            'Monoisotopic Mass [Da]': round(mono_mass, 4)
        }
        
        # Add oxidized mass only if calculating oxidized variants
        if calculate_oxidized:
            result['Oxidized Mass [Da]'] = round(oxidized_mass, 4)
        
        # Add all charge states
        result.update(charge_states)
        if calculate_oxidized:
            result.update(oxidized_charge_states)
        
        return result
        
    except Exception as e:
        print(f"Error analyzing peptide {peptide_id}: {e}")
        # Try fallback method
        try:
            return analyze_single_peptide_fallback(sequence, peptide_id, fragment_type, calculate_oxidized)
        except Exception as e2:
            print(f"Fallback method also failed for {peptide_id}: {e2}")
            return None

def analyze_single_peptide_fallback(sequence, peptide_id, fragment_type, calculate_oxidized=True):
    """
    Fallback method using manual mass calculation.
    """
    amino_acid_masses = {
        'A': 71.03711, 'R': 156.10111, 'N': 114.04293, 'D': 115.02694,
        'C': 103.00919, 'E': 129.04259, 'Q': 128.05858, 'G': 57.02146,
        'H': 137.05891, 'I': 113.08406, 'L': 113.08406, 'K': 128.09496,
        'M': 131.04049, 'F': 147.06841, 'P': 97.05276, 'S': 87.03203,
        'T': 101.04768, 'W': 186.07931, 'Y': 163.06333, 'V': 99.06841
    }
    
    # Calculate monoisotopic mass
    mono_mass = sum(amino_acid_masses.get(aa, 0) for aa in sequence) + 18.01056  # Add water
    met_count = sequence.count('M')
    
    # Calculate m/z values for different charge states
    charge_states = {}
    for z in range(1, 6):
        mz = (mono_mass + z * 1.007276) / z
        charge_states[f'm/z [z={z}]'] = round(mz, 4)
    
    # Calculate oxidized variant
    oxidized_mass = mono_mass + (met_count * 15.9949) if calculate_oxidized and met_count > 0 else mono_mass
    oxidized_charge_states = {}
    
    if calculate_oxidized:
        for z in range(1, 6):
            mz = (oxidized_mass + z * 1.007276) / z
            oxidized_charge_states[f'Oxidized m/z [z={z}]'] = round(mz, 4)
    
    result = {
        'ID': peptide_id,
        'Fragment_Type': fragment_type,
        'Sequence': sequence,
        'Length': len(sequence),
        '#Met': met_count,
        'Monoisotopic Mass [Da]': round(mono_mass, 4)
    }
    
    # Add oxidized mass only if calculating oxidized variants
    if calculate_oxidized:
        result['Oxidized Mass [Da]'] = round(oxidized_mass, 4)
    
    result.update(charge_states)
    if calculate_oxidized:
        result.update(oxidized_charge_states)
    
    return result

def generate_cut_after_fragments(sequence, cut_after_residues):
    """
    Generate all possible combinations of fragments by cutting after specified residues.
    For sequence with multiple cut sites, generates all combinations (2^n fragments sets).
    Excludes cuts at C-terminus as they don't generate new fragments.
    """
    from itertools import combinations
    
    fragments_sets = []
    cut_positions = []
    
    # Find all cut positions (excluding the last position - C-terminus)
    for i, aa in enumerate(sequence):
        if aa in cut_after_residues and i < len(sequence) - 1:  # Don't cut after last amino acid
            cut_positions.append(i + 1)  # Cut after this position
    
    if not cut_positions:
        return []
    
    print(f"DEBUG: Found internal cut positions: {cut_positions} in sequence: {sequence}")
    
    # Generate all possible combinations of cuts (including no cuts and all cuts)
    # We'll generate combinations from 1 to len(cut_positions) cuts
    all_fragments = []
    
    for r in range(1, len(cut_positions) + 1):  # r = number of cuts to make
        for cut_combination in combinations(cut_positions, r):
            # Generate fragments for this specific combination of cuts
            fragments_for_this_combination = []
            sorted_cuts = sorted(cut_combination)
            
            start = 0
            for cut_pos in sorted_cuts:
                if cut_pos > start:
                    fragment = sequence[start:cut_pos]
                    if fragment:  # Only add non-empty fragments
                        fragments_for_this_combination.append(fragment)
                start = cut_pos
            
            # Add remaining fragment if any
            if start < len(sequence):
                remaining = sequence[start:]
                if remaining:
                    fragments_for_this_combination.append(remaining)
            
            # Add all fragments from this combination
            all_fragments.extend(fragments_for_this_combination)
            
            print(f"DEBUG: Cut combination {cut_combination} produced fragments: {fragments_for_this_combination}")
    
    # Remove duplicates while preserving order
    unique_fragments = []
    seen = set()
    for frag in all_fragments:
        if frag not in seen:
            unique_fragments.append(frag)
            seen.add(frag)
    
    print(f"DEBUG: All unique fragments: {unique_fragments}")
    return unique_fragments

def generate_c_terminal_cuts(sequence, c_terminal_cuts):
    """
    Generate fragments by progressively removing C-terminal residues.
    For example: APAKPEED -> [APAKPEE, APAKPE, APAKP] if cutting D and E
    This will continue removing amino acids from C-terminus as long as they are in the cut list.
    """
    fragments = []
    current_sequence = sequence
    
    print(f"DEBUG: Starting C-terminal cuts for {sequence}")
    print(f"DEBUG: Cut residues: {c_terminal_cuts}")
    
    # Progressively remove C-terminal residues
    while len(current_sequence) > 0:
        last_aa = current_sequence[-1]
        print(f"DEBUG: Current sequence: {current_sequence}, last AA: {last_aa}")
        
        if last_aa in c_terminal_cuts:
            current_sequence = current_sequence[:-1]  # Remove last amino acid
            print(f"DEBUG: Removed {last_aa}, new sequence: {current_sequence}")
            if current_sequence:  # Only add if not empty
                fragments.append(current_sequence)
                print(f"DEBUG: Added fragment: {current_sequence}")
        else:
            print(f"DEBUG: {last_aa} not in cut list, stopping")
            break  # Stop if last amino acid is not in cut list
    
    print(f"DEBUG: Final fragments: {fragments}")
    return fragments

if __name__ == "__main__":
    # Configuration
    fasta_file = "peptides.fasta"  # Change this to your input file path
    output_file = "protein_masses_enhanced.xlsx"  # Change this to your desired output file
    
    # Define cutting parameters
    cut_after_residues = ['D']  # Cut after aspartic acid
    #c_terminal_cuts = ['D', 'E']  # Remove D and E from C-terminus progressively
    c_terminal_cuts = []  # No C-terminal processing
    

    # Define whether to calculate oxidized methionine variants
    calculate_oxidized = False  # Set to False to skip oxidized mass calculations
    
    print(f"Processing FASTA file: {fasta_file}")
    print(f"Cut after residues: {cut_after_residues}")
    print(f"C-terminal cuts: {c_terminal_cuts}")
    print(f"Calculate oxidized variants: {calculate_oxidized}")
    print("-" * 50)
    
    df = calculate_masses(fasta_file, output_file, cut_after_residues, c_terminal_cuts, calculate_oxidized)
    
    if df is not None:
        print(f"\nTotal fragments analyzed: {len(df)}")
        print(f"Fragment types: {df['Fragment_Type'].value_counts().to_dict()}")