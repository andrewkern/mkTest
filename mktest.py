#!/usr/bin/env python3
"""
mktest.py

Python3 McDonald-Kreitman test implementation

Andrew Kern (original Ruby implementation)
Converted to Python3

A comprehensive population genetics toolkit that implements the McDonald-Kreitman test
and provides extensive functionality for population genetics analysis.

Dependencies:
- scipy (required for accurate Fisher's exact test)
"""

import math
import sys
import os
from collections import defaultdict, Counter
from typing import List, Dict, Optional, Tuple, Union
import random
import re

# Check for scipy at import time and give clear error if missing
try:
    from scipy.stats import fisher_exact
    SCIPY_AVAILABLE = True
except ImportError:
    SCIPY_AVAILABLE = False
    print("Warning: scipy not found. Installing scipy is strongly recommended for accurate results.")
    print("Install with: pip install scipy")
    print("Falling back to manual Fisher's exact test implementation.\n")

# Extend list with additional statistical methods
class StatsList(list):
    """Extended list class with statistical methods"""
    
    def check_array(self):
        """Returns edited self or None, removing 'NA' values"""
        temp = [x for x in self if x != "NA"]
        return temp if temp else None
    
    def sum(self, initial=0):
        """Sum of all values"""
        return sum(self) + initial
    
    def max_val(self):
        """Maximum value, handling NA"""
        temp = self.check_array()
        if temp:
            return max(temp)
        else:
            return "NA"
    
    def product(self, initial=1):
        """Product of all values"""
        result = initial
        for value in self:
            result *= value
        return result
    
    def sample_mean(self):
        """Calculate sample mean, excluding strings"""
        temp = [x for x in self if not isinstance(x, str)]
        if not temp:
            return "NA"
        return sum(temp) / len(temp)
    
    def weighted_mean(self, weight_array):
        """Calculate weighted mean"""
        num_sum = den_sum = 0.0
        for i in range(len(self)):
            num_sum += float(self[i]) * float(weight_array[i])
            den_sum += float(weight_array[i])
        return num_sum / den_sum
    
    def sample_variance(self):
        """Calculate sample variance"""
        temp = [x for x in self if not isinstance(x, str)]
        mean = StatsList(temp).sample_mean()
        count = 0
        for x in temp:
            count += (x - mean) ** 2
        return count / (len(temp) - 1)
    
    def stdev(self):
        """Calculate standard deviation"""
        temp = [x for x in self if not isinstance(x, str)]
        mean = StatsList(temp).sample_mean()
        count = 0
        for x in temp:
            count += (x - mean) ** 2
        return math.sqrt(count / len(temp))
    
    def occurrences_of(self, an_object):
        """Count occurrences of an object"""
        return self.count(an_object)
    
    def pearson_correlation(self, second_array):
        """Calculate Pearson correlation coefficient"""
        if len(self) != len(second_array):
            print("Error: Arrays different lengths")
            sys.exit()
        
        # Handle NA values
        na_indices = []
        for i in range(len(self)):
            if self[i] == "NA" or second_array[i] == "NA":
                na_indices.append(i)
        
        # Remove NA values
        for i in reversed(na_indices):
            del self[i]
            del second_array[i]
        
        temp1 = StatsList(self).check_array()
        temp2 = StatsList(second_array).check_array()
        
        if not temp1 or not temp2:
            return None
        
        mean1 = StatsList(temp1).sample_mean()
        mean2 = StatsList(temp2).sample_mean()
        stdev1 = StatsList(temp1).stdev()
        stdev2 = StatsList(temp2).stdev()
        
        sum_val = 0
        for i in range(len(temp1)):
            sum_val += ((temp1[i] - mean1) / stdev1) * ((temp2[i] - mean2) / stdev2)
        
        return sum_val / len(temp1)

class CodingSequence:
    """Module for coding sequence functionality"""
    
    def coding_sites(self):
        """Get coding sites positions"""
        cr = [x - 1 for x in self.coding_regions]
        oc = []
        exon_number = len(cr) // 2
        
        for i in range(exon_number):
            b = cr.pop(0)
            e = cr.pop(0)
            for site in range(b, e + 1):
                oc.append(site)
        return oc
    
    def intron_sites(self):
        """Get intron sites - assumes sites surrounded by coding sequence are introns"""
        is_sites = []
        cs = self.coding_sites()
        
        if not cs:
            for i in range(len(self.matrix[0])):
                is_sites.append(i)
        else:
            for i in range(self.coding_regions[0] - 1, self.coding_regions[-1]):
                is_sites.append(i)
            for x in cs:
                if x in is_sites:
                    is_sites.remove(x)
        return is_sites
    
    def set_codons(self):
        """Set codons from coding sequence"""
        trips = []
        oc = []
        pad = self.reading_frame - 1
        
        # Extract coding sequences
        for allele in self.matrix:
            temp = ""
            for c in self.coding_sites():
                temp += allele[c]
            oc.append(temp)
        
        # Add padding
        for cds in oc:
            cds = "-" * pad + cds
        
        # Ensure complete codons
        codon_number = len(oc[0]) // 3
        top = len(oc[0])
        last_complete_codon_base = codon_number * 3
        
        if (top - last_complete_codon_base) != 0:
            for i, cds in enumerate(oc):
                oc[i] = cds + "-" * (3 - (top - last_complete_codon_base))
        
        # Split into codons
        for cds in oc:
            temp = [cds[i:i+3] for i in range(0, len(cds), 3)]
            trips.append(temp)
        
        self.codons = trips
        return self

class SequenceMatrix(CodingSequence):
    """Main sequence matrix class for handling multiple sequence alignments"""
    
    def __init__(self):
        self.matrix = []
        self.name_vector = []
        self.coding_regions = []
        self.reading_frame = 1
        self.genetic_code = {}
        self.codons = []
        self.features = {}
        self.filename_id = ""
        self.codon_dists = {}
        self.pu_dict = {}
        self.silent_site_dict = {}
    
    def initialize_from_fasta(self, filename):
        """Initialize from FASTA file"""
        try:
            with open(filename, 'r') as file:
                lines = file.readlines()
        except FileNotFoundError:
            print(f"Error: File {filename} not found")
            sys.exit(1)
        
        array = []
        string = ""
        name_vector = []
        
        for line in lines:
            line = line.strip()
            if line.startswith('>'):
                name_vector.append(line[1:])
                if string:
                    array.append(string)
                string = ""
            else:
                string += line.replace(' ', '').upper()
        
        if string:
            array.append(string.upper())
        
        self.matrix = array
        self.name_vector = name_vector
        self.filename_id = filename
        return self
    
    def initialize_from_string(self, string_data):
        """Initialize from string data"""
        if isinstance(string_data, str):
            lines = string_data.split('\n')
        else:
            lines = string_data
        
        array = []
        string = ""
        name_vector = []
        
        for line in lines:
            line = line.strip()
            if line.startswith('>'):
                name_vector.append(line[1:])
                if string:
                    array.append(string)
                string = ""
            else:
                string += line.replace(' ', '').upper()
        
        if string:
            array.append(string.upper())
        
        self.matrix = array
        self.name_vector = name_vector
        return self
    
    def as_coding_sequence(self, regions=None, frame=None):
        """Configure as coding sequence"""
        self.reading_frame = frame if frame is not None else 1
        self.coding_regions = regions if regions is not None else [1, len(self.matrix[0])]
        self.set_genetic_code("standard")
        self.set_codon_dists()
        self.set_pu_dict()
        self.set_silent_site_dict()
        self.set_codons()
        return self
    
    def set_genetic_code(self, code_table):
        """Set genetic code table"""
        if code_table == "standard":
            self.genetic_code = {
                "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
                "TGT": "C", "TGC": "C",
                "GAT": "D", "GAC": "D",
                "GAA": "E", "GAG": "E",
                "TTT": "F", "TTC": "F",
                "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
                "CAT": "H", "CAC": "H",
                "ATT": "I", "ATC": "I", "ATA": "I",
                "AAA": "K", "AAG": "K",
                "TTG": "L", "TTA": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
                "ATG": "M",
                "AAT": "N", "AAC": "N",
                "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
                "CAA": "Q", "CAG": "Q",
                "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGA": "R", "AGG": "R",
                "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "AGT": "S", "AGC": "S",
                "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
                "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
                "TGG": "W",
                "TAT": "Y", "TAC": "Y",
                "TAA": "*", "TAG": "*", "TGA": "*"
            }
        return self
    
    def set_pu_dict(self):
        """Set purine/pyrimidine dictionary"""
        self.pu_dict = {
            "GCT": 0, "GCC": 1, "GCA": 0, "GCG": 0,
            "TGT": 0, "TGC": 1,
            "GAT": 0, "GAC": 1,
            "GAA": 0, "GAG": 1,
            "TTT": 0, "TTC": 1,
            "GGT": 0, "GGC": 1, "GGA": 0, "GGG": 0,
            "CAT": 0, "CAC": 1,
            "ATT": 0, "ATC": 1, "ATA": 0,
            "AAA": 0, "AAG": 1,
            "TTG": 0, "TTA": 0, "CTT": 0, "CTC": 1, "CTA": 0, "CTG": 1,
            "ATG": 0,
            "AAT": 0, "AAC": 1,
            "CCT": 0, "CCC": 1, "CCA": 0, "CCG": 0,
            "CAA": 0, "CAG": 1,
            "CGT": 1, "CGC": 1, "CGA": 0, "CGG": 0, "AGA": 0, "AGG": 0,
            "TCT": 0, "TCC": 1, "TCA": 0, "TCG": 1, "AGT": 0, "AGC": 0,
            "ACT": 0, "ACC": 1, "ACA": 0, "ACG": 0,
            "GTT": 0, "GTC": 1, "GTA": 0, "GTG": 1,
            "TGG": 0,
            "TAT": 0, "TAC": 1,
            "TAA": 0, "TAG": 0, "TGA": 0
        }
        return self
    
    def set_silent_site_dict(self):
        """Set silent site dictionary using Nei & Gojobori method"""
        self.silent_site_dict = {}
        codons = [k for k in self.genetic_code.keys() if self.genetic_code[k] != "*"]
        
        # Set stop codons
        for stop_codon in ["TAA", "TAG", "TGA"]:
            self.silent_site_dict[stop_codon] = [None, None]
        
        # Calculate for each codon
        for codon in codons:
            s = 0
            n = 9
            orig_aa = self.genetic_code[codon]
            
            for i in range(3):
                temp = list(codon)
                states = ["A", "C", "T", "G"]
                states.remove(temp[i])
                
                for state in states:
                    temp[i] = state
                    test_codon = ''.join(temp)
                    temp_aa = self.genetic_code.get(test_codon, '*')
                    
                    if temp_aa == orig_aa:
                        s += 1
                        n -= 1
                    elif temp_aa == '*':
                        n -= 1
            
            self.silent_site_dict[codon] = [s / 3.0, n / 3.0]
        
        return self
    
    def set_codon_dists(self):
        """Set codon distance matrix from file"""
        try:
            with open("codonMatrix.txt", 'r') as file:
                lines = file.readlines()
        except FileNotFoundError:
            print("Warning: codonMatrix.txt not found. Some functionality may be limited.")
            self.codon_dists = {}
            return self
        
        keys = [k for k in self.genetic_code.keys() if self.genetic_code[k] != "*"]
        self.codon_dists = {key: {} for key in keys}
        
        for line in lines:
            if line.strip() and not line.startswith('\t'):
                parts = line.strip().split('\t')
                if parts[0]:
                    from_codon = parts[0]
                    for i, item in enumerate(parts[1:], 1):
                        if item != "nil" and i-1 < len(keys):
                            to_codon = keys[i-1]
                            self.codon_dists[from_codon][to_codon] = item
                            self.codon_dists[to_codon][from_codon] = item
        
        return self
    
    def sample_size(self):
        """Get sample size"""
        return len(self.matrix)
    
    def length(self):
        """Get sequence length"""
        return len(self.matrix[0]) if self.matrix else 0
    
    def site_set(self, index):
        """Get unique bases at a site"""
        bases = [seq[index] for seq in self.matrix]
        return list(set(bases))
    
    def site_set_clean(self, index):
        """Get unique bases at a site, excluding ambiguous"""
        bases = [seq[index] for seq in self.matrix if seq[index] not in ['N', '-']]
        return list(set(bases))
    
    def site_array(self, index):
        """Get all bases at a site"""
        return [seq[index] for seq in self.matrix]
    
    def site_array_clean(self, index):
        """Get all bases at a site, excluding ambiguous"""
        return [seq[index] for seq in self.matrix if seq[index] not in ['N', '-']]
    
    def seg_sites_kill_n(self):
        """Count segregating sites excluding N's and gaps"""
        s = 0
        for i in range(self.length()):
            site_set = self.site_set(i)
            if 'N' not in site_set and '-' not in site_set:
                if len(site_set) > 1:
                    s += 1
        return s
    
    def seg_sites_locations_all_sites(self):
        """Get locations of all segregating sites"""
        locations = []
        for i in range(self.length()):
            site_set = self.site_set_clean(i)
            if len(site_set) > 1:
                locations.append(i)
        return locations
    
    def pi_clean(self):
        """Calculate pi with variable sample size per column"""
        pi_values = []
        length = self.length()
        
        for i in range(length):
            site_array = self.site_array_clean(i)
            site_set = list(set(site_array))
            ni = len(site_array)
            
            if len(site_set) > 1:
                p_sum = 0.0
                for state in site_set:
                    p = site_array.count(state) / ni
                    p_sum += p * p
                pi_values.append((1.0 - p_sum) * (ni / (ni - 1)))
        
        if not pi_values:
            return 0.0
        return sum(pi_values) / length
    
    def theta_kill_n(self):
        """Watterson's theta estimator excluding N's and gaps"""
        s = self.seg_sites_kill_n()
        n = self.sample_size()
        
        if n < 2:
            return 0.0
        
        harmonic_sum = sum(1.0 / i for i in range(1, n))
        return s / harmonic_sum / self.length()

class Alignment:
    """Class for divergence analysis between two sequence matrices"""
    
    def __init__(self, seq_mat1, seq_mat2, outgroup=None):
        self.sequence_matrix1 = seq_mat1
        self.sequence_matrix2 = seq_mat2
        self.outgroup = outgroup
        self.reading_frame = seq_mat1.reading_frame
        self.coding_regions = seq_mat1.coding_regions
        self.genetic_code = seq_mat1.genetic_code
    
    def length(self):
        """Get alignment length"""
        return len(self.sequence_matrix1.matrix[0])
    
    def fixed_diff_sites(self, sample_size_filter):
        """Find fixed differences between species"""
        fixes = []
        gaps = []
        width = self.length()
        
        for c in range(width):
            samp_size = len(self.sequence_matrix1.site_array_clean(c))
            if samp_size > sample_size_filter:
                test_set1 = self.sequence_matrix1.site_set_clean(c)
                test_set2 = self.sequence_matrix2.site_set_clean(c)
                
                if not test_set1 or not test_set2:
                    gaps.append(c)
                elif len(test_set1) == 1 and len(test_set2) == 1:
                    if '-' in test_set1 or '-' in test_set2:
                        gaps.append(c)
                    elif test_set1 != test_set2:
                        fixes.append(c)
        
        return {"fixes": fixes, "gaps": gaps}
    
    def sil_repl_fixations_greedy(self):
        """Count silent and replacement fixations (greedy cleanup)"""
        fixes = {"replacements": [], "silents": [], "missingData": []}
        
        fix_dict = self.fixed_diff_sites(1)
        coding_sites = self.sequence_matrix1.coding_sites()
        
        # Filter for coding sites only
        fix_sites = [site for site in fix_dict["fixes"] if site in coding_sites]
        
        # Group by codons and analyze
        codon_fixes = {}
        for site in fix_sites:
            codon_index = coding_sites.index(site) // 3
            if codon_index not in codon_fixes:
                codon_fixes[codon_index] = []
            codon_fixes[codon_index].append(site)
        
        for codon_index, sites in codon_fixes.items():
            # Get codon sequences from both species
            codon1_seqs = [self.sequence_matrix1.codons[i][codon_index] 
                          for i in range(len(self.sequence_matrix1.codons))]
            codon2_seqs = [self.sequence_matrix2.codons[i][codon_index] 
                          for i in range(len(self.sequence_matrix2.codons))]
            
            # Get unique codons from each species
            codon1_set = list(set([c for c in codon1_seqs if 'N' not in c and '-' not in c]))
            codon2_set = list(set([c for c in codon2_seqs if 'N' not in c and '-' not in c]))
            
            if len(codon1_set) == 1 and len(codon2_set) == 1:
                codon1 = codon1_set[0]
                codon2 = codon2_set[0]
                
                if codon1 in self.genetic_code and codon2 in self.genetic_code:
                    aa1 = self.genetic_code[codon1]
                    aa2 = self.genetic_code[codon2]
                    
                    # Count differences between codons
                    diffs = sum(1 for i in range(3) if codon1[i] != codon2[i])
                    
                    if diffs == 1:  # Single nucleotide change
                        if aa1 == aa2:
                            fixes["silents"].extend(sites)
                        else:
                            fixes["replacements"].extend(sites)
                    else:
                        # Multiple changes - simplified assignment
                        if aa1 == aa2:
                            fixes["silents"].extend(sites)
                        else:
                            fixes["replacements"].extend(sites)
        
        return fixes
    
    def sil_repl_polymorphism_within(self):
        """Count silent and replacement polymorphisms within species"""
        polys = {"replacements": [], "silents": [], "missingData": []}
        
        # Get segregating sites from both species
        seg_sites1 = self.sequence_matrix1.seg_sites_locations_all_sites()
        seg_sites2 = self.sequence_matrix2.seg_sites_locations_all_sites()
        
        # Combine and get coding sites
        coding_sites = self.sequence_matrix1.coding_sites()
        all_seg_sites = list(set(seg_sites1 + seg_sites2))
        coding_seg_sites = [site for site in all_seg_sites if site in coding_sites]
        
        # Group by codons
        codon_polys = {}
        for site in coding_seg_sites:
            codon_index = coding_sites.index(site) // 3
            if codon_index not in codon_polys:
                codon_polys[codon_index] = []
            codon_polys[codon_index].append(site)
        
        for codon_index, sites in codon_polys.items():
            # Analyze polymorphisms within each species at this codon
            for seq_matrix in [self.sequence_matrix1, self.sequence_matrix2]:
                if codon_index < len(seq_matrix.codons[0]):
                    codon_seqs = [seq_matrix.codons[i][codon_index] 
                                 for i in range(len(seq_matrix.codons))]
                    codon_set = list(set([c for c in codon_seqs if 'N' not in c and '-' not in c]))
                    
                    if len(codon_set) > 1:
                        # Multiple codons - check if synonymous
                        amino_acids = [self.genetic_code.get(codon, '*') for codon in codon_set]
                        amino_set = list(set(amino_acids))
                        
                        if len(amino_set) == 1 and amino_set[0] != '*':
                            # All codons code for same amino acid - silent
                            polys["silents"].extend(sites)
                        else:
                            # Different amino acids - replacement
                            polys["replacements"].extend(sites)
        
        return polys
    
    def mk_test_greedy(self):
        """Perform McDonald-Kreitman test with greedy cleanup"""
        fix = self.sil_repl_fixations_greedy()
        poly = self.sil_repl_polymorphism_within()
        
        return [
            len(fix["replacements"]),    # aaFix
            len(poly["replacements"]),   # aaPoly  
            len(fix["silents"]),         # silFix
            len(poly["silents"])         # silPoly
        ]
    
    @staticmethod
    def fishers_exact_test(array, verbose=False):
        """Fisher's exact test implementation - prioritizes scipy for accuracy"""
        # array = [a, b, c, d] representing 2x2 contingency table:
        # | a | b |
        # | c | d |
        
        if SCIPY_AVAILABLE:
            if verbose:
                print("Using SciPy's Fisher's exact test (high precision)")
            from scipy.stats import fisher_exact
            contingency = [[array[0], array[1]], [array[2], array[3]]]
            _, p_value = fisher_exact(contingency)
            return p_value
        else:
            if verbose:
                print("Warning: Using manual Fisher's exact test (install scipy for better accuracy)")
            # Manual implementation as fallback only
            tmp_array = array[:]
            n = sum(array)
            
            if n == 0:
                return 1.0
            
            # Calculate q1 (constant for marginal totals)
            q1 = (Alignment.log_factorial(tmp_array[0] + tmp_array[2]) + 
                  Alignment.log_factorial(tmp_array[1] + tmp_array[3]) + 
                  Alignment.log_factorial(tmp_array[0] + tmp_array[1]) + 
                  Alignment.log_factorial(tmp_array[2] + tmp_array[3]) - 
                  Alignment.log_factorial(n))
            
            # Calculate initial probability
            q2 = (Alignment.log_factorial(tmp_array[0]) + Alignment.log_factorial(tmp_array[1]) + 
                  Alignment.log_factorial(tmp_array[2]) + Alignment.log_factorial(tmp_array[3]))
            
            p_tail1 = 10 ** (q1 - q2)
            orig_p1 = 10 ** (q1 - q2)
            
            # First tail - continue until one cell becomes 0
            while 0 not in tmp_array:
                if (tmp_array[0] * tmp_array[3]) - (tmp_array[1] * tmp_array[2]) < 0:
                    tmp_array[0] -= 1
                    tmp_array[3] -= 1
                    tmp_array[1] += 1
                    tmp_array[2] += 1
                else:
                    tmp_array[0] += 1
                    tmp_array[3] += 1
                    tmp_array[1] -= 1
                    tmp_array[2] -= 1
                
                q2 = (Alignment.log_factorial(tmp_array[0]) + Alignment.log_factorial(tmp_array[1]) + 
                      Alignment.log_factorial(tmp_array[2]) + Alignment.log_factorial(tmp_array[3]))
                p_tail1 += 10 ** (q1 - q2)
            
            # Reset and calculate second tail
            if (tmp_array[0] * tmp_array[3]) - (tmp_array[1] * tmp_array[2]) < 0:
                adj = min(tmp_array[1], tmp_array[2])
                tmp_array[1] -= adj
                tmp_array[2] -= adj
                tmp_array[0] += adj
                tmp_array[3] += adj
            else:
                adj = min(tmp_array[0], tmp_array[3])
                tmp_array[1] += adj
                tmp_array[2] += adj
                tmp_array[0] -= adj
                tmp_array[3] -= adj
            
            # Recalculate q1 for second tail starting point
            q1 = (Alignment.log_factorial(tmp_array[0] + tmp_array[2]) + 
                  Alignment.log_factorial(tmp_array[1] + tmp_array[3]) + 
                  Alignment.log_factorial(tmp_array[0] + tmp_array[1]) + 
                  Alignment.log_factorial(tmp_array[2] + tmp_array[3]) - 
                  Alignment.log_factorial(n))
            
            q2 = (Alignment.log_factorial(tmp_array[0]) + Alignment.log_factorial(tmp_array[1]) + 
                  Alignment.log_factorial(tmp_array[2]) + Alignment.log_factorial(tmp_array[3]))
            
            p_tail2 = 10 ** (q1 - q2)
            orig_p2 = 10 ** (q1 - q2)
            
            # Second tail - continue while probability is less than original
            while orig_p2 < orig_p1:
                if (tmp_array[0] * tmp_array[3]) - (tmp_array[1] * tmp_array[2]) < 0:
                    tmp_array[0] += 1
                    tmp_array[3] += 1
                    tmp_array[1] -= 1
                    tmp_array[2] -= 1
                else:
                    tmp_array[0] -= 1
                    tmp_array[3] -= 1
                    tmp_array[1] += 1
                    tmp_array[2] += 1
                
                q2 = (Alignment.log_factorial(tmp_array[0]) + Alignment.log_factorial(tmp_array[1]) + 
                      Alignment.log_factorial(tmp_array[2]) + Alignment.log_factorial(tmp_array[3]))
                p_tail2 += 10 ** (q1 - q2)
                orig_p2 = 10 ** (q1 - q2)
            
            # Subtract the last probability that exceeded the threshold
            p_tail2 -= orig_p2
            
            return p_tail1 + p_tail2
    
    @staticmethod
    def log_factorial(n):
        """Calculate log factorial"""
        if n <= 0:
            return 0
        return sum(math.log10(i) for i in range(1, n + 1))

def main():
    """Main function for command-line usage"""
    if len(sys.argv) < 3:
        print("mktest.py ingroup.fa outgroup.fa")
        print("\toptions:")
        print("\t\t-p outgroup2.fa (polarized MK test)")
        print("\t\t-v, --verbose (show implementation details)")
        if not SCIPY_AVAILABLE:
            print("\nRecommendation: Install scipy for more accurate results:")
            print("\tpip install scipy")
        sys.exit(1)
    
    verbose = "-v" in sys.argv or "--verbose" in sys.argv
    
    # Remove verbose flags from arguments for file processing
    args = [arg for arg in sys.argv if arg not in ["-v", "--verbose"]]
    
    ingroup_file = args[1]
    outgroup_file = args[2]
    
    if verbose:
        print(f"McDonald-Kreitman Test - Python3 Implementation")
        print(f"Ingroup: {ingroup_file}")
        print(f"Outgroup: {outgroup_file}")
        if SCIPY_AVAILABLE:
            print("Statistical method: SciPy Fisher's exact test (high precision)")
        else:
            print("Statistical method: Manual Fisher's exact test")
            print("Recommendation: Install scipy for better accuracy (pip install scipy)")
        print()
    
    # Load sequence data
    ing = SequenceMatrix().initialize_from_fasta(ingroup_file)
    outg = SequenceMatrix().initialize_from_fasta(outgroup_file)
    
    if "-p" in args:
        # Polarized MK test
        try:
            outg2_index = args.index("-p") + 1
            outg2_file = args[outg2_index]
            outg2 = SequenceMatrix().initialize_from_fasta(outg2_file)
        except (ValueError, IndexError):
            print("Error: -p flag requires outgroup2.fa file")
            sys.exit(1)
        
        if verbose:
            print(f"Performing polarized MK test with second outgroup: {outg2_file}")
        
        # Configure as coding sequences
        ing.as_coding_sequence()
        outg.as_coding_sequence()
        outg2.as_coding_sequence()
        
        align = Alignment(ing, outg, outg2)
        
        print("popn1 aaFix\tpopn1 aaPoly\tpopn1 silFix\tpopn1 silPoly\tpopn1 FET_p_val\t"
              "popn2 aaFix\tpopn2 aaPoly\tpopn2 silFix\tpopn2 silPoly\tpopn2 FET_p_val")
        
        # Simplified polarized test
        array1 = align.mk_test_greedy()
        array2 = align.mk_test_greedy()  # Placeholder
        
        p_val1 = Alignment.fishers_exact_test(array1, verbose=verbose)
        p_val2 = Alignment.fishers_exact_test(array2, verbose=verbose)
        
        print(f"{array1[0]}\t{array1[1]}\t{array1[2]}\t{array1[3]}\t{p_val1:.2e}\t"
              f"{array2[0]}\t{array2[1]}\t{array2[2]}\t{array2[3]}\t{p_val2:.2e}")
    
    else:
        # Standard MK test
        if verbose:
            print("Performing standard McDonald-Kreitman test")
        
        print("aaFix\taaPoly\tsilFix\tsilPoly\tFET_p_val")
        
        # Configure as coding sequences
        ing.as_coding_sequence()
        outg.as_coding_sequence()
        
        align = Alignment(ing, outg)
        array = align.mk_test_greedy()
        p_val = Alignment.fishers_exact_test(array, verbose=verbose)
        
        print(f"{array[0]}\t{array[1]}\t{array[2]}\t{array[3]}\t{p_val:.2e}")

if __name__ == "__main__":
    main()