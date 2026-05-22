#!/usr/bin/env python3
"""
SnapGene DNA File Parser
Reads SnapGene format (.dna) files and extracts amino acid sequences from CDS features.

SnapGene files are binary files with embedded XML metadata containing sequence information,
features (promoters, CDS, markers), and translations.

Author: Simon (Letarte Scientific)
Usage: python snapgene_parser.py <filename.dna> [--output output.fasta] [--feature-name NAME]
"""

import sys
import re
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Optional, Dict, List, Tuple


class SnapGeneParser:
    """Parse SnapGene format (.dna) files and extract sequence information."""
    
    def __init__(self, filepath: str):
        """
        Initialize parser with SnapGene file path.
        
        Args:
            filepath: Path to .dna file
            
        Raises:
            FileNotFoundError: If file doesn't exist
            ValueError: If file is not a valid SnapGene file
        """
        self.filepath = Path(filepath)
        
        if not self.filepath.exists():
            raise FileNotFoundError(f"File not found: {filepath}")
        
        self.file_data = self._read_file()
        self.dna_sequence = None
        self.features = []
        self.plasmid_name = None
        self.description = None
        self.notes = {}
        
        self._parse_snapgene()
    
    def _read_file(self) -> bytes:
        """Read binary file content."""
        with open(self.filepath, 'rb') as f:
            return f.read()
    
    def _extract_xml_sections(self) -> Dict[str, str]:
        """
        Extract XML sections from binary SnapGene file.
        
        SnapGene files contain multiple XML sections embedded in the binary data.
        This method finds and extracts them.
        
        Returns:
            Dictionary mapping XML type to content
        """
        xml_sections = {}
        
        # Look for XML declarations and extract content until closing tag
        xml_pattern = rb'<\?xml[^>]*\?>(.*?)(?=<\?xml|$)'
        
        matches = re.finditer(xml_pattern, self.file_data, re.DOTALL)
        
        for i, match in enumerate(matches):
            xml_content = b'<?xml' + match.group(0)[5:]  # Re-add opening tag
            try:
                # Try to identify the XML type
                if b'<Features' in xml_content:
                    xml_sections['features'] = xml_content
                elif b'<Primers' in xml_content:
                    xml_sections['primers'] = xml_content
                elif b'<EnzymeVisibilities' in xml_content or b'<CustomEnzymeSets' in xml_content:
                    xml_sections['enzymes'] = xml_content
            except:
                pass
        
        # Also try to extract text sections that might contain DNA/AA sequences
        text_section = self.file_data.decode('utf-8', errors='ignore')
        
        # Look for long DNA sequences (ATGC...)
        dna_matches = re.findall(r'[ATGCatgc]{100,}', text_section)
        if dna_matches:
            xml_sections['dna_sequence'] = dna_matches[0]
        
        return xml_sections
    
    def _parse_snapgene(self):
        """Parse SnapGene file and extract sequence/feature information."""
        xml_sections = self._extract_xml_sections()
        
        # Get plasmid name from filename
        self.plasmid_name = self.filepath.stem
        
        # Parse features XML if available
        if 'features' in xml_sections:
            self._parse_features_xml(xml_sections['features'])
        
        # Parse DNA sequence if found
        if 'dna_sequence' in xml_sections:
            self.dna_sequence = xml_sections['dna_sequence'].upper()
    
    def _parse_features_xml(self, xml_data: bytes):
        """
        Parse Features XML section.
        
        Args:
            xml_data: XML content as bytes
        """
        try:
            # Find the complete Features XML element
            match = re.search(rb'<Features[^>]*>(.*?)</Features>', xml_data, re.DOTALL)
            if not match:
                return
            
            features_xml = b'<Features>' + match.group(1) + b'</Features>'
            
            root = ET.fromstring(features_xml)
            
            # Extract all features
            for feature in root.findall('Feature'):
                feature_info = self._parse_feature_element(feature)
                if feature_info:
                    self.features.append(feature_info)
                    
        except ET.ParseError as e:
            print(f"Warning: Could not parse features XML: {e}", file=sys.stderr)
    
    def _parse_feature_element(self, feature_elem) -> Optional[Dict]:
        """
        Parse individual Feature element.
        
        Args:
            feature_elem: ElementTree Element for a Feature
            
        Returns:
            Dictionary with feature information or None
        """
        feature_info = {
            'name': feature_elem.get('name', 'Unknown'),
            'type': feature_elem.get('type', 'misc_feature'),
            'direction': feature_elem.get('directionality', '1'),
            'segments': [],
            'translation': None,
            'notes': {}
        }
        
        # Extract segment ranges
        for segment in feature_elem.findall('Segment'):
            range_str = segment.get('range', '')
            feature_info['segments'].append(range_str)
        
        # Extract properties (Q elements)
        for q_elem in feature_elem.findall('Q'):
            q_name = q_elem.get('name', '')
            
            # Get value from V child elements
            values = []
            for v_elem in q_elem.findall('V'):
                if 'text' in v_elem.attrib:
                    values.append(v_elem.get('text'))
                elif v_elem.text:
                    values.append(v_elem.text)
            
            if values:
                value = ' '.join(values)
                if q_name == 'translation':
                    feature_info['translation'] = value
                else:
                    feature_info['notes'][q_name] = value
        
        return feature_info if feature_info['name'] != 'Unknown' else None
    
    def get_cds_translations(self) -> Dict[str, str]:
        """
        Extract all CDS (Coding Sequence) translations from features.
        
        Returns:
            Dictionary mapping feature name to amino acid sequence
        """
        translations = {}
        
        for feature in self.features:
            if feature['type'] == 'CDS' and feature['translation']:
                # Remove stop codon (*) if present
                aa_seq = feature['translation'].rstrip('*')
                translations[feature['name']] = aa_seq
        
        return translations
    
    def get_all_features(self) -> List[Dict]:
        """
        Get all parsed features.
        
        Returns:
            List of feature dictionaries
        """
        return self.features
    
    def print_summary(self):
        """Print summary of parsed file."""
        print(f"\n{'='*70}")
        print(f"SnapGene File: {self.filepath.name}")
        print(f"{'='*70}\n")
        
        if self.plasmid_name:
            print(f"Plasmid Name: {self.plasmid_name}")
        
        print(f"\nFeatures Found: {len(self.features)}\n")
        
        for feature in self.features:
            print(f"  • {feature['name']}")
            print(f"    Type: {feature['type']}")
            if feature['segments']:
                print(f"    Location(s): {', '.join(feature['segments'])}")
            if feature['translation']:
                aa_seq = feature['translation'].rstrip('*')
                print(f"    Protein: {len(aa_seq)} aa")
                print(f"    Sequence: {aa_seq[:60]}{'...' if len(aa_seq) > 60 else ''}")
            if feature['notes']:
                for note_key, note_val in feature['notes'].items():
                    note_display = note_val[:50] + '...' if len(note_val) > 50 else note_val
                    print(f"    {note_key}: {note_display}")
            print()
    
    def export_fasta(self, output_file: Optional[str] = None, feature_name: Optional[str] = None) -> str:
        """
        Export amino acid sequences to FASTA format.
        
        Args:
            output_file: Output file path (default: {input_name}.fasta)
            feature_name: Specific feature to export (default: all CDS features)
            
        Returns:
            Output file path
        """
        translations = self.get_cds_translations()
        
        if not translations:
            print("No CDS translations found in file.", file=sys.stderr)
            return None
        
        if feature_name:
            if feature_name not in translations:
                print(f"Feature '{feature_name}' not found.", file=sys.stderr)
                return None
            translations = {feature_name: translations[feature_name]}
        
        if output_file is None:
            output_file = self.filepath.with_suffix('.fasta').name
        
        output_path = Path(output_file)
        
        with open(output_path, 'w') as f:
            for name, seq in translations.items():
                # Create FASTA header
                header = f">{self.plasmid_name}|{name}|{len(seq)}aa"
                f.write(header + '\n')
                
                # Write sequence in 60-character lines
                for i in range(0, len(seq), 60):
                    f.write(seq[i:i+60] + '\n')
        
        return str(output_path)


def main():
    """Command-line interface."""
    if len(sys.argv) < 2:
        print("Usage: python snapgene_parser.py <filename.dna> [--output output.fasta] [--feature FEATURE_NAME]")
        print("\nExamples:")
        print("  python snapgene_parser.py plasmid.dna")
        print("  python snapgene_parser.py plasmid.dna --output my_sequences.fasta")
        print("  python snapgene_parser.py plasmid.dna --feature VEGF165")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = None
    feature_name = None
    
    # Parse optional arguments
    i = 2
    while i < len(sys.argv):
        if sys.argv[i] == '--output' and i + 1 < len(sys.argv):
            output_file = sys.argv[i + 1]
            i += 2
        elif sys.argv[i] == '--feature' and i + 1 < len(sys.argv):
            feature_name = sys.argv[i + 1]
            i += 2
        else:
            i += 1
    
    try:
        # Parse SnapGene file
        parser = SnapGeneParser(input_file)
        
        # Print summary
        parser.print_summary()
        
        # Export to FASTA
        fasta_file = parser.export_fasta(output_file, feature_name)
        if fasta_file:
            print(f"✓ FASTA file exported: {fasta_file}\n")
        
    except FileNotFoundError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
