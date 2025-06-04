#!/usr/bin/env python3
"""
Manual Tests for Generated Gene Descriptions

This script validates that generated gene descriptions contain expected data categories
for well-known, well-studied genes across different data providers.

The tests are designed to be flexible and account for data changes over time while
ensuring that key description categories are present for important genes.
"""

import json
import os
import re
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass


@dataclass
class TestGene:
    """Represents a gene to test with expected description patterns."""
    gene_id: str
    gene_symbol: str
    expected_patterns: List[str]
    description_categories: List[str]
    description: str = ""


class DescriptionValidator:
    """Validates gene descriptions against expected patterns and categories."""
    
    def __init__(self, generated_descriptions_dir: str = "generated_descriptions"):
        self.descriptions_dir = generated_descriptions_dir
        self.results = {}
        
    def load_descriptions(self, provider: str) -> Dict[str, Dict]:
        """Load descriptions from JSON file for a given provider."""
        json_file = os.path.join(self.descriptions_dir, f"{provider}.json")
        if not os.path.exists(json_file):
            raise FileNotFoundError(f"Description file not found: {json_file}")
            
        with open(json_file, 'r') as f:
            data = json.load(f)
        
        # Convert to dict keyed by gene_id for easy lookup
        descriptions = {}
        for gene_data in data.get('data', []):
            descriptions[gene_data['gene_id']] = gene_data
        
        return descriptions
    
    def check_patterns(self, description: str, patterns: List[str]) -> Tuple[List[str], List[str]]:
        """Check which patterns are found in the description."""
        if not description:
            return [], patterns
            
        found_patterns = []
        missing_patterns = []
        
        for pattern in patterns:
            # Use case-insensitive regex matching
            if re.search(pattern, description, re.IGNORECASE):
                found_patterns.append(pattern)
            else:
                missing_patterns.append(pattern)
        
        return found_patterns, missing_patterns
    
    def check_categories(self, gene_data: Dict, categories: List[str]) -> Tuple[List[str], List[str]]:
        """Check which description categories have content."""
        found_categories = []
        missing_categories = []
        
        for category in categories:
            if gene_data.get(category) and gene_data[category].strip():
                found_categories.append(category)
            else:
                missing_categories.append(category)
        
        return found_categories, missing_categories
    
    def test_gene(self, gene: TestGene, gene_data: Dict) -> Dict:
        """Test a single gene against expected patterns and categories."""
        result = {
            'gene_id': gene.gene_id,
            'gene_symbol': gene.gene_symbol,
            'has_description': bool(gene_data.get('description')),
            'description': gene_data.get('description', 'No description available'),
            'found_patterns': [],
            'missing_patterns': [],
            'found_categories': [],
            'missing_categories': [],
            'passed': False
        }
        
        # Check patterns in main description
        description = gene_data.get('description', '')
        if description:
            found_patterns, missing_patterns = self.check_patterns(description, gene.expected_patterns)
            result['found_patterns'] = found_patterns
            result['missing_patterns'] = missing_patterns
        
        # Check description categories
        found_categories, missing_categories = self.check_categories(gene_data, gene.description_categories)
        result['found_categories'] = found_categories
        result['missing_categories'] = missing_categories
        
        # Gene passes if it has a description and either:
        # 1. Some expected patterns are found, OR
        # 2. Some expected categories have content
        result['passed'] = (
            result['has_description'] and 
            (len(result['found_patterns']) > 0 or len(result['found_categories']) > 0)
        )
        
        return result
    
    def run_tests(self, provider: str, test_genes: List[TestGene]) -> Dict:
        """Run tests for all genes in a provider."""
        print(f"\n{'='*60}")
        print(f"Testing {provider} descriptions")
        print(f"{'='*60}")
        
        try:
            descriptions = self.load_descriptions(provider)
        except FileNotFoundError as e:
            print(f"❌ Error: {e}")
            return {'provider': provider, 'error': str(e)}
        
        results = {
            'provider': provider,
            'total_genes': len(test_genes),
            'passed': 0,
            'failed': 0,
            'gene_results': []
        }
        
        for gene in test_genes:
            gene_data = descriptions.get(gene.gene_id, {})
            if not gene_data:
                print(f"❌ Gene {gene.gene_symbol} ({gene.gene_id}) not found in descriptions")
                result = {
                    'gene_id': gene.gene_id,
                    'gene_symbol': gene.gene_symbol,
                    'error': 'Gene not found in descriptions',
                    'passed': False
                }
            else:
                result = self.test_gene(gene, gene_data)
            
            results['gene_results'].append(result)
            
            if result['passed']:
                results['passed'] += 1
                status = "✅ PASS"
            else:
                results['failed'] += 1
                status = "❌ FAIL"
            
            print(f"{status} {gene.gene_symbol} ({gene.gene_id})")
            if result.get('description'):
                print(f"    Description: {result['description'][:100]}...")
            if result.get('found_patterns'):
                print(f"    Found patterns: {', '.join(result['found_patterns'])}")
            if result.get('found_categories'):
                print(f"    Found categories: {', '.join(result['found_categories'])}")
            if result.get('missing_patterns') and not result['passed']:
                print(f"    Missing patterns: {', '.join(result['missing_patterns'])}")
            if result.get('missing_categories') and not result['passed']:
                print(f"    Missing categories: {', '.join(result['missing_categories'])}")
            print()
        
        print(f"Results: {results['passed']}/{results['total_genes']} tests passed")
        return results


def define_test_genes() -> Dict[str, List[TestGene]]:
    """Define well-known genes to test for each data provider."""
    
    # Common patterns that might appear in descriptions
    go_patterns = [
        r"Predicted to enable",
        r"Involved in",
        r"Located in",
        r"Part of",
        r"Acts upstream"
    ]
    
    expression_patterns = [
        r"Is expressed in",
        r"Expressed in"
    ]
    
    # Disease description patterns - covering all three categories
    disease_biomarker_patterns = [
        r"biomarker.*for",
        r"marker.*for",
        r"indicator.*of",
        r"associated.*with.*disease"
    ]
    
    disease_used_to_study_patterns = [
        r"Used to study",
        r"Model.*for",
        r"Research.*model",
        r"Study.*of"
    ]
    
    disease_orthology_patterns = [
        r"Human ortholog.*implicated in",
        r"ortholog.*associated.*with",
        r"Orthologous to human.*gene.*implicated",
        r"implicated.*via.*orthology"
    ]
    
    # Combined disease patterns
    disease_patterns = disease_biomarker_patterns + disease_used_to_study_patterns + disease_orthology_patterns
    
    orthology_patterns = [
        r"Orthologous to",
        r"Human ortholog"
    ]
    
    # Common description categories to check - including all three disease categories
    common_categories = [
        'go_description',
        'go_function_description', 
        'go_process_description',
        'go_component_description',
        'do_description',  # Direct disease annotations
        'do_biomarker_description',  # Disease biomarker annotations
        'do_orthology_description',  # Disease via orthology annotations
        'orthology_description',
        'tissue_expression_description'
    ]
    
    return {
        'WB': [  # C. elegans - WormBase
            # Disease via orthology - insulin receptor, aging/longevity research
            TestGene(
                gene_id='WB:WBGene00000898',
                gene_symbol='daf-2',
                expected_patterns=go_patterns + disease_orthology_patterns + orthology_patterns,
                description_categories=common_categories
            ),
            # Used to study - behavioral phenotypes
            TestGene(
                gene_id='WB:WBGene00001179',
                gene_symbol='egl-10',
                expected_patterns=go_patterns + disease_used_to_study_patterns + expression_patterns,
                description_categories=common_categories
            ),
            # Used to study - FGFR signaling, development
            TestGene(
                gene_id='WB:WBGene00001184',
                gene_symbol='egl-15',
                expected_patterns=go_patterns + disease_used_to_study_patterns,
                description_categories=common_categories
            ),
            # Disease via orthology - calcium channel, epilepsy research
            TestGene(
                gene_id='WB:WBGene00001187',
                gene_symbol='egl-19',
                expected_patterns=go_patterns + disease_orthology_patterns + disease_used_to_study_patterns,
                description_categories=common_categories
            ),
            # Disease via orthology - Notch pathway, cancer research
            TestGene(
                gene_id='WB:WBGene00003001',
                gene_symbol='lin-12',
                expected_patterns=go_patterns + disease_orthology_patterns + expression_patterns,
                description_categories=common_categories
            ),
            # Used to study - spectrin, developmental disorders
            TestGene(
                gene_id='WB:WBGene00004855',
                gene_symbol='sma-1',
                expected_patterns=go_patterns + expression_patterns + disease_used_to_study_patterns,
                description_categories=common_categories
            ),
            # Biomarker potential - APP homolog, Alzheimer's research
            TestGene(
                gene_id='WB:WBGene00000094',
                gene_symbol='apl-1',
                expected_patterns=go_patterns + disease_biomarker_patterns + disease_used_to_study_patterns,
                description_categories=common_categories
            ),
        ],
        
        'MGI': [  # Mouse - MGI
            # Used to study cancer - oncogene
            TestGene(
                gene_id='MGI:96680',
                gene_symbol='Kras',
                expected_patterns=go_patterns + disease_used_to_study_patterns,
                description_categories=common_categories
            ),
            # Biomarker potential - MYC pathway component
            TestGene(
                gene_id='MGI:96921',
                gene_symbol='Max',
                expected_patterns=go_patterns + disease_biomarker_patterns + expression_patterns,
                description_categories=common_categories
            ),
            # Used to study neuroblastoma - oncogene
            TestGene(
                gene_id='MGI:97357', 
                gene_symbol='Mycn',
                expected_patterns=go_patterns + disease_used_to_study_patterns + expression_patterns,
                description_categories=common_categories
            ),
            # Used to study cancer, biomarker - tumor suppressor
            TestGene(
                gene_id='MGI:98834',
                gene_symbol='Trp53',
                expected_patterns=go_patterns + disease_used_to_study_patterns + disease_biomarker_patterns,
                description_categories=common_categories
            ),
            # Biomarker for breast cancer - BRCA1
            TestGene(
                gene_id='MGI:104537',
                gene_symbol='Brca1',
                expected_patterns=go_patterns + disease_biomarker_patterns + disease_used_to_study_patterns,
                description_categories=common_categories
            ),
        ],
        
        'SGD': [  # Yeast - SGD
            # Disease via orthology - cytochrome P450, drug metabolism
            TestGene(
                gene_id='SGD:S000002810',
                gene_symbol='DIT2', 
                expected_patterns=go_patterns + disease_orthology_patterns + orthology_patterns,
                description_categories=common_categories
            ),
            # Used to study mitochondrial diseases - ribosomal protein
            TestGene(
                gene_id='SGD:S000001486',
                gene_symbol='MRP17',
                expected_patterns=go_patterns + disease_used_to_study_patterns + orthology_patterns,
                description_categories=common_categories
            ),
            # Disease via orthology - NAD transporter, metabolic disorders
            TestGene(
                gene_id='SGD:S000001268',
                gene_symbol='YIA6',
                expected_patterns=go_patterns + disease_orthology_patterns + orthology_patterns,
                description_categories=common_categories
            ),
            # Used to study cancer - DNA repair, BRCA1 ortholog
            TestGene(
                gene_id='SGD:S000000304',
                gene_symbol='RAD51',
                expected_patterns=go_patterns + disease_used_to_study_patterns + orthology_patterns,
                description_categories=common_categories
            ),
        ],
        
        'RGD': [  # Rat - RGD
            # Used to study cancer, biomarker - tumor suppressor
            TestGene(
                gene_id='RGD:2894',
                gene_symbol='Tp53',
                expected_patterns=go_patterns + disease_used_to_study_patterns + disease_biomarker_patterns,
                description_categories=common_categories
            ),
            # Biomarker for cardiovascular disease - renin
            TestGene(
                gene_id='RGD:3555',
                gene_symbol='Ren',
                expected_patterns=go_patterns + disease_biomarker_patterns,
                description_categories=common_categories
            ),
        ],
        
        'ZFIN': [  # Zebrafish - ZFIN
            # Used to study cancer - tumor suppressor
            TestGene(
                gene_id='ZFIN:ZDB-GENE-990415-8',
                gene_symbol='tp53',
                expected_patterns=go_patterns + disease_used_to_study_patterns + expression_patterns,
                description_categories=common_categories
            ),
            # Used to study development and disease - sonic hedgehog
            TestGene(
                gene_id='ZFIN:ZDB-GENE-980526-166',
                gene_symbol='shha',
                expected_patterns=go_patterns + disease_used_to_study_patterns + expression_patterns,
                description_categories=common_categories
            ),
        ],
        
        'FB': [  # Drosophila - FlyBase
            # Used to study cancer - tumor suppressor
            TestGene(
                gene_id='FB:FBgn0003996',
                gene_symbol='p53',
                expected_patterns=go_patterns + disease_used_to_study_patterns + expression_patterns,
                description_categories=common_categories
            ),
            # Used to study neurodegenerative diseases - alpha-synuclein
            TestGene(
                gene_id='FB:FBgn0283521',
                gene_symbol='Snca',
                expected_patterns=go_patterns + disease_used_to_study_patterns + expression_patterns,
                description_categories=common_categories
            ),
        ]
    }


def main():
    """Run all manual tests."""
    print("Manual Gene Description Validation Tests")
    print("=" * 60)
    print("This script tests that well-known genes have appropriate descriptions")
    print("with expected patterns and data categories.")
    print()
    
    validator = DescriptionValidator()
    test_genes = define_test_genes()
    
    all_results = {}
    total_passed = 0
    total_tests = 0
    
    # Run tests for each provider
    for provider, genes in test_genes.items():
        result = validator.run_tests(provider, genes)
        all_results[provider] = result
        
        if 'error' not in result:
            total_passed += result['passed']
            total_tests += result['total_genes']
    
    # Summary
    print(f"\n{'='*60}")
    print("SUMMARY")
    print(f"{'='*60}")
    print(f"Total tests: {total_tests}")
    print(f"Total passed: {total_passed}")
    print(f"Total failed: {total_tests - total_passed}")
    print(f"Success rate: {total_passed/total_tests*100:.1f}%" if total_tests > 0 else "N/A")
    
    # Provider breakdown
    print("\nBy provider:")
    for provider, result in all_results.items():
        if 'error' in result:
            print(f"  {provider}: ERROR - {result['error']}")
        else:
            success_rate = result['passed']/result['total_genes']*100 if result['total_genes'] > 0 else 0
            print(f"  {provider}: {result['passed']}/{result['total_genes']} ({success_rate:.1f}%)")
    
    print(f"\n{'='*60}")
    print("Notes:")
    print("- Tests check for flexible patterns that may evolve over time")
    print("- 'Predicted' annotations may become 'experimental' with new data")
    print("- Missing patterns don't always indicate problems if other categories are present")
    print("- Genes without descriptions are expected for some less-studied genes")


if __name__ == "__main__":
    main()
