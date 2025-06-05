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
class CoverageThreshold:
    """Represents expected coverage thresholds for a data provider."""
    provider: str
    description_threshold: float  # Percentage (0-100)
    go_threshold: float
    disease_threshold: float
    expression_threshold: float
    orthology_threshold: float


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

    def test_coverage_thresholds(self, provider: str, thresholds: CoverageThreshold) -> Dict:
        """Test that description coverage meets minimum thresholds."""
        print(f"\n{'='*60}")
        print(f"Testing {provider} coverage thresholds")
        print(f"{'='*60}")

        try:
            descriptions = self.load_descriptions(provider)
        except FileNotFoundError as e:
            print(f"❌ Error: {e}")
            return {'provider': provider, 'error': str(e)}

        total_genes = len(descriptions)

        # Calculate actual coverage
        actual_coverage = {
            'description': sum(1 for gene in descriptions.values() if gene.get('description')) / total_genes * 100,
            'go': sum(1 for gene in descriptions.values() if gene.get('go_description')) / total_genes * 100,
            'disease': sum(1 for gene in descriptions.values() if gene.get('do_description')) / total_genes * 100,
            'expression': (
                sum(1 for gene in descriptions.values() if gene.get('tissue_expression_description')) /
                total_genes * 100
            ),
            'orthology': (
                sum(1 for gene in descriptions.values() if gene.get('orthology_description')) /
                total_genes * 100
            )
        }

        results = {
            'provider': provider,
            'total_genes': total_genes,
            'passed_thresholds': 0,
            'failed_thresholds': 0,
            'coverage_tests': []
        }

        # Test each threshold
        thresholds_to_test = [
            ('Description', actual_coverage['description'], thresholds.description_threshold),
            ('GO', actual_coverage['go'], thresholds.go_threshold),
            ('Disease', actual_coverage['disease'], thresholds.disease_threshold),
            ('Expression', actual_coverage['expression'], thresholds.expression_threshold),
            ('Orthology', actual_coverage['orthology'], thresholds.orthology_threshold)
        ]

        for category, actual, threshold in thresholds_to_test:
            passed = actual >= threshold
            if passed:
                results['passed_thresholds'] += 1
                status = "✅ PASS"
            else:
                results['failed_thresholds'] += 1
                status = "❌ FAIL"

            result = {
                'category': category,
                'actual': actual,
                'threshold': threshold,
                'passed': passed
            }
            results['coverage_tests'].append(result)

            print(f"{status} {category}: {actual:.1f}% (threshold: {threshold:.1f}%)")

        success_rate = results['passed_thresholds'] / len(thresholds_to_test) * 100
        print(
            f"\nThreshold tests: {results['passed_thresholds']}/{len(thresholds_to_test)} "
            f"passed ({success_rate:.1f}%)"
        )

        return results


def define_coverage_thresholds() -> Dict[str, CoverageThreshold]:
    """Define minimum coverage thresholds (20% lower than current levels)."""
    return {
        'WB': CoverageThreshold('WB', 24.9, 23.5, 5.9, 9.6, 12.0),
        'MGI': CoverageThreshold('MGI', 22.1, 20.5, 7.0, 15.1, 19.8),
        'SGD': CoverageThreshold('SGD', 61.3, 61.2, 17.0, 0.0, 32.6),
        'RGD': CoverageThreshold('RGD', 31.9, 29.7, 9.3, 0.0, 26.9),
        'ZFIN': CoverageThreshold('ZFIN', 48.1, 43.8, 15.0, 21.9, 38.7),
        'FB': CoverageThreshold('FB', 34.1, 31.7, 10.2, 22.0, 20.7),
        'HUMAN': CoverageThreshold('HUMAN', 38.5, 38.5, 0.0, 0.0, 0.0),
        'XBXL': CoverageThreshold('XBXL', 71.5, 56.8, 24.5, 8.5, 66.5),
        'XBXT': CoverageThreshold('XBXT', 75.8, 71.1, 21.8, 2.1, 60.5),
    }


def define_test_genes() -> Dict[str, List[TestGene]]:
    """Define well-known genes to test for each data provider."""

    # GO patterns from config file - meaningful prefixes only
    go_patterns = [
        r"exhibits",
        r"enables",
        r"contributes to",
        r"contributes as a",
        r"contributes as an",
        r"colocalizes with",
        r"colocalizes with a",
        r"colocalizes with an",
        r"predicted to have",
        r"predicted to be a",
        r"predicted to be an",
        r"predicted to enable",
        r"predicted to contribute to",
        r"predicted to contribute as a",
        r"predicted to contribute as an",
        r"predicted to colocalize with",
        r"predicted to colocalize with a",
        r"predicted to colocalize with an",
        r"involved in",
        r"predicted to be involved in",
        r"acts upstream of with a positive effect on",
        r"acts upstream of with a negative effect on",
        r"acts upstream of or within",
        r"acts upstream of or within with a negative effect on",
        r"acts upstream of or within with a positive effect on",
        r"predicted to act upstream of with a positive effect on",
        r"predicted to act upstream of with a negative effect on",
        r"predicted to act upstream of or within",
        r"predicted to act upstream of or within with a negative effect on",
        r"predicted to act upstream of or within with a positive effect on",
        r"localizes to",
        r"located in",
        r"part of",
        r"is active in",
        r"predicted to localize to",
        r"predicted to be located in",
        r"predicted to be part of",
        r"predicted to be active in"
    ]

    expression_patterns = [
        r"is expressed in",
        r"is enriched in"
    ]

    # Disease biomarker patterns from config file
    disease_biomarker_patterns = [
        r"biomarker of"
    ]

    disease_used_to_study_patterns = [
        r"used to study"
    ]

    disease_orthology_patterns = [
        r"human ortholog\(s\) of this gene implicated in",
        r"ortholog\(s\) of this gene implicated in"
    ]

    # Human-specific disease patterns from config file
    disease_human_patterns = [
        r"implicated in"
    ]

    # Combined disease patterns
    disease_patterns = (disease_biomarker_patterns + disease_used_to_study_patterns +
                        disease_orthology_patterns + disease_human_patterns)

    orthology_patterns = [
        r"human ortholog\(s\) of this gene implicated in",
        r"ortholog\(s\) of this gene implicated in"
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
        'WB': [  # C. elegans - WormBase (4 genes - balanced)
            # Disease via orthology - insulin receptor, aging/longevity research
            TestGene(
                gene_id='WB:WBGene00000898',
                gene_symbol='daf-2',
                expected_patterns=(go_patterns + disease_orthology_patterns +
                                   orthology_patterns),
                description_categories=common_categories
            ),
            # Used to study - behavioral phenotypes
            TestGene(
                gene_id='WB:WBGene00001179',
                gene_symbol='egl-10',
                expected_patterns=(go_patterns + disease_used_to_study_patterns +
                                   expression_patterns),
                description_categories=common_categories
            ),
            # Disease via orthology - Notch pathway, cancer research
            TestGene(
                gene_id='WB:WBGene00003001',
                gene_symbol='lin-12',
                expected_patterns=(go_patterns + disease_orthology_patterns +
                                   expression_patterns),
                description_categories=common_categories
            ),
            # Biomarker potential - APP homolog, Alzheimer's research
            TestGene(
                gene_id='WB:WBGene00000149',
                gene_symbol='apl-1',
                expected_patterns=(go_patterns + disease_biomarker_patterns +
                                   disease_used_to_study_patterns),
                description_categories=common_categories
            ),
        ],

        'MGI': [  # Mouse - MGI
            # Used to study cancer - oncogene
            TestGene(
                gene_id='MGI:96680',
                gene_symbol='Kras',
                expected_patterns=(go_patterns + disease_used_to_study_patterns),
                description_categories=common_categories
            ),
            # Biomarker potential - MYC pathway component
            TestGene(
                gene_id='MGI:96921',
                gene_symbol='Max',
                expected_patterns=(go_patterns + disease_biomarker_patterns +
                                   expression_patterns),
                description_categories=common_categories
            ),
            # Used to study neuroblastoma - oncogene
            TestGene(
                gene_id='MGI:97357',
                gene_symbol='Mycn',
                expected_patterns=(go_patterns + disease_used_to_study_patterns +
                                   expression_patterns),
                description_categories=common_categories
            ),
            # Used to study cancer, biomarker - tumor suppressor
            TestGene(
                gene_id='MGI:98834',
                gene_symbol='Trp53',
                expected_patterns=(go_patterns + disease_used_to_study_patterns +
                                   disease_biomarker_patterns),
                description_categories=common_categories
            ),
            # Biomarker for breast cancer - BRCA1
            TestGene(
                gene_id='MGI:104537',
                gene_symbol='Brca1',
                expected_patterns=(go_patterns + disease_biomarker_patterns +
                                   disease_used_to_study_patterns),
                description_categories=common_categories
            ),
        ],

        'SGD': [  # Yeast - SGD
            # Disease via orthology - cytochrome P450, drug metabolism
            TestGene(
                gene_id='SGD:S000002810',
                gene_symbol='DIT2',
                expected_patterns=(
                    go_patterns + disease_orthology_patterns + orthology_patterns
                ),
                description_categories=common_categories
            ),
            # Used to study mitochondrial diseases - ribosomal protein
            TestGene(
                gene_id='SGD:S000001486',
                gene_symbol='MRP17',
                expected_patterns=(
                    go_patterns + disease_used_to_study_patterns + orthology_patterns
                ),
                description_categories=common_categories
            ),
            # Disease via orthology - NAD transporter, metabolic disorders
            TestGene(
                gene_id='SGD:S000001268',
                gene_symbol='YIA6',
                expected_patterns=(
                    go_patterns + disease_orthology_patterns + orthology_patterns
                ),
                description_categories=common_categories
            ),
            # Used to study cancer - DNA repair, BRCA1 ortholog
            TestGene(
                gene_id='SGD:S000000897',
                gene_symbol='RAD51',
                expected_patterns=(
                    go_patterns + disease_used_to_study_patterns + orthology_patterns
                ),
                description_categories=common_categories
            ),
        ],

        'RGD': [  # Rat - RGD (4 genes - diverse pathways)
            # Biomarker for cardiovascular disease - renin
            TestGene(
                gene_id='RGD:3555',
                gene_symbol='Ren',
                expected_patterns=go_patterns + disease_biomarker_patterns,
                description_categories=common_categories
            ),
            # FGFR1 - growth factor receptor, development/disease
            TestGene(
                gene_id='RGD:620713',
                gene_symbol='Fgfr1',
                expected_patterns=(
                    go_patterns + disease_biomarker_patterns + disease_used_to_study_patterns
                ),
                description_categories=common_categories
            ),
            # DNA methylation - epigenetics pathway
            TestGene(
                gene_id='RGD:1303274',
                gene_symbol='Dnmt3b',
                expected_patterns=(
                    go_patterns + disease_used_to_study_patterns + orthology_patterns
                ),
                description_categories=common_categories
            ),
            # Transcription factor - gene regulation
            TestGene(
                gene_id='RGD:11414885',
                gene_symbol='Hoxa1',
                expected_patterns=go_patterns + expression_patterns + orthology_patterns,
                description_categories=common_categories
            ),
        ],

        'ZFIN': [  # Zebrafish - ZFIN (4 genes - diverse pathways)
            # Used to study development and disease - sonic hedgehog
            TestGene(
                gene_id='ZFIN:ZDB-GENE-980526-166',
                gene_symbol='shha',
                expected_patterns=(
                    go_patterns + disease_used_to_study_patterns + expression_patterns
                ),
                description_categories=common_categories
            ),
            # PDX1 - pancreatic development, diabetes research
            TestGene(
                gene_id='ZFIN:ZDB-GENE-990415-122',
                gene_symbol='pdx1',
                expected_patterns=(
                    go_patterns + disease_orthology_patterns + expression_patterns + 
                    orthology_patterns
                ),
                description_categories=common_categories
            ),
            # TCF7L2 - Wnt signaling, diabetes/metabolic disease
            TestGene(
                gene_id='ZFIN:ZDB-GENE-991110-8',
                gene_symbol='tcf7l2',
                expected_patterns=(
                    go_patterns + disease_orthology_patterns + expression_patterns + 
                    orthology_patterns
                ),
                description_categories=common_categories
            ),
            # Notch signaling - developmental patterning
            TestGene(
                gene_id='ZFIN:ZDB-GENE-011128-3',
                gene_symbol='jag2b',
                expected_patterns=go_patterns + expression_patterns + orthology_patterns,
                description_categories=common_categories
            ),
        ],

        'FB': [  # Drosophila - FlyBase (4 genes - diverse pathways)
            # Armadillo/beta-catenin - Wnt signaling, cell adhesion
            TestGene(
                gene_id='FB:FBgn0000117',
                gene_symbol='arm',
                expected_patterns=go_patterns + expression_patterns + orthology_patterns,
                description_categories=common_categories
            ),
            # Delta - Notch signaling ligand, development
            TestGene(
                gene_id='FB:FBgn0000463',
                gene_symbol='Delta',
                expected_patterns=go_patterns + expression_patterns + orthology_patterns,
                description_categories=common_categories
            ),
            # Chromatin regulation - development
            TestGene(
                gene_id='FB:FBgn0032030',
                gene_symbol='Wdr82',
                expected_patterns=go_patterns + expression_patterns + orthology_patterns,
                description_categories=common_categories
            ),
            # Used to study neurodegenerative diseases - alpha-synuclein
            TestGene(
                gene_id='FB:FBgn0283521',
                gene_symbol='Snca',
                expected_patterns=(
                    go_patterns + disease_used_to_study_patterns + expression_patterns
                ),
                description_categories=common_categories
            ),
        ],

        'HUMAN': [  # Human - HGNC
            # Used to study + biomarker cancer - "guardian of the genome"
            TestGene(
                gene_id='RGD:HGNC:11998',
                gene_symbol='TP53',
                expected_patterns=(
                    go_patterns + disease_used_to_study_patterns + disease_biomarker_patterns
                ),
                description_categories=common_categories
            ),
            # Biomarker for breast cancer - BRCA1
            TestGene(
                gene_id='RGD:HGNC:1100',
                gene_symbol='BRCA1',
                expected_patterns=(
                    go_patterns + disease_biomarker_patterns + disease_used_to_study_patterns
                ),
                description_categories=common_categories
            ),
            # Used to study cancer - oncogene
            TestGene(
                gene_id='RGD:HGNC:6407',
                gene_symbol='KRAS',
                expected_patterns=go_patterns + disease_used_to_study_patterns,
                description_categories=common_categories
            ),
            # Used to study cancer + biomarker - transcription factor
            TestGene(
                gene_id='RGD:HGNC:7553',
                gene_symbol='MYC',
                expected_patterns=(
                    go_patterns + disease_used_to_study_patterns + disease_biomarker_patterns
                ),
                description_categories=common_categories
            ),
            # Used to study cancer - growth factor receptor
            TestGene(
                gene_id='RGD:HGNC:3236',
                gene_symbol='EGFR',
                expected_patterns=go_patterns + disease_used_to_study_patterns,
                description_categories=common_categories
            ),
        ],

        'XBXL': [  # Xenopus laevis - XenBase
            # Disease via orthology - TBX1, DiGeorge syndrome
            TestGene(
                gene_id='Xenbase:XB-GENE-478088',
                gene_symbol='tbx1.L',
                expected_patterns=(
                    go_patterns + disease_orthology_patterns + expression_patterns + 
                    orthology_patterns
                ),
                description_categories=common_categories
            ),
            # Used to study development - brachyury/T
            TestGene(
                gene_id='Xenbase:XB-GENE-6255736',
                gene_symbol='tbxt.L',
                expected_patterns=(
                    go_patterns + disease_orthology_patterns + expression_patterns + 
                    orthology_patterns
                ),
                description_categories=common_categories
            ),
            # Disease via orthology - casein kinase, Alzheimer's
            TestGene(
                gene_id='Xenbase:XB-GENE-478125',
                gene_symbol='csnk1a1.L',
                expected_patterns=(
                    go_patterns + disease_orthology_patterns + orthology_patterns
                ),
                description_categories=common_categories
            ),
            # Disease via orthology - hemoglobin alpha
            TestGene(
                gene_id='Xenbase:XB-GENE-865144',
                gene_symbol='hba1.S',
                expected_patterns=(
                    go_patterns + disease_orthology_patterns + expression_patterns + 
                    orthology_patterns
                ),
                description_categories=common_categories
            ),
            # Used to study development - homeobox gene
            TestGene(
                gene_id='Xenbase:XB-GENE-17337706',
                gene_symbol='hoxc6.L',
                expected_patterns=go_patterns + expression_patterns + orthology_patterns,
                description_categories=common_categories
            ),
        ],

        'XBXT': [  # Xenopus tropicalis - XenBase
            # Disease via orthology - TBX1, DiGeorge syndrome
            TestGene(
                gene_id='Xenbase:XB-GENE-478084',
                gene_symbol='tbx1',
                expected_patterns=(
                    go_patterns + disease_orthology_patterns + expression_patterns + 
                    orthology_patterns
                ),
                description_categories=common_categories
            ),
            # Used to study development - brachyury/T
            TestGene(
                gene_id='Xenbase:XB-GENE-478789',
                gene_symbol='tbxt',
                expected_patterns=(
                    go_patterns + disease_orthology_patterns + expression_patterns + 
                    orthology_patterns
                ),
                description_categories=common_categories
            ),
            # Disease via orthology - casein kinase, Alzheimer's
            TestGene(
                gene_id='Xenbase:XB-GENE-478121',
                gene_symbol='csnk1a1',
                expected_patterns=(
                    go_patterns + disease_orthology_patterns + orthology_patterns
                ),
                description_categories=common_categories
            ),
            # Used to study development - sonic hedgehog
            TestGene(
                gene_id='Xenbase:XB-GENE-488039',
                gene_symbol='shh',
                expected_patterns=(
                    go_patterns + disease_orthology_patterns + expression_patterns + 
                    orthology_patterns
                ),
                description_categories=common_categories
            ),
            # Used to study development - master eye regulator
            TestGene(
                gene_id='Xenbase:XB-GENE-484088',
                gene_symbol='pax6',
                expected_patterns=(
                    go_patterns + disease_orthology_patterns + expression_patterns + 
                    orthology_patterns
                ),
                description_categories=common_categories
            ),
        ]
    }


def main():
    """Run all manual tests and coverage threshold tests."""
    print("Manual Gene Description Validation Tests")
    print("=" * 60)
    print("This script tests that well-known genes have appropriate descriptions")
    print("with expected patterns and data categories, plus coverage thresholds.")
    print()

    validator = DescriptionValidator()
    test_genes = define_test_genes()
    coverage_thresholds = define_coverage_thresholds()

    # Run gene-specific tests
    all_results = {}
    total_passed = 0
    total_tests = 0

    print("\n" + "="*80)
    print("PART 1: GENE-SPECIFIC VALIDATION TESTS")
    print("="*80)

    for provider, genes in test_genes.items():
        result = validator.run_tests(provider, genes)
        all_results[provider] = result

        if 'error' not in result:
            total_passed += result['passed']
            total_tests += result['total_genes']

    # Run coverage threshold tests
    coverage_results = {}
    total_coverage_passed = 0
    total_coverage_tests = 0

    print("\n" + "="*80)
    print("PART 2: COVERAGE THRESHOLD TESTS")
    print("="*80)

    for provider, thresholds in coverage_thresholds.items():
        if provider in test_genes:  # Only test providers we have gene tests for
            result = validator.test_coverage_thresholds(provider, thresholds)
            coverage_results[provider] = result

            if 'error' not in result:
                total_coverage_passed += result['passed_thresholds']
                total_coverage_tests += len(result['coverage_tests'])

    # Combined Summary
    print(f"\n{'='*80}")
    print("OVERALL SUMMARY")
    print(f"{'='*80}")

    print("\nGene-Specific Tests:")
    print(f"  Total tests: {total_tests}")
    print(f"  Total passed: {total_passed}")
    print(f"  Success rate: {total_passed/total_tests*100:.1f}%" if total_tests > 0 else "N/A")

    print("\nCoverage Threshold Tests:")
    print(f"  Total threshold tests: {total_coverage_tests}")
    print(f"  Total passed: {total_coverage_passed}")
    print(
        f"  Success rate: {total_coverage_passed/total_coverage_tests*100:.1f}%"
        if total_coverage_tests > 0 else "N/A"
    )

    # Provider breakdown
    print("\nBy provider:")
    for provider in test_genes.keys():
        gene_result = all_results.get(provider, {})
        coverage_result = coverage_results.get(provider, {})

        if 'error' in gene_result:
            print(f"  {provider}: GENE ERROR - {gene_result['error']}")
        elif 'error' in coverage_result:
            print(f"  {provider}: COVERAGE ERROR - {coverage_result['error']}")
        else:
            gene_rate = gene_result.get('passed', 0) / gene_result.get('total_genes', 1) * 100
            coverage_rate = (
                coverage_result.get('passed_thresholds', 0) /
                len(coverage_result.get('coverage_tests', [1])) * 100
            )
            print(
                f"  {provider}: Genes {gene_result.get('passed', 0)}/"
                f"{gene_result.get('total_genes', 0)} ({gene_rate:.1f}%), "
                f"Coverage {coverage_result.get('passed_thresholds', 0)}/"
                f"{len(coverage_result.get('coverage_tests', []))} ({coverage_rate:.1f}%)"
            )

    print(f"\n{'='*80}")
    print("Notes:")
    print("- Gene tests check for flexible patterns that may evolve over time")
    print("- Coverage tests ensure minimum description levels are maintained")
    print("- Thresholds are set 20% below current coverage levels")
    print("- Some categories (expression, disease) vary by MOD data availability")


if __name__ == "__main__":
    main()
