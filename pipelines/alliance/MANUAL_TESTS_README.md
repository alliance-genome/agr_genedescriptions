# Manual Tests for Generated Gene Descriptions

This document describes the manual validation tests for generated gene descriptions across different data providers.

## Overview

The `manual_test_descriptions.py` script validates that well-known, well-studied genes have appropriate descriptions with expected data categories and patterns. These tests are designed to be **manual** and are not part of the automated unit test suite.

## Purpose

- **Quality Assurance**: Ensure that important genes have meaningful descriptions
- **Regression Detection**: Catch when descriptions unexpectedly lose content
- **Pattern Validation**: Verify that expected description patterns are present
- **Data Category Coverage**: Check that key description categories have content

## Design Principles

### Flexibility Over Exactness
- Tests use **flexible pattern matching** rather than exact string comparisons
- Account for data evolution over time (e.g., "Predicted" → "Experimental")
- Allow for minor wording changes in descriptions
- Focus on presence of data categories rather than specific content

### Well-Known Gene Selection
- Selected genes are **famous and well-studied** in their respective organisms
- **Comprehensive annotation requirement**: Good candidate genes have annotations across **all (or most) major categories** including GO, expression, orthology, and disease data
- These genes should consistently have rich descriptions across multiple data types
- If these genes lack descriptions, it indicates a systemic issue

### Data Category Focus
- Tests check for presence of key description categories:
  - `go_description` - Gene Ontology annotations
  - `go_function_description` - Molecular function
  - `go_process_description` - Biological process
  - `go_component_description` - Cellular component
  - `do_description` - Direct disease ontology annotations
  - `do_biomarker_description` - Disease biomarker annotations
  - `do_orthology_description` - Disease via orthology annotations
  - `orthology_description` - Orthology information
  - `tissue_expression_description` - Expression data

## Test Gene Selection

Genes are selected to ensure comprehensive coverage across **all major data categories** with **balanced representation** across MODs:

### Primary Selection Criteria:
1. **Gene Ontology (GO)**: Well-annotated genes with molecular function, biological process, and cellular component data
2. **Disease Associations**: Genes with disease annotations across three categories (where available by MOD):
   - **Biomarker**: Genes that serve as indicators or markers for diseases
   - **Used to Study**: Genes used as research models for human diseases
   - **Disease via Orthology**: Genes linked to disease through human orthologs
3. **Expression Data**: Genes with tissue/anatomical expression information
4. **Orthology**: Genes with well-established cross-species orthologous relationships

### Balanced Selection Strategy:
- **4-5 genes per MOD**: Ensures fair representation across all data providers
- **Pathway diversity**: Avoids over-representation of single gene families (e.g., TP53 orthologs)
- **MOD-appropriate selection**: Accounts for data availability differences (e.g., some MODs don't capture biomarker data)
- **Biological relevance**: Focus on genes famous in their respective research domains

### Test Structure:
The test suite now includes **two complementary validation approaches**:
1. **Gene-Specific Tests**: Validate that well-known genes have rich, appropriate descriptions
2. **Coverage Threshold Tests**: Ensure minimum description coverage levels are maintained across all categories

### C. elegans (WB) - 4 genes
- **daf-2**: Insulin receptor (disease via orthology - aging/longevity)
- **egl-10**: RGS protein (used to study - behavioral phenotypes)
- **lin-12**: Notch receptor (disease via orthology - cancer)
- **apl-1**: APP homolog (biomarker - Alzheimer's disease)

### Mouse (MGI) - 5 genes
- **Kras**: Oncogene (used to study - cancer)
- **Max**: MYC binding partner (biomarker - MYC pathway)
- **Mycn**: Oncogene (used to study - neuroblastoma)
- **Trp53**: Tumor suppressor (used to study + biomarker - cancer)
- **Brca1**: Tumor suppressor (biomarker - breast cancer)

### Yeast (SGD) - 4 genes
- **DIT2**: Cytochrome P450 (disease via orthology - drug metabolism)
- **MRP17**: Mitochondrial ribosome (used to study - mitochondrial diseases)
- **YIA6**: NAD transporter (disease via orthology - metabolic disorders)
- **RAD51**: DNA repair (used to study - cancer, BRCA1 ortholog)

### Rat (RGD) - 4 genes (diverse pathways)
- **Ren**: Renin (biomarker - cardiovascular disease)
- **Fgfr1**: Growth factor receptor (development/disease research)
- **Dnmt3b**: DNA methyltransferase (epigenetics pathway)
- **Hoxa1**: Homeobox transcription factor (gene regulation/development)

### Zebrafish (ZFIN) - 4 genes (diverse pathways)
- **shha**: Sonic hedgehog (used to study - development and disease)
- **pdx1**: Pancreatic development (diabetes research)
- **tcf7l2**: Wnt signaling (diabetes/metabolic disease)
- **jag2b**: Notch signaling (developmental patterning)

### Drosophila (FB) - 4 genes (diverse pathways)
- **arm**: Armadillo/β-catenin (Wnt signaling, cell adhesion)
- **Delta**: Notch signaling ligand (development, cell fate)
- **Wdr82**: Chromatin regulation (development)
- **Snca**: Alpha-synuclein (used to study - neurodegenerative diseases)

### Human (HUMAN) - 5 genes
- **TP53**: "Guardian of the genome" tumor suppressor (used to study + biomarker - cancer)
- **BRCA1**: DNA repair gene (biomarker + used to study - breast/ovarian cancer)
- **KRAS**: Major oncogene (used to study - cancer, RAS signaling)
- **MYC**: Master transcriptional regulator (used to study + biomarker - cancer)
- **EGFR**: Growth factor receptor (used to study - cancer, drug target)

### Xenopus laevis (XBXL)
- **tbx1.L**: T-box transcription factor (disease via orthology - DiGeorge syndrome)
- **tbxt.L**: Brachyury/T (used to study - mesoderm formation, axis specification)
- **csnk1a1.L**: Casein kinase (disease via orthology - Alzheimer's disease)
- **hba1.S**: Hemoglobin alpha (disease via orthology - thalassemia, anemia)
- **hoxc6.L**: Homeobox gene (used to study - development, body patterning)

### Xenopus tropicalis (XBXT)
- **tbx1**: T-box transcription factor (disease via orthology - DiGeorge syndrome)
- **tbxt**: Brachyury/T (used to study - mesoderm formation, axis specification)
- **csnk1a1**: Casein kinase (disease via orthology - Alzheimer's disease)
- **shh**: Sonic hedgehog (used to study - development, morphogenesis)
- **pax6**: Paired box 6 (used to study - eye development, neural development)

## Expected Patterns

The tests look for these patterns in descriptions, extracted directly from the `config_alliance.yml` configuration file to ensure 100% consistency with the actual sentence generation logic:

### Gene Ontology Patterns (from config)
- `"exhibits"` - Function exhibition
- `"enables"` - Function enablement
- `"contributes to"` - Functional contribution
- `"contributes as a"` / `"contributes as an"` - Specific contribution types
- `"colocalizes with"` / `"colocalizes with a"` / `"colocalizes with an"` - Colocalization
- `"predicted to have"` - Predicted functions
- `"predicted to be a"` / `"predicted to be an"` - Predicted types
- `"predicted to enable"` - Predicted function enablement
- `"predicted to contribute to"` - Predicted contributions
- `"predicted to contribute as a"` / `"predicted to contribute as an"` - Predicted specific contributions
- `"predicted to colocalize with"` - Predicted colocalization
- `"involved in"` - Process involvement
- `"predicted to be involved in"` - Predicted process involvement
- `"acts upstream of with a positive effect on"` - Positive upstream regulation
- `"acts upstream of with a negative effect on"` - Negative upstream regulation
- `"acts upstream of or within"` - General upstream or within regulation
- `"acts upstream of or within with a negative effect on"` - Negative upstream or within regulation
- `"acts upstream of or within with a positive effect on"` - Positive upstream or within regulation
- `"predicted to act upstream of..."` - Predicted regulatory relationships
- `"localizes to"` - Cellular localization
- `"located in"` - Cellular location
- `"part of"` - Component membership
- `"is active in"` - Activity location
- `"predicted to localize to"` - Predicted localization
- `"predicted to be located in"` - Predicted location
- `"predicted to be part of"` - Predicted component membership
- `"predicted to be active in"` - Predicted activity location

### Expression Patterns (from config)
- `"is expressed in"` - Tissue/cell expression
- `"is enriched in"` - Expression enrichment

### Disease Patterns (from config)

#### Biomarker Patterns
- `"biomarker of"` - Disease biomarker

#### Used to Study Patterns
- `"used to study"` - Research model applications

#### Disease via Orthology Patterns
- `"human ortholog(s) of this gene implicated in"` - Disease via human orthology
- `"ortholog(s) of this gene implicated in"` - Disease via orthology (for humans)

#### Human-Specific Disease Patterns
- `"implicated in"` - Direct disease implication (for human genes)

### Orthology Patterns (from config)
- `"human ortholog(s) of this gene implicated in"` - Human orthology with disease link
- `"ortholog(s) of this gene implicated in"` - General orthology with disease link

## Coverage Threshold Tests

In addition to gene-specific validation, the test suite includes **automated coverage threshold tests** that ensure minimum description coverage levels are maintained across all data categories.

### Threshold Methodology
- **Baseline calculation**: non-null description counts for each category
- **Threshold setting**: 20% below initially calculated coverage levels to allow for normal fluctuation
- **Categories tested**: Description, GO, Disease, Expression, Orthology coverage
- **MOD-specific**: Accounts for data availability differences (e.g., some MODs don't have expression data)

### Coverage Validation
The tests validate that each MOD maintains:
- **Overall descriptions**: Genes with non-null main descriptions
- **GO annotations**: Gene Ontology coverage (function, process, component)
- **Disease annotations**: Disease ontology coverage (direct, biomarker, orthology-based)
- **Expression data**: Tissue/anatomical expression information
- **Orthology data**: Cross-species orthologous relationships

### Threshold Examples
- **WB**: 24.9% descriptions, 23.5% GO, 5.9% disease, 9.6% expression, 12.0% orthology
- **ZFIN**: 48.1% descriptions, 43.8% GO, 15.0% disease, 21.9% expression, 38.7% orthology
- **XBXT**: 75.8% descriptions, 71.1% GO, 21.8% disease, 2.1% expression, 60.5% orthology

## Pattern Updates (2025)

**Important**: The test patterns have been updated to match exactly the prefixes defined in `config_alliance.yml`. This ensures:

1. **100% consistency** with actual sentence generation logic
2. **No false negatives** due to pattern mismatches
3. **Easier maintenance** - patterns sync with config changes
4. **Accurate validation** - tests what the system actually produces

The patterns now exclude very short prefixes (like "a", "an", "is") that are too generic for meaningful testing, while including all substantive prefixes that indicate specific biological relationships or evidence levels.

## Running the Tests

### Basic Usage
```bash
cd /path/to/agr_genedescriptions
python3 manual_test_descriptions.py
```

### Requirements
- Python 3.6+
- JSON description files in `generated_descriptions/` directory
- No additional dependencies required

### Output Format

**Gene-Specific Test Output:**
```
Testing WB descriptions
============================================================
✅ PASS daf-2 (WB:WBGene00000898)
    Description: Enables PTB domain binding activity...
    Found patterns: Involved in, Located in, Acts upstream...
    Found categories: go_description, do_description...

❌ FAIL some-gene (WB:WBGene00000123)
    Description: No description available
    Missing patterns: Predicted to enable, Involved in...
    Missing categories: go_description, do_description...
```

**Coverage Threshold Test Output:**
```
Testing WB coverage thresholds
============================================================
✅ PASS Description: 30.8% (threshold: 24.9%)
✅ PASS GO: 29.4% (threshold: 23.5%)
✅ PASS Disease: 7.4% (threshold: 5.9%)
✅ PASS Expression: 12.0% (threshold: 9.6%)
✅ PASS Orthology: 15.0% (threshold: 12.0%)

Threshold tests: 5/5 passed (100.0%)
```

## Test Success Criteria

A gene **passes** if:
1. It has a description (not "No description available"), AND
2. Either:
   - Some expected patterns are found in the description, OR
   - Some expected description categories have content

## Interpreting Results

### Expected Outcomes
- **High success rate** (>80%) indicates healthy description generation
- **Low success rate** may indicate:
  - Data pipeline issues
  - Missing or incomplete source data
  - Configuration problems

### Common Failure Modes
1. **"No description available"**: Gene has no generated description
2. **Missing patterns**: Description exists but lacks expected content types
3. **Missing categories**: Individual description categories are empty

### Not Necessarily Problems
- **Missing specific patterns**: Other patterns may be present
- **Empty categories**: Some genes may not have all data types
- **Prediction vs Experimental**: Annotations may evolve over time

## Maintenance

### Adding New Test Genes
To add genes for testing:

1. Find well-known genes in your organism of interest
2. Look up their correct gene IDs in the generated descriptions
3. Add them to the `define_test_genes()` function
4. Choose appropriate expected patterns based on the gene's known biology

### Updating Patterns
- Add new patterns as description formats evolve
- Remove obsolete patterns that are no longer used
- Keep patterns flexible using regex when possible

### Handling Data Changes
- **Predicted → Experimental**: Update expected patterns as predictions become confirmed
- **New annotation sources**: Add patterns for new data types
- **Terminology changes**: Update patterns to match current vocabulary

## File Structure

```
agr_genedescriptions/pipelines/alliance/
├── manual_test_descriptions.py    # Main test script
├── MANUAL_TESTS_README.md         # This documentation
└── generated_descriptions/        # Description files to test
    ├── WB.json                   # C. elegans descriptions
    ├── MGI.json                  # Mouse descriptions
    ├── SGD.json                  # Yeast descriptions
    ├── RGD.json                  # Rat descriptions
    ├── ZFIN.json                 # Zebrafish descriptions
    ├── FB.json                   # Drosophila descriptions
    ├── HUMAN.json                # Human descriptions
    ├── XBXL.json                 # Xenopus laevis descriptions
    ├── XBXT.json                 # Xenopus tropicalis descriptions
    └── ...                       # Other data providers
```

## Future Enhancements

### Potential Additions
- **Coverage metrics**: Track what percentage of genes have descriptions
- **Category completeness**: Monitor which description categories are most/least populated across all data types
- **Historical tracking**: Compare description quality over time
- **Provider-specific patterns**: Customize expected patterns per data provider
- **Automated alerts**: Flag significant drops in description quality
- **Cross-category validation**: Ensure genes have balanced coverage across GO, disease, expression, and orthology

### Integration Possibilities
- **CI/CD integration**: Run tests after description generation
- **Quality dashboards**: Visualize description quality metrics
- **Data provider feedback**: Report quality issues back to source databases
- **Category completeness metrics**: Track which data types are most/least populated per provider
