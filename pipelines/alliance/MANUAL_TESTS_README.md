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
- These genes should consistently have rich descriptions
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

Genes are selected to cover all three disease annotation categories:
- **Biomarker**: Genes that serve as indicators or markers for diseases
- **Used to Study**: Genes used as research models for human diseases
- **Disease via Orthology**: Genes linked to disease through human orthologs

Some genes may have multiple categories, providing comprehensive coverage.

### C. elegans (WB)
- **daf-2**: Insulin receptor (disease via orthology - aging/longevity)
- **egl-10**: RGS protein (used to study - behavioral phenotypes)
- **egl-15**: FGFR (used to study - development)
- **egl-19**: Calcium channel (disease via orthology - epilepsy)
- **lin-12**: Notch receptor (disease via orthology - cancer)
- **sma-1**: Spectrin (used to study - developmental disorders)
- **apl-1**: APP homolog (biomarker - Alzheimer's disease)

### Mouse (MGI)
- **Kras**: Oncogene (used to study - cancer)
- **Max**: MYC binding partner (biomarker - MYC pathway)
- **Mycn**: Oncogene (used to study - neuroblastoma)
- **Trp53**: Tumor suppressor (used to study + biomarker - cancer)
- **Brca1**: Tumor suppressor (biomarker - breast cancer)

### Yeast (SGD)
- **DIT2**: Cytochrome P450 (disease via orthology - drug metabolism)
- **MRP17**: Mitochondrial ribosome (used to study - mitochondrial diseases)
- **YIA6**: NAD transporter (disease via orthology - metabolic disorders)
- **RAD51**: DNA repair (used to study - cancer, BRCA1 ortholog)

### Rat (RGD) 
- **Tp53**: Tumor suppressor (used to study + biomarker - cancer)
- **Ren**: Renin (biomarker - cardiovascular disease)

### Zebrafish (ZFIN)
- **tp53**: Tumor suppressor (used to study - cancer)
- **shha**: Sonic hedgehog (used to study - development and disease)

### Drosophila (FB)
- **p53**: Tumor suppressor (used to study - cancer)
- **Snca**: Alpha-synuclein (used to study - neurodegenerative diseases)

### Human (HUMAN)
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

The tests look for these flexible patterns in descriptions:

### Gene Ontology Patterns
- `"Predicted to enable"` - Molecular function predictions
- `"Involved in"` - Biological process involvement
- `"Located in"` - Cellular component location
- `"Part of"` - Component membership
- `"Acts upstream"` - Pathway relationships

### Expression Patterns
- `"Is expressed in"` - Tissue/cell expression
- `"Expressed in"` - Alternative expression phrasing

### Disease Patterns

#### Biomarker Patterns
- `"biomarker.*for"` - Biomarker for disease
- `"marker.*for"` - Marker for condition
- `"indicator.*of"` - Indicator of disease
- `"associated.*with.*disease"` - Disease association

#### Used to Study Patterns
- `"Used to study"` - Research applications
- `"Model.*for"` - Disease model
- `"Research.*model"` - Research model
- `"Study.*of"` - Study applications

#### Disease via Orthology Patterns
- `"Human ortholog.*implicated in"` - Disease via orthology
- `"ortholog.*associated.*with"` - Ortholog disease association
- `"Orthologous to human.*gene.*implicated"` - Orthology-based disease link
- `"implicated.*via.*orthology"` - Explicit orthology connection

### Orthology Patterns
- `"Orthologous to"` - General orthology
- `"Human ortholog"` - Specific human orthology

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