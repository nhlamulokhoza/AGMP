# GWAS Data Processing Pipeline

## Overview
This Metaflow-based pipeline processes GWAS (Genome-Wide Association Studies) data by integrating information from multiple sources, including BioMart queries for gene annotations, UniProt data, and geographic classifications. The pipeline handles data cleaning, transformation, and enrichment to produce a comprehensive dataset suitable for further analysis.

## Features
- Automated processing of GWAS catalog data
- BioMart integration for gene information retrieval
- UniProt data integration
- Geographic region classification based on African Union standards
- Concurrent processing for API queries
- Automatic retry mechanism for failed API calls
- Progress tracking and logging
- Data validation and cleaning
- Export to multiple formats (CSV and Excel)

## Prerequisites
- Python 3.7 or higher
- pip (Python package installer)

### Required Python Packages
```bash
pip install metaflow==2.9.1
pip install pandas==2.0.3
pip install requests==2.31.0
pip install tqdm==4.65.0
pip install openpyxl  # for Excel file handling
```

## Input Requirements
The pipeline expects the following input file:
- `gwas_africa_assoc_MERGED.xls`: Initial GWAS data file containing association data

## Installation
1. Clone this repository:
```bash
git clone <repository-url>
cd gwas-processing-pipeline
```

2. Install required packages:
```bash
pip install -r requirements.txt
```

## Usage
Run the pipeline using the following command:
```bash
python gwas_processing_flow.py run
```

### Pipeline Steps
1. **Start**: Loads and initializes GWAS data
2. **Prepare rsIDs**: Cleans and formats rsID data
3. **Query BioMart Genes**: Retrieves gene information from BioMart
4. **Query BioMart UniProt**: Retrieves UniProt and function information
5. **Process Results**: Combines and processes query results
6. **Process Geographic Data**: Handles regional classifications
7. **Finalise Data**: Performs final data cleaning and saves results

### Output Files
The pipeline generates two output files:
- `gwas_catalog_data.xlsx`: Processed data in Excel format
- `gwas_catalog_data.csv`: Processed data in CSV format

## Output Data Structure
The final dataset includes the following key columns:
- `id`: SNP identifier
- `PUBMEDID`: PubMed reference
- `phenotype`: Associated phenotype
- `geographical_region`: African Union regional classification
- `mixed_population`: Population mixture status
- `curated_gene_symbol`: Validated gene symbols
- `gene_name`: Ensembl gene names
- `uniprot`: UniProt identifiers
- `function`: Gene function descriptions

## Troubleshooting
Common issues and solutions:

1. **BioMart Connection Issues**
```python
# The pipeline automatically retries failed API calls
# You can modify retry parameters in the decorator:
@retry(times=3)
```

2. **Memory Issues**
- Increase available memory or
- Reduce chunk_size in data processing steps

3. **Missing Data**
- Check input file format
- Verify BioMart service status
- Ensure internet connectivity

## Development
To modify the pipeline:

1. Edit step functions in `gwas_processing_flow.py`
2. Add new steps using the `@step` decorator
3. Maintain data flow between steps using self.next()

### Adding New Features
Example of adding a new processing step:
```python
@step
def new_step(self):
    """
    Documentation for new step
    """
    # Process data
    self.next(self.following_step)
```

## Contributing
1. Fork the repository
2. Create a feature branch
3. Commit changes
4. Push to the branch
5. Create a Pull Request

## Performance Considerations
- The pipeline uses concurrent processing for API calls
- Batch sizes are optimized for BioMart queries
- Progress tracking is implemented for long-running operations

## Limitations
- Dependent on BioMart service availability
- Rate limited by external API constraints
- Memory usage scales with input data size

## Contact
Nhlamulo Khoza
Email: nhlamulokhoza0@gmail.com

## Acknowledgments
- Africa Genomics Data Hub
    - Developed for the African Genomic Medicine Portal (AGMP)
- GWAS Catalog for providing phenotype associations
- BioMart for providing genomic data services
- Metaflow developers for the workflow framework
- African Union for regional classification standards
