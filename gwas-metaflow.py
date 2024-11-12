# Author: Nhlamulo Khoza
# Email: nhlamulokhoza0@gmail.com
# Developed for the Developed for the African Genomic Medicine Portal (AGMP), part of the African Genomics Data Hub.

from metaflow import FlowSpec, step, conda_base, batch, resources, retry
import pandas as pd
import requests
import string
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm
import numpy as np
from io import StringIO

@conda_base(libraries={'pandas': '2.0.3', 'requests': '2.31.0', 'tqdm': '4.65.0'})
class GWASProcessingFlow(FlowSpec):
    
    @step
    def start(self):
        """
        Start the flow and load initial GWAS data
        """
        self.df = pd.read_excel('gwas_africa_assoc_MERGED.xls')
        self.df = self.df.drop(['gene_name', 'function', 'uniprot'], axis=1, errors='ignore')
        self.next(self.prepare_rsids)

    @step
    def prepare_rsids(self):
        """
        Prepare and clean rsIDs for BioMart query
        """
        unique_id_count = self.df['id'].nunique()
        total_id_count = len(self.df['id'])
        print(f"Unique rsIDs: {unique_id_count} / {total_id_count}")

        # Clean and write rsIDs
        with open('rsID.txt', 'w') as outfile:
            for rsid in self.df['id'].unique():
                cleaned_rsid = rsid.translate(str.maketrans('', '', string.punctuation)).replace('x', '').strip()
                if cleaned_rsid.startswith('rs'):
                    outfile.write(f"{cleaned_rsid}\n")
                else:
                    for part in cleaned_rsid.split('rs'):
                        if part.startswith('rs'):
                            outfile.write(f"rs{part.strip()}\n")

        # Additional cleaning for split rsIDs
        with open('rsID.txt', 'r') as infile, open('rsID_split.txt', 'w') as outfile:
            for line in infile:
                line = line.strip()
                if 'rs' in line:
                    parts = line.split('rs')
                    for i in range(1, len(parts)):
                        if parts[i]:
                            outfile.write(f"rs{parts[i].strip()}\n")

        # Read final cleaned rsIDs
        with open('rsID_split.txt', 'r') as file:
            self.rsids = file.read().splitlines()
            
        self.next(self.query_biomart_genes)

    @retry(times=3)
    @step
    def query_biomart_genes(self):
        """
        Query BioMart for gene information using rsIDs
        """
        def create_biomart_query(rsids):
            query_template = '''<?xml version="1.0" encoding="UTF-8"?>
            <!DOCTYPE Query>
            <Query virtualSchemaName="default" formatter="TSV" header="0" uniqueRows="0" count="" datasetConfigVersion="0.6">
                <Dataset name="hsapiens_snp" interface="default">
                    <Filter name="snp_filter" value="{rsids}"/>
                    <Attribute name="refsnp_id" />
                    <Attribute name="refsnp_source" />
                    <Attribute name="chr_name" />
                    <Attribute name="chrom_start" />
                    <Attribute name="chrom_end" />
                    <Attribute name="ensembl_gene_name" />
                </Dataset>
            </Query>'''
            return query_template.format(rsids=','.join(rsids))

        def query_chunk(rsids_chunk):
            xml_query = create_biomart_query(rsids_chunk)
            response = requests.post('http://www.ensembl.org/biomart/martservice', 
                                  data={'query': xml_query})
            response.raise_for_status()
            df = pd.read_csv(StringIO(response.text), sep='\t', header=None)
            df.columns = ["refsnp_id", "refsnp_source", "chr_name", "chrom_start", 
                         "chrom_end", "ensembl_gene_name"]
            return df

        # Process in chunks
        chunk_size = 200
        chunks = [self.rsids[i:i+chunk_size] for i in range(0, len(self.rsids), chunk_size)]
        
        results = []
        with ThreadPoolExecutor() as executor:
            futures = [executor.submit(query_chunk, chunk) for chunk in chunks]
            for future in tqdm(as_completed(futures), total=len(chunks)):
                results.append(future.result())

        self.biomart_df = pd.concat(results, ignore_index=True)
        self.biomart_df = self.biomart_df.drop_duplicates(subset=['refsnp_id', 'ensembl_gene_name'])
        
        self.next(self.query_biomart_uniprot)

    @retry(times=3)
    @step
    def query_biomart_uniprot(self):
        """
        Query BioMart for UniProt IDs and function information
        """
        # Merge previous results with GWAS data
        merged_df = self.df.merge(self.biomart_df[['refsnp_id', 'ensembl_gene_name']], 
                                left_on='id', right_on='refsnp_id', how='left')
        merged_df.drop('refsnp_id', axis=1, inplace=True)
        
        # Extract gene names for UniProt query
        gene_names = merged_df['ensembl_gene_name'].dropna().unique().tolist()

        def create_uniprot_query(gene_names_batch):
            query_xml = '''<?xml version="1.0" encoding="UTF-8"?>
            <!DOCTYPE Query>
            <Query virtualSchemaName="default" formatter="TSV" header="1" uniqueRows="0" count="" datasetConfigVersion="0.6">
                <Dataset name="hsapiens_gene_ensembl" interface="default">
                    <Filter name="ensembl_gene_id" value="{gene_names}"/>
                    <Attribute name="ensembl_gene_id" />
                    <Attribute name="name_1006" />
                    <Attribute name="definition_1006" />
                    <Attribute name="uniprotswissprot" />
                </Dataset>
            </Query>'''
            return query_xml.format(gene_names=','.join(gene_names_batch))

        def query_batch(gene_names_batch):
            query_xml = create_uniprot_query(gene_names_batch)
            response = requests.post('https://www.ensembl.org/biomart/martservice', 
                                  data={'query': query_xml})
            return response.text if response.status_code == 200 else None

        # Process in batches
        batch_size = 500
        batches = [gene_names[i:i+batch_size] for i in range(0, len(gene_names), batch_size)]
        
        results = []
        with ThreadPoolExecutor(max_workers=10) as executor:
            futures = [executor.submit(query_batch, batch) for batch in batches]
            for future in tqdm(as_completed(futures), total=len(batches)):
                result = future.result()
                if result:
                    results.append(result)

        # Process results
        combined_results = ''.join(results)
        self.uniprot_df = pd.read_csv(StringIO(combined_results), sep='\t')
        
        self.merged_df = merged_df
        self.next(self.process_results)

    @step
    def process_results(self):
        """
        Process and combine all results
        """
        # Process UniProt results
        uniprot_df = self.uniprot_df.drop_duplicates(subset=['Gene stable ID', 'GO term definition', 'UniProtKB/Swiss-Prot ID'])
        uniprot_df = uniprot_df.drop('GO term name', axis=1)
        
        # Aggregate UniProt data
        uniprot_agg = uniprot_df.groupby('Gene stable ID').agg({
            'UniProtKB/Swiss-Prot ID': lambda x: '|'.join(set(x.dropna().astype(str))),
            'GO term definition': lambda x: '|'.join(set(x.dropna().astype(str)))
        }).reset_index()
        
        # Rename columns
        uniprot_agg = uniprot_agg.rename(columns={
            'Gene stable ID': 'gene_name',
            'UniProtKB/Swiss-Prot ID': 'uniprot',
            'GO term definition': 'function'
        })
        
        # Clean up and merge
        uniprot_agg = uniprot_agg[~uniprot_agg['gene_name'].str.contains('Gene stable ID', na=False)]
        self.merged_df = self.merged_df.rename(columns={'ensembl_gene_name': 'gene_name'})
        self.merged_df = self.merged_df.merge(uniprot_agg, on='gene_name', how='left')
        
        self.next(self.process_geographic_data)

    @step
    def process_geographic_data(self):
        """
        Process geographic classifications and mixed populations
        """
        # Define regional classifications
        libraries = {
            'Southern Africa': ['Angola', 'Botswana', 'Eswatini', 'Swaziland', 'Lesotho', 'Malawi', 
                              'Mozambique', 'Namibia', 'South Africa', 'Zambia', 'Zimbabwe'],
            'East Africa': ['Comoros', 'Djibouti', 'Eritrea', 'Ethiopia', 'Kenya', 'Madagascar', 
                          'Mauritius', 'Rwanda', 'Seychelles', 'Somalia', 'South Sudan', 'Sudan', 
                          'Tanzania', 'Uganda'],
            # ... [other regions as defined in your original code]
            'Africa unspecified': ['NR']
        }

        def find_libraries(countries):
            matched_libraries = []
            for library, country_list in libraries.items():
                for country in countries:
                    if country in country_list:
                        matched_libraries.append(library)
            return ', '.join(matched_libraries)

        def determine_mixed_population(row):
            if pd.isnull(row['origin_of_participants']) or row['origin_of_participants'] == 'NR':
                return 'NR'
            elif ',' in row['origin_of_participants'] or ';' in row['origin_of_participants']:
                return 'TRUE'
            return 'FALSE'

        # Apply geographic classifications
        self.merged_df['libraries'] = self.merged_df['origin_of_participants'].apply(
            lambda x: find_libraries(x.split(', ')))
        self.merged_df['mixed_population'] = self.merged_df.apply(determine_mixed_population, axis=1)
        
        self.next(self.finalise_data)

    @step
    def finalise_data(self):
        """
        Finalise data processing and save results
        """
        # Final column renaming and ordering
        self.merged_df = self.merged_df.rename(columns={
            'libraries': 'geographical_region',
            'Pubmed_id': 'PUBMEDID'
        })

        # Process gene symbols
        def match_gene_symbol(row):
            curated_genes = set(str(row['curated_gene_symbol']).replace(',', ' ').split())
            mined_genes = set(str(row['gene_symbol_mined']).replace(',', ' ').split())
            matched_genes = curated_genes & mined_genes
            return ', '.join(matched_genes) if matched_genes else row['curated_gene_symbol']

        self.merged_df['validate'] = self.merged_df.apply(match_gene_symbol, axis=1)
        
        # Final cleaning and formatting
        self.merged_df = self.merged_df[~self.merged_df['validate'].str.contains(',|-|--', na=False)]
        self.merged_df = self.merged_df.dropna(subset=['validate'])
        self.merged_df = self.merged_df.rename(columns={'validate': 'curated_gene_symbol'})
        
        # Clean phenotype field
        self.merged_df['phenotype'] = self.merged_df['phenotype'].str.replace(r'[\\/]', '-', regex=True)
        self.merged_df['curated_gene_symbol'] = self.merged_df['curated_gene_symbol'].replace(
            'intergenic', 'Intergenic', regex=True)
        
        # Save results
        self.merged_df.to_excel('gwas_catalog_data.xlsx', index=False)
        self.merged_df.to_csv('gwas_catalog_data.csv', index=False)
        
        self.next(self.end)

    @step
    def end(self):
        """
        End the flow
        """
        print("Flow completed successfully!")

if __name__ == '__main__':
    GWASProcessingFlow()
