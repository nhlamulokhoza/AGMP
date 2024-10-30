# -*- coding: utf-8 -*-
"""
GWAS catalog data wrangling, step-up-step
AGMP Project
Author: Nhlamulo Khoza (nhlamulokhoza0@gmail.com)
"""

"""# Dependencies"""

!pip install pandas requests tqdm
import requests
import pandas as pd
import tqdm
import string
import os
import numpy as np
import gzip
import shutil
import xml.etree.ElementTree as ET
from time import sleep
from concurrent.futures import ThreadPoolExecutor, as_completed

"""# Data"""

df = pd.read_excel('gwas_africa_assoc_MERGED.xls')

"""# Data Preparation"""

df = df.drop(['gene_name', 'function', 'uniprot'], axis=1, errors='ignore')

# df.columns

"""# Query BioMart for Gene Names

Prepare BioMart input file:
"""

# Count and print unique rsIDs
unique_id_count = df['id'].nunique()
total_id_count = len(df['id'])
print(f"Unique rsIDs: {unique_id_count} / {total_id_count}")

# Write and clean unique rsIDs in a single step
with open('rsID.txt', 'w') as outfile:
    for rsid in df['id'].unique():
        # Remove punctuation, replace 'x', and clean the rsID
        cleaned_rsid = rsid.translate(str.maketrans('', '', string.punctuation)).replace('x', '').strip()
        if cleaned_rsid.startswith('rs'):
            outfile.write(f"{cleaned_rsid}\n")
        else:
            # Handle possible split by 'rs' in case of malformed entries
            for part in cleaned_rsid.split('rs'):
                if part.startswith('rs'):
                    outfile.write(f"rs{part.strip()}\n")

# Additional cleaning step: Remove entries that don't start with 'rs'
with open('rsID.txt', 'r') as infile, open('rsID_cleaned.txt', 'w') as outfile:
    for line in infile:
        for rsid in line.split():
            if rsid.startswith('rs'):
                outfile.write(f"{rsid}\n")

# The file rsID_cleaned.txt has lines in this form "rs145494474rs60366222", read through the lines and if the is rs* in the middle, put rs* and it following numbers in a new line. do this for all rows
with open('rsID_cleaned.txt', 'r') as infile, open('rsID_split.txt', 'w') as outfile:
    for line in infile:
        line = line.strip()
        if 'rs' in line:
            parts = line.split('rs')
            for i in range(1, len(parts)):
                if parts[i]:
                    outfile.write(f"rs{parts[i].strip()}\n")

!mv rsID_split.txt rsID.txt

"""BioMart query:"""

# Reads rsIDs from a given file and returns them as a list
def read_rsids(file_path):
    with open(file_path, 'r') as file:
        rsids = file.read().splitlines()
    return rsids

# Splits the list of rsIDs into chunks of specified size (200) to not overwhelm BioMart server
def chunk_rsids(rsids, chunk_size=200):
    for i in range(0, len(rsids), chunk_size):
        yield rsids[i:i + chunk_size]

# Creates a BioMart XML query with the given rsIDs.
def create_biomart_query(rsids):
    rsids_str = ','.join(rsids)
    query_template = '''<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >
    <Dataset name = "hsapiens_snp" interface = "default" >
        <Filter name = "snp_filter" value = "{rsids}"/>
        <Attribute name = "refsnp_id" />
        <Attribute name = "refsnp_source" />
        <Attribute name = "chr_name" />
        <Attribute name = "chrom_start" />
        <Attribute name = "chrom_end" />
        <Attribute name = "ensembl_gene_name" />
    </Dataset>
</Query>'''
    return query_template.format(rsids=rsids_str)

# Sends the BioMart query for a chunk of rsIDs and returns the results as a pandas DataFrame
def query_biomart(rsids_chunk):
    xml_query = create_biomart_query(rsids_chunk)
    url = 'http://www.ensembl.org/biomart/martservice'
    response = requests.post(url, data={'query': xml_query})
    response.raise_for_status()  # Raise an error for bad status
    data = response.text

    # Convert the TSV response to a DataFrame
    from io import StringIO
    df = pd.read_csv(StringIO(data), sep='\t', header=None)
    df.columns = ["refsnp_id", "refsnp_source", "chr_name", "chrom_start", "chrom_end", "ensembl_gene_name"]
    return df

def main():
    rsid_file_path = 'rsID.txt'
    rsids = read_rsids(rsid_file_path)
    rsid_chunks = list(chunk_rsids(rsids))

    with ThreadPoolExecutor() as executor:
        results = executor.map(query_biomart, rsid_chunks)

    # Combine all the results into a single DataFrame and save to file called 'combined_biomart_results1.csv'
    combined_df = pd.concat(results, ignore_index=True)
    combined_df.to_csv('biomart_results1.csv', index=False)

if __name__ == '__main__':
    main()

"""Process BioMart output:"""

# Keep only unique rows based on 'refsnp_id' and 'ensembl_gene_name'
df2 = pd.read_csv('biomart_results1.csv')
df2 = df2.drop_duplicates(subset=['refsnp_id', 'ensembl_gene_name'])

"""Add Gene Names to GWAS file:"""

# Merge 'combined_df' with 'merged_df' based on 'refsnp_id' and 'id' and drop the extra 'refsnp_id' colummn in the output -- "merged_df"
merged_df = df.merge(df2[['refsnp_id', 'ensembl_gene_name']],
                            left_on='id', right_on='refsnp_id',
                            how='left')
merged_df.drop('refsnp_id', axis=1, inplace=True)

"""# Query BioMart for UniProt IDs and Function

Prepare BioMart input file:
"""

# Extract 'ensembl_gene_name' values from 'merged_df' and write them to a file called 'gene_names.txt' -- this will be the input for BioMart Query/API #2
gene_names = merged_df['ensembl_gene_name'].tolist()

with open('gene_names.txt', 'w') as file:
    for gene_name in gene_names:
        if isinstance(gene_name, str):
            file.write(gene_name.translate(str.maketrans('', '', string.punctuation)) + '\n')

"""BioMart query:"""

# Function to read gene names from file
def read_gene_names(file_path):
    with open(file_path, 'r') as file:
        gene_names = file.read().splitlines()
    return gene_names

# Function to create a BioMart query XML for a batch of gene names
def create_query_xml(gene_names_batch):
    query_xml = '''<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query  virtualSchemaName = "default" formatter = "TSV" header = "1" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >
	<Dataset name = "hsapiens_gene_ensembl" interface = "default" >
		<Filter name = "ensembl_gene_id" value = "{gene_names}"/>
		<Attribute name = "ensembl_gene_id" />
        <Attribute name = "name_1006" />
        <Attribute name = "definition_1006" />
        <Attribute name = "uniprotswissprot" />
	</Dataset>
</Query>'''.format(gene_names=','.join(gene_names_batch))
    return query_xml

# Query BioMart using POST request
def query_biomart(query_xml):
    url = 'https://www.ensembl.org/biomart/martservice'
    response = requests.post(url, data={'query': query_xml})
    if response.status_code == 200:
        return response.text
    else:
        response.raise_for_status()

# Handle gene names in batches (500 gene names at a time) to now overwhelm BioMart server
def handle_batch(gene_names_batch):
    query_xml = create_query_xml(gene_names_batch)
    result = query_biomart(query_xml)
    return result

# Main function to process gene names in batches and combine the results thereafter
def main(file_path, batch_size=500):
    gene_names = read_gene_names(file_path)
    total_batches = len(gene_names) // batch_size + (1 if len(gene_names) % batch_size > 0 else 0)

    results = []
    with ThreadPoolExecutor(max_workers=10) as executor:
        futures = []
        for i in range(0, len(gene_names), batch_size):
            gene_names_batch = gene_names[i:i + batch_size]
            futures.append(executor.submit(handle_batch, gene_names_batch))

        for future in tqdm.tqdm(as_completed(futures), total=total_batches, desc="Processing Batches"):
            results.append(future.result())

    combined_results = ''.join(results)

    with open('biomart_results2.tsv', 'w') as file:
        file.write(combined_results)

if __name__ == "__main__":
    main('gene_names.txt')

"""Process BioMart output:"""

df3 = pd.read_csv('biomart_results2.tsv', sep='\t')       # BioMart Output
df3 = df3.drop_duplicates(subset=['Gene stable ID', 'GO term definition', 'UniProtKB/Swiss-Prot ID'])
df3 = df3.drop('GO term name', axis=1)                    # Delete 'Go term name' column, we dont need it

# Aggregate df3 so that for each "Gene stable ID", only one row is retained that includes all the 'UniProtKB/Swiss-Prot ID' and 'GO term definition' entries matching to that "Gene stable ID". Separate them by '|' but only keep the unique entries
df3_agg = df3.groupby('Gene stable ID').agg({
    'UniProtKB/Swiss-Prot ID': lambda x: '|'.join(set(x.dropna().astype(str))),
    'GO term definition': lambda x: '|'.join(set(x.dropna().astype(str)))
}).reset_index()

# Rename columns
df3_agg = df3_agg.rename(columns={
    'Gene stable ID': 'gene_name',
    'UniProtKB/Swiss-Prot ID': 'uniprot',
    'GO term definition': 'function'})

# Remove rows with "Gene stable ID" in column "gene_name"
df3_agg = df3_agg[~df3_agg['gene_name'].str.contains('Gene stable ID', na=False)]

# Change column name in GWAS data in preparation for merge
merged_df = merged_df.rename(columns={'ensembl_gene_name': 'gene_name'})

"""Add UniProt IDs and Functions to GWAS file"""

# Perform the merge and drop duplicates based on a subset of crucial columns
merged_df = merged_df.merge(df3_agg[['gene_name', 'uniprot', 'function']], on='gene_name', how='left')
merged_df.drop_duplicates(inplace=True)

merged_df.drop_duplicates(subset=['id', 'phenotype', 'gene_name', 'uniprot'], inplace=True)

"""# Regional Classifications as per African Union"""

# Define the libraries/AU clasifications then match country in 'origin_of_participants' column to an AU classification
libraries = {
    'Southern Africa': ['Angola', 'Botswana', 'Eswatini', 'Swaziland', 'Lesotho', 'Malawi', 'Mozambique', 'Namibia', 'South Africa', 'Zambia', 'Zimbabwe'],
    'East Africa': ['Comoros', 'Djibouti', 'Eritrea', 'Ethiopia', 'Kenya', 'Madagascar', 'Mauritius', 'Rwanda', 'Seychelles', 'Somalia', 'South Sudan', 'Sudan', 'Tanzania', 'Uganda'],
    'West Africa': ['Benin', 'Burkina Faso', 'Cape Verde', 'Cabo Verde', 'Ivory Coast', 'Gambia', 'Ghana', 'Guinea', 'Guinea-Bissau', 'Liberia', 'Mali', 'Niger', 'Nigeria', 'Senegal', 'Sierra Leone', 'Togo'],
    'North Africa': ['Algeria', 'Egypt', 'Libya', 'Mauritania', 'Morocco', 'Tunisia'],
    'Central Africa': ['Burundi', 'Cameroon', 'Central African Republic', 'Chad', 'Congo Republic', 'Democratic Republic of the Congo', 'Equatorial Guinea', 'Gabon', 'São Tomé and Príncipe'],
    'Oceania': ['Australia', 'Fiji', 'Kiribati', 'Marshall Islands', 'Micronesia', 'Nauru', 'New Zealand', 'Palau', 'Papua New Guinea', 'Samoa', 'Solomon Islands', 'Tonga', 'Tuvalu', 'Vanuatu'],
    'African-American / Afro-Caribbean': ['Canada', 'United States', 'USA', 'U.S.', 'Mexico', 'Guatemala', 'Belize', 'El Salvador', 'Honduras', 'Nicaragua', 'Costa Rica', 'Panama', 'Antigua and Barbuda', 'Bahamas', 'Barbados', 'Cuba', 'Dominica', 'Dominican Republic', 'Grenada', 'Haiti', 'Jamaica', 'Saint Kitts and Nevis', 'Saint Lucia', 'Saint Vincent and the Grenadines', 'Trinidad and Tobago', 'Guatemala', 'Belize', 'El Salvador', 'Honduras', 'Nicaragua', 'Costa Rica', 'Panama', 'Colombia', 'Venezuela', 'Guyana', 'Suriname', 'Brazil', 'Ecuador', 'Peru', 'Bolivia', 'Chile', 'Argentina', 'Paraguay', 'Uruguay'],
    'Europe': ['Albania', 'Andorra', 'Austria', 'Belarus', 'Belgium', 'Bosnia and Herzegovina', 'Bulgaria', 'Croatia', 'Cyprus', 'Czech Republic (Czechia)', 'Czechia', 'Czech Republic', 'Denmark', 'Estonia', 'Finland', 'France', 'Germany', 'Greece', 'Hungary', 'Iceland', 'Ireland', 'Italy', 'Kosovo', 'Latvia', 'Liechtenstein', 'Lithuania', 'Luxembourg', 'Malta', 'Moldova', 'Monaco', 'Montenegro', 'Netherlands', 'North Macedonia', 'Norway', 'Poland', 'Portugal', 'Romania', 'Russia', 'San Marino', 'Serbia', 'Slovakia', 'Slovenia', 'Spain', 'Sweden', 'Switzerland', 'Ukraine', 'United Kingdom (UK)', 'United Kingdom', 'UK', 'Vatican City'],
    'Asia': ['Afghanistan', 'Armenia', 'Azerbaijan', 'Bahrain', 'Bangladesh', 'Bhutan', 'Brunei', 'Cambodia', 'China', 'Cyprus', 'Georgia', 'India', 'Indonesia', 'Iran', 'Iraq', 'Israel', 'Japan', 'Jordan', 'Kazakhstan', 'Kuwait', 'Kyrgyzstan', 'Laos', 'Lebanon', 'Malaysia', 'Maldives', 'Mongolia', 'Myanmar (Burma)', 'Burma', 'Myanmar', 'Nepal', 'North Korea', 'Oman', 'Pakistan', 'Palestine', 'Philippines', 'Qatar', 'Saudi Arabia', 'Singapore', 'South Korea', 'Sri Lanka', 'Syria', 'Taiwan', 'Tajikistan', 'Thailand', 'Timor-Leste (East Timor)', 'East Timor', 'Timor-Leste', 'Turkey', 'Turkmenistan', 'United Arab Emirates (UAE)', 'UAE', 'United Arab Emirates', 'Uzbekistan', 'Vietnam', 'Yemen'],
    'Africa unspecified':['NR']
}

def find_libraries(countries):
    matched_libraries = []
    for library, country_list in libraries.items():
        for country in countries:
            if country in country_list:
                matched_libraries.append(library)
    return ', '.join(matched_libraries)

merged_df['libraries'] = merged_df['origin_of_participants'].apply(lambda x: find_libraries(x.split(', ')))
merged_df.drop_duplicates(inplace=True)

"""# Additional checks and edits to the GWAS Catalog data"""

# Check if the columns exist before deleting them
if 'geographical_region' in merged_df.columns:
    del merged_df['geographical_region']
if 'function_x' in merged_df.columns:
    del merged_df['function_x']
if 'uniprot_x' in merged_df.columns:
    del merged_df['uniprot_x']

# Rename columns
merged_df = merged_df.rename(columns={'function_y': 'function'})
merged_df = merged_df.rename(columns={'libraries': 'geographical_region'})
merged_df = merged_df.rename(columns={'uniprot_y': 'uniprot'})
merged_df = merged_df.rename(columns={'Pubmed_id': 'PUBMEDID'})

column_order = ['id', 'PUBMEDID', 'phenotype', 'origin_of_participants', 'geographical_region',
                 'Ethnicity', 'mixed_population', 'p-value', 'curated_gene_symbol', 'Notes',
                 'study_type', 'gene_symbol_mined', 'gene_name', 'chromosome', 'uniprot',
                 'function', 'data_ac', 'publication_type', 'title', 'publication',
                 'publication_year', 'source', 'id_in_source', 'variant_type', 'drug_name',
                 'ID Drug bank', 'state', 'Indication', 'IUPAC_name']

merged_df = merged_df.reindex(columns=column_order)
merged_df.drop_duplicates(inplace=True)

# Desplaying enteries matching to multiple AU classification
merged_df['geographical_region'] = merged_df['geographical_region'].apply(lambda x: ', '.join(set(x.split(', '))))
# ...and functions
merged_df['function'] = merged_df['function'].apply(lambda x: '; '.join(set(str(x).split('; '))) if isinstance(x, str) else x)

"""# Resolving Mixed Population Column"""

# Resolve mixed_population column
def determine_mixed_population(row):
  if pd.isnull(row['origin_of_participants']) or row['origin_of_participants'] == 'NR':
    return 'NR'
  elif ',' in row['origin_of_participants'] or ';' in row['origin_of_participants']:
    return 'TRUE'
  else:
    return 'FALSE'

merged_df['mixed_population'] = merged_df.apply(determine_mixed_population, axis=1)

"""# Dropping dupplicates

Dropping duplicates based on 'id', 'gene_name' and 'uniprot' column
"""

merged_df = merged_df.drop_duplicates(subset=['id', 'gene_name', 'uniprot'], keep='first')
merged_df = merged_df.reset_index(drop=True)

"""In the case of duplicated rs IDs, keep only the ones with 'gene_name' entries"""

# Group by 'id' and count the occurrences of each row
counts = merged_df.groupby('id').size().reset_index(name='count')
# Merge the counts back into the original dataframe
merged_df_with_counts = pd.merge(merged_df, counts, on='id')
# Sort the DataFrame by 'id' and 'count' in descending order
merged_df_sorted = merged_df_with_counts.sort_values(['id', 'count'], ascending=[True, False])
# Drop duplicates based on 'id', keeping only the first occurrence (which will be the most populated row)
merged_df_unique = merged_df_sorted.drop_duplicates(subset=['id'], keep='first')
# Remove the 'count' column if it's no longer needed and reset the index
merged_df_unique = merged_df_unique.drop('count', axis=1)
merged_df_unique = merged_df_unique.reset_index(drop=True)

# df_unique contains only unique 'id' entries, with duplicates removed by keeping the most populated rows.
print(f"The new df has {merged_df_unique.shape[0]} rows.")
merged_df = merged_df_unique

"""# Gene Symbols"""

def match_gene_symbol(row):
    # Handle potential float values by converting to string first
    curated_genes = set(str(row['curated_gene_symbol']).replace(',', ' ').split())
    mined_genes = set(str(row['gene_symbol_mined']).replace(',', ' ').split())

    # Find intersection
    matched_genes = curated_genes & mined_genes

    # Return matched gene if found, otherwise return curated_gene_symbol
    if matched_genes:
        return ', '.join(matched_genes)  # Join matched genes with a comma
    else:
        return row['curated_gene_symbol']

# Apply the function to each row in 'merged_df'
merged_df['validate'] = merged_df.apply(match_gene_symbol, axis=1)

# Drop duplicates considering all columns except 'curated_gene_symbol' and 'gene_symbol_mined'
merged_df = merged_df.drop_duplicates(subset=[col for col in merged_df.columns if col not in ['curated_gene_symbol', 'gene_symbol_mined']])

"""# Preparing the dataframe for ingest"""

# Copy values from 'PUBMEDID' to 'publication'
merged_df['publication'] = merged_df['PUBMEDID']

# # Resolve 'PUBMEDID' by copying the values of "publication" and unpopulate "Ethnicity" therefater
# merged_df['PUBMEDID'] = merged_df.groupby('id')['publication'].transform('first')

# Resolve 'study_type' column
merged_df['study_type'] = 'GWAS'

# Remove rows with commas or '-' or '--'  or NaN anywhere in the 'validate' column and reset the index
merged_df = merged_df[~merged_df['validate'].str.contains(',|-|--', na=False)]
merged_df = merged_df.dropna(subset=['validate'])
merged_df = merged_df.reset_index(drop=True)

# Delete the 'curated_gene_symbol' column and Rename the 'validate' column to 'curated_gene_symbol'
if 'curated_gene_symbol' in merged_df.columns:
  merged_df = merged_df.drop('curated_gene_symbol', axis=1)

merged_df = merged_df.rename(columns={'validate': 'curated_gene_symbol'})

column_order = ['id', 'PUBMEDID', 'phenotype', 'origin_of_participants', 'geographical_region',
                 'Ethnicity', 'mixed_population', 'p-value', 'curated_gene_symbol', 'Notes',
                 'study_type', 'gene_symbol_mined', 'gene_name', 'chromosome', 'uniprot',
                 'function', 'data_ac', 'publication_type', 'title', 'publication',
                 'publication_year', 'source', 'id_in_source', 'variant_type', 'drug_name',
                 'ID Drug bank', 'state', 'Indication', 'IUPAC_name']

merged_df = merged_df.reindex(columns=column_order)
merged_df.drop_duplicates(inplace=True)

# print(f"Rows with NaN gene_names: {merged_df['gene_name'].isnull().sum()}")
# print(f"Rows without NaN in gene_name: {merged_df.dropna(subset=['gene_name']).shape[0]}")

# Remove rows with no Gene Name enteries
new_df = merged_df.dropna(subset=['gene_name']).reset_index(drop=True)

# Replace '\' with '-' in phenotype (this is a requirement for the portal search functionality)
new_df['phenotype'] = new_df['phenotype'].str.replace(r'[\\/]', '-', regex=True)

new_df['curated_gene_symbol'] = new_df['curated_gene_symbol'].replace('intergenic', 'Intergenic', regex=True)

# Save the DataFrame to Excel and CSV files
new_df.to_excel('gwas_catalog_data.xlsx', index=False)
new_df.to_csv('gwas_catalog_data.csv', index=False)
