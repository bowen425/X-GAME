import pandas as pd

# File paths
input_file_path = "/X-GAME/data/merged_inter_data.xlsx"
output_file_path = "/X-GAME/data/merged_inter_data_with_context.xlsx"

# Read the data
try:
    data = pd.read_excel(input_file_path)
except Exception as e:
    print(f"Error reading the input file: {e}")
    raise

# Step 1: Find all distinct cross-linking sites for each PPI
ppi_sites = data.groupby(['Protein1', 'Protein2']).apply(
    lambda x: {
        'protein1_sites': set(x['Site1']),  # All distinct sites of Protein1
        'protein2_sites': set(x['Site2']),  # All distinct sites of Protein2
        'total_sites': len(x)               # Total number of cross-linking site pairs
    }
).reset_index(name='sites_info')

# Step 2: Identify Rich PPIs (Criteria: At least two cross-linking pairs, and each protein has at least two distinct sites)
rich_ppis = []
for _, row in ppi_sites.iterrows():
    protein1, protein2 = row['Protein1'], row['Protein2']
    sites_info = row['sites_info']
    
    # Check if the criteria are met:
    # 1. At least two cross-linking pairs
    # 2. Protein1 has at least two distinct sites
    # 3. Protein2 has at least two distinct sites
    if (sites_info['total_sites'] >= 2 and 
        len(sites_info['protein1_sites']) >= 2 and 
        len(sites_info['protein2_sites']) >= 2):
        rich_ppis.append((protein1, protein2))

# Convert rich_ppis to a DataFrame
rich_ppis_df = pd.DataFrame(rich_ppis, columns=['Protein1', 'Protein2'])

# Step 3: Add Context column to the original data
data['Context'] = 'poor'  # Default value is "poor"
for idx, row in data.iterrows():
    protein1, protein2 = row['Protein1'], row['Protein2']
    # Check if this PPI is in rich_ppis
    if ((rich_ppis_df['Protein1'] == protein1) & (rich_ppis_df['Protein2'] == protein2)).any():
        data.at[idx, 'Context'] = 'rich'

# Step 4: Save the results to a new Excel file
try:
    data.to_excel(output_file_path, index=False)
    print(f"Processed data saved to: {output_file_path}")
except Exception as e:
    print(f"Error saving the output file: {e}")
    raise
