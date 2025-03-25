import pandas as pd
import os
import argparse

def extract_mut_type(df):
  df = df.rename(columns={'Syn': 'MUT_TYPE'})
  mut_type_df = df[['POSITION:NT_CHANGE', 'MUT_TYPE']].drop_duplicates()
  mut_types_count = mut_type_df['POSITION:NT_CHANGE'].value_counts().reset_index()
  mut_types_count.columns = ['POSITION:NT_CHANGE', 'Freq']
  mut_type_count_ambig = mut_types_count[mut_types_count['Freq'] > 1]
  ambig_mut_all = mut_type_df[mut_type_df['POSITION:NT_CHANGE'].isin(mut_type_count_ambig['POSITION:NT_CHANGE'])]
  
  ambig_mut_meta = pd.DataFrame({
    'POSITION:NT_CHANGE': ambig_mut_all['POSITION:NT_CHANGE'].unique(),
    'MUT_TYPE': None
  })

  for i in range(len(ambig_mut_meta)):
    mut_id = ambig_mut_meta.loc[i, 'POSITION:NT_CHANGE']
    mut_types = ambig_mut_all.loc[ambig_mut_all['POSITION:NT_CHANGE'] == mut_id, 'MUT_TYPE'].sort_values().unique()
    new_mut_type = ' or '.join(mut_types)
    ambig_mut_meta.loc[i, 'MUT_TYPE'] = new_mut_type

  unambig_mut = mut_types_count[mut_types_count['Freq'] == 1]
  unambig_mut_meta = mut_type_df[mut_type_df['POSITION:NT_CHANGE'].isin(unambig_mut['POSITION:NT_CHANGE'])]
  unambig_mut
  mut_types_all = pd.concat([unambig_mut_meta, ambig_mut_meta], ignore_index=True)
  return mut_types_all

def extract_mat_peptide(df):
  mat_pep = df[df['MatPeptide'] != '-']
  if mat_pep.empty:
    return None
  
  mat_pep = mat_pep[['POSITION:NT_CHANGE', 'MatPeptide']].drop_duplicates().reset_index(drop=True)
  mat_pep[['MAT_PEPTIDE', 'MAT_PEP_AA_CHANGE', 'MAT_PEP_NT_CHANGE']] = mat_pep['MatPeptide'].str.split(r"[:;] ", regex=True, expand = True) #.tolist()
  mat_pep.drop(columns=['MatPeptide'], inplace=True)
  return mat_pep


def extract_snv_metadata(df):
  meta = df[["POSITION:NT_CHANGE", "NucCorrect", "AminoCorrect", "Position", "Protein"]]

  meta = meta.rename(columns={
    'NucCorrect': 'NT_CHANGE',
    'AminoCorrect': 'AA_CHANGE',
    'Position': 'POSITION',
    'Protein': 'PROTEIN',
  })
  meta.drop_duplicates(inplace=True)
  
  # extracting mutation type information
  mut_types = extract_mut_type(df[['POSITION:NT_CHANGE', 'Syn']])
  meta = pd.merge(meta, mut_types, on='POSITION:NT_CHANGE', how='outer')

  # parsing MatPeptide information
  mat_peptide_info = extract_mat_peptide(df[['POSITION:NT_CHANGE', 'MatPeptide']])
  if mat_peptide_info is not None:
    meta = pd.merge(meta, mat_peptide_info, on='POSITION:NT_CHANGE', how='left')

  return meta
  


def long_to_wide(input_csv, headers, output_dir):
  df = pd.read_csv(input_csv)
  df['Position'] = df['Position'].astype(int)
  headers_df = pd.read_csv(headers)

  df['Sample'] = df['Sample'].str.replace('.fastq.gz', '')  # refactor: extract into a separate function (getting wide table)
  df = df[df['Sample'].isin(headers_df['SAMPLE_ID'])]
  df = pd.merge(df, headers_df, left_on='Sample', right_on='SAMPLE_ID', how='left')

  df['Syn'] = df['Syn'].str.replace(' SNV', '')
  df['POSITION:NT_CHANGE'] = df['Position'].astype(str) + ':' + df['NucCorrect']
  af_df_long = df[['NEW_HEADER', 'POSITION:NT_CHANGE', 'AF']]
  af_df_wide = af_df_long.pivot(index='POSITION:NT_CHANGE', columns='NEW_HEADER', values='AF')
  
  # reordering columns based on headers_df
  af_df_wide = af_df_wide.reindex(columns=headers_df['NEW_HEADER'].tolist())
  
  max_vals = af_df_wide.max(axis=1, skipna=True)
  count_vals = af_df_wide.count(axis=1)
  af_df_wide['MAX_AF'] = max_vals
  af_df_wide['COUNT_AF'] = count_vals   # end of getting wide table function
  
  # extracting metadata for each SNV:
  snv_metadata = extract_snv_metadata(df)

  # merge AF wide dataframe with SNV metadata
  af_df_wide = pd.merge(af_df_wide, snv_metadata, on='POSITION:NT_CHANGE', how='outer')
  
  # Save the wide format DataFrame to a CSV file
  output_file_path = os.path.join(output_dir, 'af_df_wide.xlsx')   # refactor: extract into a separate function
  try:
    os.makedirs(output_dir, exist_ok=True)
  except OSError as e:
    raise ValueError(f"Failed to create output directory: {output_dir}. Error: {e}")
  af_df_wide.to_excel(output_file_path, index=False, engine='openpyxl')

def main():
  parser = argparse.ArgumentParser(description='Transform long RAVA visualization.csv into wide format',
                                usage='Usage: ./main.py <input_csv> <headers_csv> --output_dir=<output_dir>',
                                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('input_csv', type=str, help="Path to RAVA's visualization.csv file")
  parser.add_argument(
    'headers', 
    type=str, 
    help=('Path to CSV file with column "SAMPLE_ID" AND "NEW_HEADER". '
          'Row order will be used to order wide columns'))
  parser.add_argument('--output_dir', type=str, default='./results', help='Path to output directory')

  args = parser.parse_args()
  print("Extracting data...")
  long_to_wide(args.input_csv, args.headers, args.output_dir)
  print("Done!")
  
if __name__ == "__main__":
  main()