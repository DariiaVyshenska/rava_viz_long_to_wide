import pandas as pd

def extract_mut_type(df):
  df = df.rename(columns={'Syn': 'MUT_TYPE'})
  mut_type_df = df[['POSITION:NT_CHANGE', 'MUT_TYPE']].drop_duplicates()
  
  # Group and collapse mutation types where more than one type across samples
  groupped = (
    mut_type_df
    .groupby('POSITION:NT_CHANGE')['MUT_TYPE']
    .agg(lambda x: ' or '.join(sorted(set(x))))
    .reset_index()
  )
  return groupped

def extract_mat_peptide(df):
  mat_pep = df[df['MatPeptide'] != '-']
  if mat_pep.empty:
    return None
  
  mat_pep = (
    mat_pep[['POSITION:NT_CHANGE', 'MatPeptide']]
    .drop_duplicates()
    .reset_index(drop=True)
  )
  mat_pep[
    ['MAT_PEPTIDE', 'MAT_PEP_AA_CHANGE', 'MAT_PEP_NT_CHANGE']
  ] = mat_pep['MatPeptide'].str.split(
    r"[:;] ", regex=True, expand = True
  )
  mat_pep.drop(columns=['MatPeptide'], inplace=True)
  return mat_pep

def extract_snv_metadata(df):
  meta = df[
    ["POSITION:NT_CHANGE", "NucCorrect", "AminoCorrect", "Position", "Protein"]
  ].drop_duplicates()

  meta.rename(columns={
    'NucCorrect': 'NT_CHANGE',
    'AminoCorrect': 'AA_CHANGE',
    'Position': 'POSITION',
    'Protein': 'PROTEIN',
  }, inplace=True)

  # extract mutation type information
  mut_types = extract_mut_type(df[['POSITION:NT_CHANGE', 'Syn']])
  meta = pd.merge(meta, mut_types, on='POSITION:NT_CHANGE', how='outer')

  # extract MatPeptide information
  mat_peptide_info = extract_mat_peptide(df[['POSITION:NT_CHANGE', 'MatPeptide']])
  if mat_peptide_info is not None:
    meta = pd.merge(meta, mat_peptide_info, on='POSITION:NT_CHANGE', how='left')

  return meta