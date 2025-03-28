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