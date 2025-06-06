import pandas as pd
import os
import argparse
from extractors import extract_snv_metadata

def long_to_wide(input_csv, headers, output_dir):
  df = pd.read_csv(input_csv)
  df['Position'] = df['Position'].astype(int)
  headers_df = pd.read_csv(headers)
  headers_df = headers_df[~headers_df['NEW_HEADER'].isna()]

  df['Sample'] = df['Sample'].str.replace('.fastq.gz', '')  # NEXT: refactor: extract into a separate function (getting wide table)
  df = df[df['Sample'].isin(headers_df['SAMPLE_ID'])]
  df = pd.merge(df, headers_df, left_on='Sample', right_on='SAMPLE_ID', how='left')

  df['Syn'] = df['Syn'].str.replace(' SNV', '')
  df['POSITION:NT_CHANGE'] = df['Position'].astype(str) + ':' + df['NucCorrect']
  af_df_long = df[['NEW_HEADER', 'POSITION:NT_CHANGE', 'AF']]
  af_df_wide = af_df_long.pivot(index='POSITION:NT_CHANGE', columns='NEW_HEADER', values='AF')
  
  # reordering columns based on headers_df
  af_df_wide = af_df_wide.reindex(columns=headers_df['NEW_HEADER'].tolist())
  os.makedirs(output_dir, exist_ok=True)
  af_df_wide.to_csv(os.path.join(output_dir, 'tmp.csv'), index=False)
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
  af_df_wide.to_csv(os.path.join(output_dir, 'af_df_wide.csv'), index=False)

def main():
  parser = argparse.ArgumentParser(description='Transform long RAVA visualization.csv into wide format',
                                usage='Usage: ./main.py <input_csv> <headers_csv> [--output_dir <output_dir>]',
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