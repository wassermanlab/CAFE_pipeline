import pandas as pd
import argparse

def process_tsv(input_file, output_file):
    df = pd.read_csv(input_file, sep='\t')
    column_name="variant"
    sum_exceptions = ["quality"]

    duplicates = df[df.duplicated(column_name, keep=False)]
    df_duplicates = df[df[column_name].isin(duplicates[column_name])]
    columns_to_sum = [col for col in df_duplicates.columns if col not in sum_exceptions and pd.api.types.is_numeric_dtype(df_duplicates[col])]
    summed_df = df_duplicates.groupby(column_name)[columns_to_sum].sum().reset_index()


    def custom_quality_aggregation(series):
        # Stub function: picks the greater quality number
        return series.max()

    summed_df['quality'] = df_duplicates.groupby(column_name)['quality'].apply(custom_quality_aggregation).reset_index(drop=True)
    summed_df.to_csv(output_file, sep='\t', index=False)

def main():
    parser = argparse.ArgumentParser(description='Produce consolidated duplicated rows from TSV file.')
    parser.add_argument('input_file', type=str, help='Path to the input TSV file')
    parser.add_argument('output_file', type=str, help='Path to the output TSV file')

    args = parser.parse_args()

    process_tsv(args.input_file, args.output_file)

if __name__ == '__main__':
    main()
