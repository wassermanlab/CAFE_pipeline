import pandas as pd
import argparse
import os
from io import StringIO

pd.set_option('future.no_silent_downcasting', True)

freq_types =         {
            'variant': str,
            'af_tot': float, 
            'af_xx': float, 
            'af_xy': float,
            'ac_tot': int,
            'ac_xx': int,
            'ac_xy': int,
            'an_tot': int,
            'an_xx': int,
            'an_xy': int,
            'hom_tot': int,
            'hom_xx': int,
            'hom_xy': int,
            'quality': float
        }
num_processed = 0
same_count = 0
diffmax_count = 0
num_triplets = 0
num_quad = 0
num_5groups = 0

def process_tsv(input_file,output_file):
    global same_count, num_processed, diffmax_count, num_triplets, num_quad, num_5groups
    print("processing file", input_file)
    df = pd.read_csv(input_file, sep='\t', dtype=freq_types, na_values=['.'])
#    df = pd.read_csv(input_file, sep='\t')
#    df.replace('.', 0.00000, inplace=True)
    df = df.astype(freq_types)
    
    duplicates = df[df.duplicated(subset=["variant"], keep=False)]
    df_duplicates = df[df["variant"].isin(duplicates["variant"])]
    
    columns_to_sum = [col for col in df_duplicates.columns if col not in ["quality", "an_tot","an_xx","an_xy"] and pd.api.types.is_numeric_dtype(df_duplicates[col])]
    grouped = df_duplicates.groupby("variant")
    
    new_rows = []
    
    for group_name, group in grouped:
#        print("group_name", group_name)
#        print("group quality", group[["quality"]])
#        print("group numeric cols", group[columns_to_sum])
#        print("idx with highest ac_tot", group["ac_tot"].idxmax())
#        print("group quality with highest ac_tot", group.loc[group["ac_tot"].idxmax()]["quality"])
#        print("summed", group[columns_to_sum].sum())
        if len(group) == 3:
            num_triplets += 1
        elif len(group) == 4:
            num_quad += 1
        elif len(group) == 5:
            num_5groups += 1
        
        new_row = {}
        new_row["variant"] = group_name
        new_row.update(group[columns_to_sum].sum())
        an_tot_max = group["an_tot"].idxmax()
        an_xx_max = group["an_xx"].idxmax()
        an_xy_max = group["an_xy"].idxmax()
        if an_tot_max == an_xx_max and an_xx_max == an_xy_max:
            pass
        else:
            diffmax_count += 1
            print("different max an_tot, an_xx, an_xy in group", group)
        
        new_row["an_tot"] = group.loc[an_tot_max]["an_tot"]
        new_row["an_xx"] = group.loc[an_xx_max]["an_xx"]
        new_row["an_xy"] = group.loc[an_xy_max]["an_xy"]
            
        new_row["quality"] = group.loc[group["ac_tot"].idxmax()]["quality"]
        if group["ac_tot"].nunique() == 1:
#            print("identical ac tot in group", group)
            new_row["quality"] = group["quality"].max()
            same_count += 1
        new_rows.append(new_row)
        num_processed += 1
        try:
            test_df = pd.DataFrame([new_row]).astype(freq_types)
        except:
            print("new row fails to enter", new_row)
            raise
#        new_df = pd.concat([new_df,new_row], ignore_index=True)
 
    if (len(new_rows) > 0):
        
        try:
            new_df = pd.DataFrame(new_rows).astype(freq_types) 
        except:
            print("unable to make dataframe from", new_rows)
            raise
#        print("new new", new_df)
    
        new_df.to_csv(output_file, sep='\t', index=False)
        print("done:", output_file)
    else:
        print("no duplicates found in", input_file)
def main():
    parser = argparse.ArgumentParser(description='Produce consolidated duplicated rows from TSV file.')
    parser.add_argument('--input_dir', type=str, help='Path to the dir containing input TSV files', default='.')
    parser.add_argument('--output_dir', type=str, help='Path to the output dir', default='./out')

    args = parser.parse_args()
    files_in_dir = os.listdir(args.input_dir)
    for file in files_in_dir:
        if file.endswith(".tsv") and not file.startswith("0_indelpatch_"):
            input_file = os.path.join(args.input_dir, file)
            output_file = os.path.join(args.output_dir, f"0_indelpatch_{file}")
            process_tsv(input_file, output_file)
    print("done creating patch TSVs. Number of rows:", num_processed)
    print("this is how many dupe indels had same allele count:", same_count)
    print("this is how many dupe indels had max allele numbers from differing copies:", diffmax_count)
    print("this is how many triplets:", num_triplets)
    print("this is how many quadruplets:", num_quad)
    print("this is how many quintuplets:", num_5groups)


if __name__ == '__main__':
    main()
