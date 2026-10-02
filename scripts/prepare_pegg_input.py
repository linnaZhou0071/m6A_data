import pandas as pd
import re
import sys
import os
import argparse

def parse_pos(pos_str):
    # Expected format: chr10_100151706_-
    if not isinstance(pos_str, str):
        return None, None, None
    m = re.fullmatch(r'chr([0-9XY]+)_(\d+)_([+-])', pos_str)
    if not m:
        return None, None, None
    return m.group(1), int(m.group(2)), m.group(3)

def convert_to_pegg_format(input_file, output_file):
    print(f"Processing {input_file}...")
    if not os.path.exists(input_file):
        print(f"File not found: {input_file}")
        raise FileNotFoundError(input_file)

    try:
        df = pd.read_csv(input_file)
    except Exception as e:
        print(f"Error reading {input_file}: {e}")
        raise RuntimeError(f"Cannot read {input_file}") from e

    pegg_rows = []
    skipped_positions = []
    
    for idx, row in df.iterrows():
        chrom, pos, strand = parse_pos(row['pos'])
        if chrom is None:
            skipped_positions.append(row['pos'])
            continue
            
        # Logic for A->G mutation
        # If +, Ref=A, Alt=G
        # If -, Ref=T, Alt=C
        
        if strand == '+':
            ref = 'A'
            alt = 'G'
        else:
            ref = 'T'
            alt = 'C'
            
        pegg_row = {
            'Chromosome': chrom,
            'Start_Position': pos,
            'End_Position': pos,
            'Variant_Type': 'SNP',
            'Reference_Allele': ref,
            'Tumor_Seq_Allele2': alt,
            'Original_Pos_ID': row['pos']
        }
        pegg_rows.append(pegg_row)
        
    pegg_df = pd.DataFrame(pegg_rows)
    
    # Save to CSV
    print(f"Saving {len(pegg_df)} rows to {output_file}...")
    pegg_df.to_csv(output_file, index=False)
    if skipped_positions:
        print(f"Skipped {len(skipped_positions)} positions with unsupported format: {', '.join(map(str, skipped_positions))}", file=sys.stderr)

def main():
    parser = argparse.ArgumentParser(description="Convert filtered m6A sites to PEGG CSV inputs")
    parser.add_argument("--input-dir", default="data/intermediate")
    parser.add_argument("--output-dir", default="data/intermediate")
    parser.add_argument("--input-file", help="Explicit filtered CSV; use with --output-file")
    parser.add_argument("--output-file", help="Explicit PEGG CSV; use with --input-file")
    args = parser.parse_args()
    if bool(args.input_file) != bool(args.output_file):
        parser.error("--input-file and --output-file must be supplied together")
    if args.input_file:
        convert_to_pegg_format(args.input_file, args.output_file)
        return
    os.makedirs(args.output_dir, exist_ok=True)
    files = [
        (os.path.join(args.input_dir, "HEK_HeLa_over80_intersection_filtered.csv"),
         os.path.join(args.output_dir, "pegg_input_intersection.csv")),
        (os.path.join(args.input_dir, "HEK_HeLa_over80_union_filtered.csv"),
         os.path.join(args.output_dir, "pegg_input_union.csv")),
    ]
    
    for infile, outfile in files:
        convert_to_pegg_format(infile, outfile)

if __name__ == "__main__":
    main()
