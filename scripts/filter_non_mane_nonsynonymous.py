import pandas as pd
import gzip
import re
import argparse
import sys
import os

def parse_mane_transcripts(gtf_path):
    mane_transcripts = set()
    print(f"Parsing MANE GTF: {gtf_path}...")
    opener = gzip.open if str(gtf_path).endswith('.gz') else open
    try:
        with opener(gtf_path, 'rt', encoding='utf-8', errors='replace') as f:
            for line in f:
                if line.startswith('#'): continue
                if 'tag "MANE Select"' not in line: continue
                m = re.search(r'transcript_id "([^"]+)"', line)
                if m: mane_transcripts.add(m.group(1))
    except Exception as e:
        print(f"Error reading GTF: {e}")
        sys.exit(1)
    print(f"Found {len(mane_transcripts)} unique MANE Select transcripts.")
    return mane_transcripts

def clean_transcript_id(tid):
    if pd.isna(tid): return ""
    tid = str(tid)
    m = re.match(r'^(.*\.\d+)(_\d+)$', tid)
    if m: return m.group(1)
    return tid

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--mane-gtf', required=True)
    parser.add_argument('--effects-csv', required=True)
    parser.add_argument('--original-csv', required=True)
    parser.add_argument('--output-csv', required=True)
    args = parser.parse_args()
    
    mane_set = parse_mane_transcripts(args.mane_gtf)
    
    print(f"Loading effects CSV: {args.effects_csv}...")
    try:
        eff_df = pd.read_csv(args.effects_csv, dtype=str)
    except Exception as e:
        print(f"Error reading effects CSV: {e}")
        sys.exit(1)

    eff_df['clean_transcript_id'] = eff_df['analysis_transcript'].apply(clean_transcript_id)
    eff_df['is_MANE'] = eff_df['clean_transcript_id'].isin(mane_set)
    
    mane_nonsyn_rows = eff_df[
        (eff_df['is_MANE']) & 
        (eff_df['effect'] == 'nonsynonymous')
    ]
    
    mane_nonsyn_sites = set(mane_nonsyn_rows['pos'].unique())
    print(f"Found {len(mane_nonsyn_sites)} unique sites that are Nonsynonymous in MANE Select transcripts.")
    
    print(f"Loading original CSV: {args.original_csv}...")
    try:
        orig_df = pd.read_csv(args.original_csv, dtype=str)
    except Exception as e:
        print(f"Error reading original CSV: {e}")
        sys.exit(1)
        
    print(f"Total sites in original file: {len(orig_df)}")
    
    filtered_df = orig_df[~orig_df['pos'].isin(mane_nonsyn_sites)].copy()
    
    print(f"Sites remaining after filtering: {len(filtered_df)}")
    print(f"Removed {len(orig_df) - len(filtered_df)} sites.")
    
    print(f"Saving to CSV: {args.output_csv}...")
    filtered_df.to_csv(args.output_csv, index=False)
    print("Done.")
    
    print("\n--- Verification Sample ---")
    if len(mane_nonsyn_sites) > 0:
        removed_site = list(mane_nonsyn_sites)[0]
        print(f"Site '{removed_site}' was removed.")
        print(f"Reason: It is Nonsynonymous in MANE transcript(s).")
        evidence = eff_df[
            (eff_df['pos'] == removed_site) & 
            (eff_df['is_MANE']) & 
            (eff_df['effect'] == 'nonsynonymous')
        ][['analysis_transcript', 'clean_transcript_id', 'effect']]
        print(evidence.head(1).to_string(index=False))

if __name__ == "__main__":
    main()
