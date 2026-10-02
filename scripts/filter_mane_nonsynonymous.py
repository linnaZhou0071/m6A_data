import pandas as pd
import gzip
import re
import argparse
import sys

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
    parser.add_argument('--input-csv', required=True)
    parser.add_argument('--output-csv', required=True)
    args = parser.parse_args()
    
    mane_set = parse_mane_transcripts(args.mane_gtf)
    
    print(f"Loading input CSV: {args.input_csv}...")
    try:
        df = pd.read_csv(args.input_csv, dtype=str)
    except Exception as e:
        print(f"Error reading CSV: {e}")
        sys.exit(1)
        
    print(f"Total rows in input: {len(df)}")
    
    if 'effect' not in df.columns:
        print("Error: 'effect' column not found.")
        sys.exit(1)
        
    nonsyn_df = df[df['effect'] == 'nonsynonymous'].copy()
    print(f"Rows with 'nonsynonymous' effect: {len(nonsyn_df)}")
    
    nonsyn_df['clean_transcript_id'] = nonsyn_df['analysis_transcript'].apply(clean_transcript_id)
    nonsyn_df['is_MANE'] = nonsyn_df['clean_transcript_id'].isin(mane_set)
    
    mane_nonsyn_df = nonsyn_df[nonsyn_df['is_MANE']].copy()
    print(f"Rows matching MANE Select transcripts: {len(mane_nonsyn_df)}")
    
    mane_nonsyn_df.to_csv(args.output_csv, index=False)
    print(f"Saved filtered results to: {args.output_csv}")
    
    print("\n--- Verification Sample ---")
    if not mane_nonsyn_df.empty:
        sample = mane_nonsyn_df.head(1)
        tid = sample.iloc[0]['analysis_transcript']
        clean_tid = sample.iloc[0]['clean_transcript_id']
        print(f"Sample Row Transcript: {tid} (Cleaned: {clean_tid})")
        print(f"Is '{clean_tid}' in MANE set? {'Yes' if clean_tid in mane_set else 'No'}")

if __name__ == "__main__":
    main()
