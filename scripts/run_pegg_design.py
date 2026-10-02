import pandas as pd
from pegg import prime
import os
import sys
import argparse

def main():
    parser = argparse.ArgumentParser(description="Run PEGG design")
    parser.add_argument("--genome", default="GCF_000001405.26_GRCh38_genomic.fna.gz", help="Path to genome fasta.gz")
    parser.add_argument("--input", required=True, help="Input CSV file (prepared for PEGG)")
    parser.add_argument("--output", required=True, help="Output CSV file for results")
    args = parser.parse_args()

    # --- Configuration ---
    genome_filename = args.genome
    input_csv = args.input
    output_csv = args.output
    
    # Custom Parameters defined by user
    rtt_lengths_param = [20, 22, 25, 27, 30, 32, 35]
    min_rha_size_param = 6
    use_sensor_param = False
    # ---------------------

    if not os.path.exists(genome_filename):
        print(f"Error: Genome file '{genome_filename}' not found.")
        sys.exit(1)
        
    if not os.path.exists(input_csv):
        print(f"Error: Input CSV '{input_csv}' not found.")
        sys.exit(1)

    print(f"Loading genome from {genome_filename}...")
    try:
        # prime.genome_loader returns (chrom_dict, val)
        chrom_dict, val = prime.genome_loader(genome_filename)
        print("Genome loaded successfully.")
    except Exception as e:
        print(f"Error loading genome: {e}")
        sys.exit(1)
        
    print(f"Loading mutations from {input_csv}...")
    mutant_input = pd.read_csv(input_csv)
    print(f"Loaded {len(mutant_input)} variants.")
    
    # --- Verification: Print Parameters ---
    print("\n" + "="*40)
    print("  PEGG RUN PARAMETERS CHECK")
    print("="*40)
    print(f"  Format       : cBioPortal")
    print(f"  Sensor       : {use_sensor_param}")
    print(f"  RTT Lengths  : {rtt_lengths_param}")
    print(f"  Min RHA Size : {min_rha_size_param}")
    print("="*40 + "\n")

    print("Running PEGG design... (This may take a while)")
    try:
        results = prime.run(
            mutant_input, 
            input_format='cBioPortal', 
            chrom_dict=chrom_dict,
            sensor=use_sensor_param,
            RTT_lengths=rtt_lengths_param,
            min_RHA_size=min_rha_size_param
        )
        
        print(f"Design complete. Generated {len(results)} pegRNAs.")
        
        print(f"Saving results to {output_csv}...")
        results.to_csv(output_csv, index=False)
        print("Done.")
        
    except Exception as e:
        print(f"Error during PEGG execution: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()