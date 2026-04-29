#!/usr/bin/env python3
"""
Flag haplotypes with excessively high counts in per-sample haplotype counts table.

This script has two commands:
1. flag - Flag haplotypes with excessively high counts from a single sample
2. merge - Combine flagged haplotypes from multiple samples into a merged list

The flag command loads a per-sample haplotype counts table output by codonwise_tally_vars_hapaware.py,
identifies haplotypes with counts per million above a threshold, and outputs the flagged
haplotype IDs and summary statistics.

The merge command combines flagged haplotypes from multiple samples, removing sample-specific
columns like total_counts and counts_per_million while preserving other haplotype information.
"""

import sys
import argparse
import pandas as pd
import numpy as np
import os


def flag_command(args):
    """Flag haplotypes with excessively high counts from a single sample"""
    # Load the haplotype counts table
    try:
        df = pd.read_csv(args.input_counts, sep='\t')
    except Exception as e:
        sys.stderr.write(f"Error loading input file {args.input_counts}: {e}\n")
        sys.exit(1)

    # Validate required columns
    required_cols = ['hap', 'total_counts']
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        sys.stderr.write(f"Error: Missing required columns: {missing_cols}\n")
        sys.stderr.write(f"Available columns: {list(df.columns)}\n")
        sys.exit(1)

    # Check if counted_codonwise column exists
    has_counted_codonwise = 'counted_codonwise' in df.columns
    if not has_counted_codonwise and not args.use_all_haps:
        sys.stderr.write("Warning: 'counted_codonwise' column not found. Using all haplotypes for sum_counts calculation.\n")
        args.use_all_haps = True


    if args.use_all_haps:
        sum_counts = df['total_counts'].sum()

        if not args.include_wildtype:
            wt_counts = df.loc[df['status']=='wildtype', 'total_counts'].sum()
            sum_counts = sum_counts - wt_counts
    else:
        sum_counts = df.loc[ df['counted_codonwise'] == True, 'total_counts'].sum()
        if args.include_wildtype:
            raise Error('must specify --use_all_haps when using --include_wildtype since it is not normally counted')

    # Calculate sum_counts based on filtering option
    if args.use_all_haps:
        usable_df = df
    else:
        usable_df = df.loc[df['counted_codonwise'] == True]

    usable_df = usable_df.loc[ usable_df['status'] != 'wildtype' ].copy()
    
    # Calculate counts per million for each haplotype
    if sum_counts > 0:
        usable_df['counts_per_million'] = (usable_df['total_counts'] / sum_counts) * 1e6
    else:
        usable_df['counts_per_million'] = 0.0

    # Identify flagged haplotypes
    flagged_mask = usable_df['counts_per_million'] > args.threshold
    flagged_df = usable_df[ (flagged_mask==True) ]

    # Write flagged haplotypes with all columns to output file
    try:
        flagged_df.to_csv(args.output_flagged, sep='\t', index=False)
    except Exception as e:
        sys.stderr.write(f"Error writing flagged haplotypes to {args.output_flagged}: {e}\n")
        sys.exit(1)

    # Calculate summary statistics
    num_flagged_haps = len(flagged_df)
    flagged_read_counts = flagged_df['total_counts'].sum()
    flagged_fraction_of_usable = flagged_read_counts / usable_df['total_counts'].sum() if usable_df['total_counts'].sum() > 0 else 0.0

    flagged_fraction_of_total = flagged_read_counts / df['total_counts'].sum() if df['total_counts'].sum() > 0 else 0.0
    
    # Create summary table
    summary_data = {
        'metric': [
            'libname',
            'num_flagged_haps',
            'flagged_read_counts', 
            'flagged_fraction_of_usable_reads',
            'flagged_fraction_of_total_reads',
            'total_usable_reads',
            'total_reads',
            'threshold_counts_per_million'
        ],
        'value': [
            args.libname,
            num_flagged_haps,
            flagged_read_counts,
            flagged_fraction_of_usable,
            flagged_fraction_of_total,
            usable_df['total_counts'].sum(),
            df['total_counts'].sum(),
            args.threshold
        ]
    }
    
    summary_df = pd.DataFrame(summary_data)

    # Write summary statistics to output file
    try:
        summary_df.to_csv(args.output_summary, sep='\t', index=False)
    except Exception as e:
        sys.stderr.write(f"Error writing summary to {args.output_summary}: {e}\n")
        sys.exit(1)

    # Print summary to stdout
    print(f"Library: {args.libname}")
    print(f"Flagged {num_flagged_haps} haplotypes with counts per million > {args.threshold}")
    print(f"Flagged haplotypes account for {flagged_read_counts} reads ({flagged_fraction_of_usable:.4f} of usable reads)")
    if not args.use_all_haps and has_counted_codonwise:
        print(f"Flagged haplotypes account for {flagged_fraction_of_total:.4f} of total reads")
    print(f"Results written to {args.output_flagged} and {args.output_summary}")


def merge_command(args):
    """Merge flagged haplotypes from multiple samples"""
    # Load the sample table
    try:
        sample_table = pd.read_csv(args.sample_table, sep='\t')
    except Exception as e:
        sys.stderr.write(f"Error loading sample table {args.sample_table}: {e}\n")
        sys.exit(1)
    
    # Validate required columns
    required_cols = ['libname', 'filt_hap_tbl']
    missing_cols = [col for col in required_cols if col not in sample_table.columns]
    if missing_cols:
        sys.stderr.write(f"Error: Missing required columns in sample table: {missing_cols}\n")
        sys.stderr.write(f"Available columns: {list(sample_table.columns)}\n")
        sys.exit(1)
    
    # Check for duplicate libnames
    if sample_table['libname'].duplicated().any():
        duplicates = sample_table['libname'][sample_table['libname'].duplicated()].unique()
        sys.stderr.write(f"Error: Duplicate libnames found: {duplicates}\n")
        sys.exit(1)
    
    # Load all flagged haplotype files
    all_flagged_dfs = {}
    all_original_dfs = {}  # Store original haplotype tables for CPM calculation
    
    for _, row in sample_table.iterrows():
        libname = row['libname']
        filt_hap_tbl = row['filt_hap_tbl']
        
        try:
            # Load flagged haplotypes
            flagged_df = pd.read_csv(filt_hap_tbl, sep='\t')
            all_flagged_dfs[libname] = flagged_df
            print(f"Loaded {len(flagged_df)} flagged haplotypes from {libname}")
            
            # For CPM calculation, we need the original haplotype table
            # Try to infer the original table path from the filtered table path
            original_tbl_path = filt_hap_tbl.replace('_flagged.txt', '.byhap.txt')
            if not os.path.exists(original_tbl_path):
                # Try alternative naming patterns
                original_tbl_path = filt_hap_tbl.replace('_flagged', '')
                if not os.path.exists(original_tbl_path):
                    sys.stderr.write(f"Warning: Could not find original haplotype table for {libname}. CPM calculation will be skipped.\n")
                    all_original_dfs[libname] = None
                    continue
            
            # Load original haplotype table for CPM calculation
            original_df = pd.read_csv(original_tbl_path, sep='\t')
            all_original_dfs[libname] = original_df
            print(f"Loaded original haplotype table for {libname}")
            
        except Exception as e:
            sys.stderr.write(f"Error loading {filt_hap_tbl} for {libname}: {e}\n")
            sys.exit(1)
    
    if not all_flagged_dfs:
        sys.stderr.write("Error: No flagged haplotype files could be loaded\n")
        sys.exit(1)
    
    # Get all unique haplotypes across all samples
    all_haps = set()
    for df in all_flagged_dfs.values():
        all_haps.update(df['hap'].unique())
    
    # Create merged dataframe with all haplotypes
    merged_df = pd.DataFrame({'hap': list(all_haps)})
    
    # Add CPM and filtered status columns for each library
    for libname in sample_table['libname']:
        # Add CPM column
        cpm_col = f"{libname}_cpm"
        filtered_col = f"{libname}_filtered"
        
        # Initialize columns
        merged_df[cpm_col] = 0.0
        merged_df[filtered_col] = False
        
        # Get flagged haplotypes for this library
        if libname in all_flagged_dfs:
            flagged_haps = set(all_flagged_dfs[libname]['hap'])
            merged_df.loc[merged_df['hap'].isin(flagged_haps), filtered_col] = True
            
            # Calculate CPM if we have the original table
            if libname in all_original_dfs and all_original_dfs[libname] is not None:
                original_df = all_original_dfs[libname]
                
                # Calculate sum_counts (same logic as flag command)
                if 'counted_codonwise' in original_df.columns and not args.use_all_haps:
                    usable_df = original_df[original_df['counted_codonwise'] == True]
                    sum_counts = usable_df['total_counts'].sum()
                else:
                    sum_counts = original_df['total_counts'].sum()
                
                # Calculate CPM for each haplotype
                for _, hap_row in original_df.iterrows():
                    hap_id = hap_row['hap']
                    counts = hap_row['total_counts']
                    cpm = (counts / sum_counts) * 1e6 if sum_counts > 0 else 0.0
                    
                    # Update CPM for this haplotype
                    merged_df.loc[merged_df['hap'] == hap_id, cpm_col] = cpm
    
    # Remove sample-specific columns from the base haplotype information
    # We'll keep the first occurrence of each haplotype for the base info
    base_info_df = None
    for libname in sample_table['libname']:
        if libname in all_flagged_dfs:
            flagged_df = all_flagged_dfs[libname]
            if base_info_df is None:
                # Use the first library's data as the base
                base_info_df = flagged_df.copy()
                # Remove sample-specific columns
                columns_to_remove = ['total_counts', 'counts_per_million']
                existing_columns_to_remove = [col for col in columns_to_remove if col in base_info_df.columns]
                if existing_columns_to_remove:
                    base_info_df = base_info_df.drop(columns=existing_columns_to_remove)
                break
    
    if base_info_df is not None:
        base_info_df['hap'] = base_info_df['hap'].astype(object)
        merged_df['hap'] = merged_df['hap'].astype(object)

        # Merge the base info with our CPM/filtered data
        final_df = base_info_df.merge(merged_df, on='hap', how='right')
    else:
        final_df = merged_df
    
    # Write merged results
    try:
        final_df.to_csv(args.output_merged, sep='\t', index=False)
        print(f"Wrote {len(final_df)} haplotypes with CPM and filtered status to {args.output_merged}")
    except Exception as e:
        sys.stderr.write(f"Error writing merged results to {args.output_merged}: {e}\n")
        sys.exit(1)


def main():
    parser = argparse.ArgumentParser(
        description='Flag haplotypes with excessively high counts',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Flag command - Basic usage with default threshold of 1000 counts per million
  flag_highcount_haps flag --input_counts sample.byhap.txt --output_flagged flagged_haps.txt --output_summary summary.txt --libname sample1
  
  # Flag command - Use custom threshold of 500 counts per million
  flag_highcount_haps flag --input_counts sample.byhap.txt --threshold 500 --output_flagged flagged_haps.txt --output_summary summary.txt --libname sample1
  
  # Flag command - Include all haplotypes in sum_counts calculation
  flag_highcount_haps flag --input_counts sample.byhap.txt --use_all_haps --output_flagged flagged_haps.txt --output_summary summary.txt --libname sample1
  
  # Merge command - Combine flagged haplotypes from multiple samples
  flag_highcount_haps merge --sample_table samples.tsv --output_merged merged_flagged.txt
        """
    )
    
    subparsers = parser.add_subparsers(dest='command', help='Available commands')
    
    # Flag command parser
    flag_parser = subparsers.add_parser('flag', help='Flag haplotypes with excessively high counts from a single sample')
    flag_parser.add_argument('--input_counts', required=True, dest='input_counts',
                             help='Input per-sample haplotype counts table (output from codonwise_tally_vars_hapaware.py)')
    flag_parser.add_argument('--output_flagged', required=True, dest='output_flagged',
                             help='Output file for flagged haplotypes with all columns (tab-separated)')
    flag_parser.add_argument('--output_summary', required=True, dest='output_summary',
                             help='Output file for summary statistics table')
    flag_parser.add_argument('--libname', required=True, dest='libname',
                             help='Unique library name for this sample')
    flag_parser.add_argument('--threshold', default=1000, type=float, dest='threshold',
                             help='Threshold for counts per million (default: 1000)')
    flag_parser.add_argument('--use_all_haps', default=False, action='store_true', dest='use_all_haps',
                             help='Use all haplotypes for sum_counts calculation (default: only use counted_codonwise == True)')
    flag_parser.add_argument('--include_wildtype', action='store_true', default=False, dest='include_wildtype',
                             help='Include wildtype haplotypes in sum_counts calculation (default: exclude wildtype)')
    flag_parser.set_defaults(func=flag_command)

    # Merge command parser
    merge_parser = subparsers.add_parser('merge', help='Merge flagged haplotypes from multiple samples')
    merge_parser.add_argument('--sample_table', required=True, dest='sample_table',
                              help='TSV table with libname and filt_hap_tbl columns')
    merge_parser.add_argument('--output_merged', required=True, dest='output_merged',
                              help='Output file for merged flagged haplotypes with CPM and filtered status')
    merge_parser.add_argument('--use_all_haps', action='store_true', dest='use_all_haps',
                              help='Use all haplotypes for sum_counts calculation (default: only use counted_codonwise == True)')
    merge_parser.set_defaults(func=merge_command)

    args = parser.parse_args()

    if args.command is None:
        parser.print_help()
        sys.exit(1)

    # Call the appropriate command function
    args.func(args)


if __name__ == '__main__':
    main()