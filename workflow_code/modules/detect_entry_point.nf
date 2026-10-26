process DETECT_ENTRY_POINT {
    input:
    path runsheet
    
    output:
    stdout emit: detected_entry_point
    
    script:
    """
    #!/usr/bin/env python3
    import pandas as pd
    import sys
    
    try:
        df = pd.read_csv("${runsheet}")
        
        # Check columns in priority order (earliest to latest step)
        if 'read1_path' in df.columns and df['read1_path'].notna().all():
            print("raw_reads", end="")
        elif 'trimmed_read1_path' in df.columns and df['trimmed_read1_path'].notna().all():
            print("trimmed_reads", end="")
        elif 'bam_path' in df.columns and df['bam_path'].notna().all():
            print("bam_files", end="")
        elif 'genes_results_path' in df.columns and df['genes_results_path'].notna().all():
            print("genes_results", end="")
        elif 'counts_table_path' in df.columns and df['counts_table_path'].notna().all():
            print("counts_table", end="")
        elif 'dge_table_path' in df.columns and df['dge_table_path'].notna().all():
            print("dge_table", end="")
        else:
            print("ERROR: No valid entry point columns found in runsheet", file=sys.stderr)
            sys.exit(1)
            
    except Exception as e:
        print(f"ERROR: Failed to parse runsheet: {e}", file=sys.stderr)
        sys.exit(1)
    """
}
