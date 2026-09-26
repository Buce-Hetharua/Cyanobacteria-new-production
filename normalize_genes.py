import pandas as pd
import numpy as np

# Load the gene run results and gene lengths
def load_data(gene_run_file, gene_length_file):
    # Load the gene run file: genes in the first column, run results in the other columns
    gene_run_df = pd.read_csv(gene_run_file, sep='\t', header=0, index_col=0)
    
    # Ensure the gene run results are numeric, replace non-numeric values with NaN
    gene_run_df = gene_run_df.apply(pd.to_numeric, errors='coerce')
    
    # Replace NaN values with 0 (for missing gene run results)
    gene_run_df.fillna(0, inplace=True)
    
    # Load the gene length file: gene names in the first column, gene lengths in the second
    gene_length_df = pd.read_csv(gene_length_file, sep='\t', header=None, names=["Gene", "GeneLength"])
    
    # Ensure gene lengths are numeric
    gene_length_df["GeneLength"] = pd.to_numeric(gene_length_df["GeneLength"], errors='coerce')
    
    # Create a dictionary for gene lengths
    gene_length_dict = dict(zip(gene_length_df["Gene"], gene_length_df["GeneLength"]))
    
    return gene_run_df, gene_length_dict

# Normalize the gene run results by gene length
def normalize_data(gene_run_df, gene_length_dict):
    # Convert the entire DataFrame to float64 to avoid dtype issues
    normalized_df = gene_run_df.copy().astype('float64')
    
    # Iterate over each gene in the gene run results
    for gene in gene_run_df.index:
        if gene in gene_length_dict:
            # Apply the correct normalization formula:
            # (Gene Read Count / Gene Length) * 1000
            normalized_df.loc[gene] = (gene_run_df.loc[gene] / gene_length_dict[gene]) * 1000
    
    # Replace infinity values with 0 (as requested)
    normalized_df.replace([np.inf, -np.inf], 0, inplace=True)
    
    return normalized_df

# Save the normalized results to a new file
def save_normalized_data(normalized_df, output_file):
    # Save the normalized results to a new TSV file
    normalized_df.to_csv(output_file, sep='\t', index=True)

# Main function to run the entire process
def main(gene_run_file, gene_length_file, output_file):
    # Load the data
    gene_run_df, gene_length_dict = load_data(gene_run_file, gene_length_file)
    
    # Normalize the data
    normalized_df = normalize_data(gene_run_df, gene_length_dict)
    
    # Save the result
    save_normalized_data(normalized_df, output_file)
    print(f"Normalized data saved to {output_file}")

if __name__ == "__main__":
    # File paths (replace these with the actual file paths)
    gene_run_file = "DATA-Raw-N95.tsv"
    gene_length_file = "N-seqs-length.tsv"
    output_file = "normalized_genelenght_N95.tsv"
    
    # Run the script
    main(gene_run_file, gene_length_file, output_file)
