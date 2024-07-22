import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr, spearmanr, kendalltau, ttest_rel, wilcoxon
from sklearn.metrics import mean_absolute_error, mean_squared_error
import numpy as np
import glob
from Bio import SeqIO

def find_quant_dir(result_dir, tool):
    result_dir = result_dir.strip()
    if not os.path.exists(result_dir):
        raise ValueError(f"Directory does not exist: {result_dir}")
    
    if tool == 'kallisto':
        pattern = os.path.join(result_dir, 'kallisto_quant_k*')
    elif tool == 'salmon':
        pattern = os.path.join(result_dir, 'salmon_quant*')
    else:
        raise ValueError(f"Unknown tool: {tool}")
    
    quant_dirs = glob.glob(pattern)
    if not quant_dirs:
        raise ValueError(f"No quantification directory found in {result_dir}")
    elif len(quant_dirs) > 1:
        print(f"Multiple quantification directories found in {result_dir}:")
        for i, dir in enumerate(quant_dirs):
            print(f"{i+1}. {os.path.basename(dir)}")
        choice = int(input("Enter the number of the directory you want to use: ")) - 1
        return quant_dirs[choice]
    else:
        return quant_dirs[0]

def find_fasta_file(result_dir):
    # Search for both .fasta and .fa files
    patterns = [
        os.path.join(result_dir, 'cleaned*.fasta'),
        os.path.join(result_dir, 'cleaned*.fa'),
        os.path.join(result_dir, '*.fasta'),
        os.path.join(result_dir, '*.fa')
    ]
    
    fasta_files = []
    for pattern in patterns:
        fasta_files.extend(glob.glob(pattern))
    
    if not fasta_files:
        # If not found, search in the parent directory
        parent_dir = os.path.dirname(result_dir)
        patterns = [
            os.path.join(parent_dir, 'cleaned*.fasta'),
            os.path.join(parent_dir, 'cleaned*.fa'),
            os.path.join(parent_dir, '*.fasta'),
            os.path.join(parent_dir, '*.fa')
        ]
        for pattern in patterns:
            fasta_files.extend(glob.glob(pattern))
    
    if not fasta_files:
        print(f"No .fasta or .fa file found in {result_dir} or its parent directory")
        fasta_file = input("Please enter the full path to the FASTA file: ").strip()
        if not os.path.exists(fasta_file):
            raise ValueError(f"Provided FASTA file does not exist: {fasta_file}")
        return fasta_file
    
    if len(fasta_files) > 1:
        print("Multiple FASTA files found:")
        for i, file in enumerate(fasta_files):
            print(f"{i+1}. {file}")
        choice = int(input("Enter the number of the file you want to use: ")) - 1
        return fasta_files[choice]
    
    return fasta_files[0]

def parse_fasta(fasta_file):
    total_transcripts = sum(1 for _ in SeqIO.parse(fasta_file, "fasta"))
    return total_transcripts
def parse_qc_report(qc_file):
    with open(qc_file, 'r') as f:
        lines = f.readlines()

    qc_data = []
    sample_data = {}
    for line in lines:
        if line.startswith("Sample:"):
            if sample_data:
                qc_data.append(sample_data)
            sample_data = {'Sample': line.split(":")[1].strip()}
        elif "Pseudoalignment rate" in line or "Reads mapped" in line:
            key, value = line.split(":")
            sample_data['MappingRate'] = float(value.strip().replace('%', ''))
    if sample_data:
        qc_data.append(sample_data)
    
    return pd.DataFrame(qc_data)

def parse_kallisto_results(kallisto_dir, qc_file):
    results = []
    qc_metrics = parse_qc_report(qc_file)
    
    for sample_dir in os.listdir(kallisto_dir):
        sample_path = os.path.join(kallisto_dir, sample_dir)
        if os.path.isdir(sample_path):
            abundance_file = os.path.join(sample_path, "abundance.tsv")
            
            if os.path.isfile(abundance_file):
                df = pd.read_csv(abundance_file, sep='\t')
                df['Sample'] = sample_dir
                results.append(df)
    
    if not results:
        raise ValueError(f"No abundance.tsv files found in {kallisto_dir}")
    
    combined_results = pd.concat(results, ignore_index=True)
    return combined_results, qc_metrics

def parse_salmon_results(salmon_dir, qc_file):
    results = []
    qc_metrics = parse_qc_report(qc_file)
    
    for sample_dir in os.listdir(salmon_dir):
        sample_path = os.path.join(salmon_dir, sample_dir)
        if os.path.isdir(sample_path):
            quant_file = os.path.join(sample_path, "quant.sf")
            
            if os.path.isfile(quant_file):
                df = pd.read_csv(quant_file, sep='\t')
                df['Sample'] = sample_dir
                results.append(df)
    
    if not results:
        raise ValueError(f"No quant.sf files found in {salmon_dir}")
    
    combined_results = pd.concat(results, ignore_index=True)
    return combined_results, qc_metrics
def normalize_results(df, tool):
    if tool == 'kallisto':
        df = df.rename(columns={'target_id': 'Transcript', 'tpm': 'TPM', 'est_counts': 'Count', 'eff_length': 'EffLength'})
    elif tool == 'salmon':
        df = df.rename(columns={'Name': 'Transcript', 'TPM': 'TPM', 'NumReads': 'Count', 'EffectiveLength': 'EffLength'})
    return df

def compare_results(df1, df2, tool1, tool2):
    df1_norm = normalize_results(df1, tool1)
    df2_norm = normalize_results(df2, tool2)
    
    merged_df = pd.merge(df1_norm[['Transcript', 'TPM', 'Count', 'EffLength', 'Sample']],
                         df2_norm[['Transcript', 'TPM', 'Count', 'EffLength', 'Sample']],
                         on=['Transcript', 'Sample'], suffixes=(f'_{tool1}', f'_{tool2}'))
    
    if merged_df.empty:
        raise ValueError("Merged dataframe is empty. Check if the transcript IDs and sample names match between the two tools.")
    
    return merged_df

def calculate_statistics(merged_df, tool1, tool2):
    stats = {}
    for sample in merged_df['Sample'].unique():
        sample_df = merged_df[merged_df['Sample'] == sample]
        if sample_df.shape[0] < 2:
            print(f"Skipping sample {sample} due to insufficient data points.")
            continue
        
        sample_stats = {}
        for metric in ['TPM', 'Count', 'EffLength']:
            col1 = f'{metric}_{tool1}'
            col2 = f'{metric}_{tool2}'
            
            pearson_corr, p_value = pearsonr(sample_df[col1], sample_df[col2])
            spearman_corr, _ = spearmanr(sample_df[col1], sample_df[col2])
            kendall_tau, _ = kendalltau(sample_df[col1], sample_df[col2])
            
            sample_stats[f'{metric}_Pearson_Correlation'] = pearson_corr
            sample_stats[f'{metric}_Pearson_P_Value'] = p_value
            sample_stats[f'{metric}_Spearman_Correlation'] = spearman_corr
            sample_stats[f'{metric}_Kendall_Tau'] = kendall_tau
            
            t_stat, t_p_value = ttest_rel(sample_df[col1], sample_df[col2])
            w_stat, w_p_value = wilcoxon(sample_df[col1], sample_df[col2])
            
            sample_stats[f'{metric}_Paired_T_Test'] = {'t_statistic': t_stat, 'p_value': t_p_value}
            sample_stats[f'{metric}_Wilcoxon_Test'] = {'w_statistic': w_stat, 'p_value': w_p_value}
            
            sample_stats[f'{metric}_MAE'] = mean_absolute_error(sample_df[col1], sample_df[col2])
            sample_stats[f'{metric}_MSE'] = mean_squared_error(sample_df[col1], sample_df[col2])
            sample_stats[f'{metric}_RMSE'] = np.sqrt(sample_stats[f'{metric}_MSE'])
        
        stats[sample] = sample_stats
    return stats

def calculate_mapping_rate_statistics(merged_qc, tool1, tool2):
    stats = {}
    if merged_qc.shape[0] < 2:
        print("Not enough samples to calculate mapping rate statistics. Skipping.")
        return stats
    
    col1 = f'MappingRate_{tool1}'
    col2 = f'MappingRate_{tool2}'
    
    pearson_corr, p_value = pearsonr(merged_qc[col1], merged_qc[col2])
    spearman_corr, _ = spearmanr(merged_qc[col1], merged_qc[col2])
    kendall_tau, _ = kendalltau(merged_qc[col1], merged_qc[col2])
    
    stats = {
        'MappingRate_Pearson_Correlation': pearson_corr,
        'MappingRate_Pearson_P_Value': p_value,
        'MappingRate_Spearman_Correlation': spearman_corr,
        'MappingRate_Kendall_Tau': kendall_tau
    }
    
    return stats

def count_transcripts(df, tool, total_transcripts):
    transcripts_tpm = df[df[f'TPM_{tool}'] > 0]['Transcript'].nunique()
    transcripts_cpm = df[df[f'Count_{tool}'] > 0]['Transcript'].nunique()
    
    return {
        f'Transcripts_TPM_{tool}': transcripts_tpm,
        f'Transcripts_CPM_{tool}': transcripts_cpm,
        f'Total_Transcripts_{tool}': total_transcripts,
        f'Transcripts_TPM_Percentage_{tool}': transcripts_tpm / total_transcripts * 100,
        f'Transcripts_CPM_Percentage_{tool}': transcripts_cpm / total_transcripts * 100
    }
def plot_scatter_and_density(sample_df, col1, col2, metric, tool1, tool2, sample, output_dir):
    # Scatter plot
    plt.figure(figsize=(10, 8))
    sns.scatterplot(x=f'{col1}_log', y=f'{col2}_log', data=sample_df, alpha=0.6)
    plt.xlabel(f'{metric} ({tool1})')
    plt.ylabel(f'{metric} ({tool2})')
    plt.title(f'Comparison of {metric} values between {tool1} and {tool2} for {sample}')
    plt.savefig(os.path.join(output_dir, f'{sample}_{metric}_scatter.png'))
    plt.close()
    
    # Density plot
    plt.figure(figsize=(10, 8))
    try:
        sns.kdeplot(data=sample_df, x=f'{col1}_log', y=f'{col2}_log', cmap="YlGnBu", fill=True)
    except ValueError as e:
        print(f"Warning: Could not create density plot for {sample}, {metric}. Falling back to hexbin plot.")
        plt.hexbin(sample_df[f'{col1}_log'], sample_df[f'{col2}_log'], gridsize=20, cmap='YlGnBu')
        plt.colorbar(label='count')
    plt.xlabel(f'{metric} ({tool1})')
    plt.ylabel(f'{metric} ({tool2})')
    plt.title(f'Density plot of {metric} values between {tool1} and {tool2} for {sample}')
    plt.savefig(os.path.join(output_dir, f'{sample}_{metric}_density.png'))
    plt.close()
def plot_bland_altman(sample_df, col1, col2, metric, tool1, tool2, sample, output_dir):
    plt.figure(figsize=(10, 8))
    mean = (sample_df[f'{col1}_log'] + sample_df[f'{col2}_log']) / 2
    diff = sample_df[f'{col1}_log'] - sample_df[f'{col2}_log']
    md = np.mean(diff)
    sd = np.std(diff, axis=0)
    plt.scatter(mean, diff, alpha=0.6)
    plt.axhline(md, color='gray', linestyle='--')
    plt.axhline(md + 1.96*sd, color='gray', linestyle='--')
    plt.axhline(md - 1.96*sd, color='gray', linestyle='--')
    plt.xlabel(f'Mean of {metric} ({tool1} and {tool2})')
    plt.ylabel(f'Difference in {metric} ({tool1} - {tool2})')
    plt.title(f'Bland-Altman Plot of {metric} for {sample}')
    plt.savefig(os.path.join(output_dir, f'{sample}_{metric}_bland_altman.png'))
    plt.close()
def plot_comparison(merged_df, tool1, tool2, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    for sample in merged_df['Sample'].unique():
        sample_df = merged_df[merged_df['Sample'] == sample]
        for metric in ['TPM', 'Count', 'EffLength']:
            col1 = f'{metric}_{tool1}'
            col2 = f'{metric}_{tool2}'
            
            # Remove zero or negative values
            sample_df_filtered = sample_df[(sample_df[col1] > 0) & (sample_df[col2] > 0)].copy()
            if sample_df_filtered.shape[0] < 2:
                print(f"No sufficient positive values for {sample}, {metric}. Skipping plot.")
                continue
            
            # Log transformation
            sample_df_filtered.loc[:, f'{col1}_log'] = np.log10(sample_df_filtered[col1])
            sample_df_filtered.loc[:, f'{col2}_log'] = np.log10(sample_df_filtered[col2])
            
            plot_scatter_and_density(sample_df_filtered, col1, col2, metric, tool1, tool2, sample, output_dir)
            plot_bland_altman(sample_df_filtered, col1, col2, metric, tool1, tool2, sample, output_dir)
def plot_mapping_rate_comparison(merged_qc, tool1, tool2, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    
    if merged_qc.shape[0] < 2:
        print(f"Not enough samples to plot mapping rate comparison. Skipping.")
        return
    
    plt.figure(figsize=(10, 8))
    sns.scatterplot(x=f'MappingRate_{tool1}', y=f'MappingRate_{tool2}', data=merged_qc, alpha=0.6)
    plt.xlabel(f'Mapping Rate ({tool1})')
    plt.ylabel(f'Mapping Rate ({tool2})')
    plt.title(f'Comparison of Mapping Rates between {tool1} and {tool2}')
    plt.savefig(os.path.join(output_dir, 'MappingRate_scatter.png'))
    plt.close()
def plot_mapping_rate_comparison(merged_qc, tool1, tool2, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    
    if merged_qc.shape[0] < 2:
        print(f"Not enough samples to plot mapping rate comparison. Skipping.")
        return
    
    plt.figure(figsize=(10, 8))
    sns.scatterplot(x=f'MappingRate_{tool1}', y=f'MappingRate_{tool2}', data=merged_qc, alpha=0.6)
    plt.xlabel(f'Mapping Rate ({tool1})')
    plt.ylabel(f'Mapping Rate ({tool2})')
    plt.title(f'Comparison of Mapping Rates between {tool1} and {tool2}')
    plt.savefig(os.path.join(output_dir, 'MappingRate_scatter.png'))
    plt.close()
def plot_correlation_heatmap(stats, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    for sample, sample_stats in stats.items():
        correlation_data = {
            'TPM': sample_stats['TPM_Pearson_Correlation'],
            'Count': sample_stats['Count_Pearson_Correlation'],
            'EffLength': sample_stats['EffLength_Pearson_Correlation']
        }
        df = pd.DataFrame([correlation_data])
        
        plt.figure(figsize=(10, 6))
        sns.heatmap(df, annot=True, cmap='coolwarm', vmin=-1, vmax=1, center=0)
        plt.title(f'Correlation Heatmap for {sample}')
        plt.savefig(os.path.join(output_dir, f'{sample}_correlation_heatmap.png'))
        plt.close()
def plot_transcript_detection(transcript_counts, tool1, tool2, output_dir):
    plt.figure(figsize=(12, 6))
    
    # Data for plotting
    tools = [tool1, tool2]
    total_transcripts = [transcript_counts[f'Total_Transcripts_{tool1}'], transcript_counts[f'Total_Transcripts_{tool2}']]
    detected_transcripts = [transcript_counts[f'Transcripts_TPM_{tool1}'], transcript_counts[f'Transcripts_TPM_{tool2}']]
    percentages = [transcript_counts[f'Transcripts_TPM_Percentage_{tool1}'], transcript_counts[f'Transcripts_TPM_Percentage_{tool2}']]
    
    x = np.arange(len(tools))
    width = 0.35
    
    fig, ax = plt.subplots(figsize=(10, 6))
    rects1 = ax.bar(x - width/2, total_transcripts, width, label='Total Transcripts', color='lightblue')
    rects2 = ax.bar(x + width/2, detected_transcripts, width, label='Detected Transcripts', color='darkblue')
    
    ax.set_ylabel('Number of Transcripts')
    ax.set_title('Transcript Detection Comparison')
    ax.set_xticks(x)
    ax.set_xticklabels(tools)
    ax.legend()
    
    # Add percentage labels
    for i, rect in enumerate(rects2):
        height = rect.get_height()
        ax.annotate(f'{percentages[i]:.2f}%',
                    xy=(rect.get_x() + rect.get_width() / 2, height),
                    xytext=(0, 3),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='bottom')
    
    fig.tight_layout()
    plt.savefig(os.path.join(output_dir, 'transcript_detection_comparison.png'))
    plt.close()
def plot_individual_mapping_rates(merged_qc, tool1, tool2, output_dir):
    plt.figure(figsize=(12, 6))
    
    samples = merged_qc['Sample']
    mapping_rates_tool1 = merged_qc[f'MappingRate_{tool1}']
    mapping_rates_tool2 = merged_qc[f'MappingRate_{tool2}']
    
    x = np.arange(len(samples))
    width = 0.35
    
    fig, ax = plt.subplots(figsize=(max(8, len(samples)), 6))
    rects1 = ax.bar(x - width/2, mapping_rates_tool1, width, label=tool1, color='lightblue')
    rects2 = ax.bar(x + width/2, mapping_rates_tool2, width, label=tool2, color='darkblue')
    
    ax.set_ylabel('Mapping Rate (%)')
    ax.set_title('Mapping Rate Comparison by Sample')
    ax.set_xticks(x)
    ax.set_xticklabels(samples, rotation=45, ha='right')
    ax.legend()
    
    # Add percentage labels
    def autolabel(rects):
        for rect in rects:
            height = rect.get_height()
            ax.annotate(f'{height:.2f}%',
                        xy=(rect.get_x() + rect.get_width() / 2, height),
                        xytext=(0, 3),  # 3 points vertical offset
                        textcoords="offset points",
                        ha='center', va='bottom')

    autolabel(rects1)
    autolabel(rects2)
    
    fig.tight_layout()
    plt.savefig(os.path.join(output_dir, 'individual_mapping_rates_comparison.png'))
    plt.close()

def main():
    while True:
        try:
            dir1 = input("Enter the path to the first tool's results directory: ").strip()
            dir2 = input("Enter the path to the second tool's results directory: ").strip()
            output_dir = input("Enter the path to save the comparison results: ").strip()

            if not os.path.exists(dir1) or not os.path.exists(dir2):
                raise ValueError("One or both of the input directories do not exist.")

            tool1 = 'kallisto' if 'kallisto' in dir1.lower() else 'salmon'
            tool2 = 'salmon' if 'salmon' in dir2.lower() else 'kallisto'

            print(f"Detected tools: {tool1} and {tool2}")

            quant_dir1 = find_quant_dir(dir1, tool1)
            quant_dir2 = find_quant_dir(dir2, tool2)

            print(f"Using quantification directories:\n{quant_dir1}\n{quant_dir2}")

            qc_file1 = os.path.join(dir1, "qc_report.txt")
            qc_file2 = os.path.join(dir2, "qc_report.txt")

            fasta_file1 = find_fasta_file(dir1)
            print(f"Using FASTA file for {tool1}: {fasta_file1}")
            total_transcripts_1 = parse_fasta(fasta_file1)

            use_same_fasta = input("Do you want to use the same FASTA file for both tools? (y/n): ").strip().lower()
            if use_same_fasta == 'y':
                fasta_file2 = fasta_file1
                total_transcripts_2 = total_transcripts_1
            else:
                fasta_file2 = find_fasta_file(dir2)
                print(f"Using FASTA file for {tool2}: {fasta_file2}")
                total_transcripts_2 = parse_fasta(fasta_file2)

            print(f"Total transcripts in {fasta_file1}: {total_transcripts_1}")
            if fasta_file1 != fasta_file2:
                print(f"Total transcripts in {fasta_file2}: {total_transcripts_2}")

            if tool1 == 'kallisto':
                df1, qc1 = parse_kallisto_results(quant_dir1, qc_file1)
            else:
                df1, qc1 = parse_salmon_results(quant_dir1, qc_file1)

            if tool2 == 'kallisto':
                df2, qc2 = parse_kallisto_results(quant_dir2, qc_file2)
            else:
                df2, qc2 = parse_salmon_results(quant_dir2, qc_file2)

            merged_df = compare_results(df1, df2, tool1, tool2)
            stats = calculate_statistics(merged_df, tool1, tool2)
            
            merged_qc = pd.merge(qc1, qc2, on='Sample', suffixes=(f'_{tool1}', f'_{tool2}'))
            mapping_rate_stats = calculate_mapping_rate_statistics(merged_qc, tool1, tool2)

            # Count transcripts for each tool
            transcript_counts = {}
            transcript_counts.update(count_transcripts(merged_df, tool1, total_transcripts_1))
            transcript_counts.update(count_transcripts(merged_df, tool2, total_transcripts_2))
            print(f"Transcript counts:\n{transcript_counts}")

            plot_comparison(merged_df, tool1, tool2, output_dir)
            plot_mapping_rate_comparison(merged_qc, tool1, tool2, output_dir)
            plot_correlation_heatmap(stats, output_dir)
            plot_transcript_detection(transcript_counts, tool1, tool2, output_dir)
            plot_individual_mapping_rates(merged_qc, tool1, tool2, output_dir)

            excel_file = os.path.join(output_dir, "comparison_results.xlsx")
            with pd.ExcelWriter(excel_file) as writer:
                merged_df.to_excel(writer, sheet_name='Comparison', index=False)
                for sample, sample_stats in stats.items():
                    pd.DataFrame.from_dict(sample_stats, orient='index').to_excel(writer, sheet_name=f'Statistics_{sample}')
                pd.DataFrame.from_dict(mapping_rate_stats, orient='index').to_excel(writer, sheet_name='MappingRate_Stats')
                pd.DataFrame.from_dict(transcript_counts, orient='index', columns=['Value']).to_excel(writer, sheet_name='Transcript_Counts')

            print(f"Comparison results and plots saved to {output_dir}")
            break  # Exit the loop if everything was successful
        except Exception as e:
            print(f"An error occurred: {str(e)}")
            print(f"Error type: {type(e).__name__}")
            print(f"Error occurred in {e.__traceback__.tb_frame.f_code.co_filename} at line {e.__traceback__.tb_lineno}")
            retry = input("Do you want to try again? (y/n): ")
            if retry.lower() != 'y':
                break

if __name__ == "__main__":
    main()