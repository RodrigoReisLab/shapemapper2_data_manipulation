import os
import re
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Patch
import glob
import logging
from collections import defaultdict

# --- Configuration ---
# Configure logging to show progress and potential issues
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Regex to find the gene name (e.g., AT3G23000)
GENE_REGEX = re.compile(r'(AT\dG\d{5})')
# File pattern to find the data file within each subdirectory
PROFILE_FILE_SUFFIX = '_profile.txt'
# Main output directory for all processed files
OUTPUT_DIR_NAME = 'gene_data_output'
# --- End Configuration ---


def find_profile_file(directory_path):
    """
    Finds the specific profile data file within a given directory.

    Args:
        directory_path (str): The path to the subdirectory.

    Returns:
        str or None: The full path to the profile file if found, otherwise None.
    """
    search_pattern = os.path.join(directory_path, f'*{PROFILE_FILE_SUFFIX}')
    files_found = glob.glob(search_pattern)
    if files_found:
        # Assuming only one profile file exists per directory
        return files_found[0]
    logging.warning(f"No '{PROFILE_FILE_SUFFIX}' file found in {directory_path}")
    return None


def extract_info_from_name(dir_name):
    """
    Extracts the gene name and treatment condition from a directory name.

    Args:
        dir_name (str): The name of the subdirectory.

    Returns:
        tuple or None: A tuple containing (gene_name, treatment_name) if successful,
                       otherwise None.
    """
    match = GENE_REGEX.search(dir_name)
    if match:
        gene_name = match.group(1)
        # Treatment is everything before the gene name
        treatment_name = dir_name[:match.start()].strip('_')
        if not treatment_name:
            logging.warning(f"Could not determine a valid treatment name for '{dir_name}'. Using 'unknown_treatment'.")
            treatment_name = 'unknown_treatment'
        return gene_name, treatment_name
    logging.warning(f"No gene name found in directory '{dir_name}'")
    return None, None


def plot_combined_data(df, gene_name, output_dir):
    """
    Generates and saves a plot from the combined dataframe, comparing Norm_profile.

    Args:
        df (pd.DataFrame): The combined dataframe to plot.
        gene_name (str): The gene name for labeling.
        output_dir (str): The directory to save the plot in.
    """
    norm_profile_cols = [col for col in df.columns if col.startswith('Norm_profile')]

    if not norm_profile_cols:
        logging.warning(f"No 'Norm_profile' columns found for {gene_name} to plot.")
        return

    plt.figure(figsize=(12, 7))
    
    for col_name in norm_profile_cols:
        # Extract treatment from column name for the legend
        treatment_name = col_name.replace('Norm_profile_', '')
        plt.plot(df.index, df[col_name], label=treatment_name)

    plt.title(f'Combined Normalized Profiles for {gene_name}')
    plt.xlabel('Nucleotide Position')
    plt.ylabel('Normalized Value (Norm_profile)')
    plt.legend(title='Treatments')
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.tight_layout()

    plot_filename = f"{gene_name}_combined_plot.png"
    plot_filepath = os.path.join(output_dir, plot_filename)

    try:
        plt.savefig(plot_filepath)
        logging.info(f"Successfully saved combined plot to {plot_filepath}")
    except Exception as e:
        logging.error(f"Failed to save plot {plot_filepath}. Reason: {e}")
    finally:
        plt.close()


def plot_box_plots(df, gene_name, output_dir):
    """
    Generates and saves a box plot for non-zero Modified_rate and Untreated_rate.

    Args:
        df (pd.DataFrame): The combined dataframe to plot.
        gene_name (str): The gene name for labeling.
        output_dir (str): The directory to save the plot in.
    """
    modified_cols = sorted([col for col in df.columns if col.startswith('Modified_rate_')])
    untreated_cols = sorted([col for col in df.columns if col.startswith('Untreated_rate_')])
    
    plot_cols = modified_cols + untreated_cols
    
    if not plot_cols:
        logging.warning(f"No rate columns found for {gene_name} to create a box plot.")
        return

    data_to_plot = [df[col].dropna() for col in plot_cols]
    
    # Check if there is any data to plot after dropping NaNs
    if not any(len(d) > 0 for d in data_to_plot):
        logging.warning(f"No non-zero rate data found for {gene_name} to create a box plot.")
        return

    plt.figure(figsize=(max(8, len(plot_cols) * 0.8), 8))

    # Create labels from column names
    labels = [col.replace('Modified_rate_', 'Mod_').replace('Untreated_rate_', 'Unt_') for col in plot_cols]

    bplot = plt.boxplot(data_to_plot, patch_artist=True, labels=labels)
    
    # Add colors
    colors = ['#A8DADC'] * len(modified_cols) + ['#F1FAEE'] * len(untreated_cols)
    for patch, color in zip(bplot['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_edgecolor('#1D3557')

    # Add median values on top of each box
    for i, line in enumerate(bplot['medians']):
        median_val = line.get_ydata()[0]
        box = bplot['boxes'][i]
        # Correctly get the top of the box using the path vertices
        top_of_box = box.get_path().vertices[:, 1].max()
        plt.text(i + 1, top_of_box, f'{median_val:.3f}',
                 horizontalalignment='center', 
                 verticalalignment='bottom',
                 fontsize=8,
                 color='#1D3557',
                 bbox=dict(facecolor='white', alpha=0.5, edgecolor='none', boxstyle='round,pad=0.2'))

    plt.title(f'Distribution of Rates for {gene_name} (values > 0)')
    plt.ylabel('probing rate')
    plt.ylim(0, 0.05)
    plt.xticks(rotation=45, ha="right")
    plt.grid(True, axis='y', linestyle='--', alpha=0.6)
    
    # Create a custom legend
    legend_elements = [Patch(facecolor='#A8DADC', edgecolor='#1D3557', label='Modified Rate'),
                       Patch(facecolor='#F1FAEE', edgecolor='#1D3557', label='Untreated Rate')]
    plt.legend(handles=legend_elements)

    plt.tight_layout()

    plot_filename = f"{gene_name}_rates_boxplot.png"
    plot_filepath = os.path.join(output_dir, plot_filename)

    try:
        plt.savefig(plot_filepath)
        logging.info(f"Successfully saved box plot to {plot_filepath}")
    except Exception as e:
        logging.error(f"Failed to save box plot {plot_filepath}. Reason: {e}")
    finally:
        plt.close()


def plot_nucleotide_boxplots(df, gene_name, output_dir):
    """
    Generates and saves a 4-panel box plot for probing rates of each nucleotide.

    Args:
        df (pd.DataFrame): The combined dataframe.
        gene_name (str): The gene name for labeling.
        output_dir (str): The directory to save the plot in.
    """
    nucleotides = ['A', 'U', 'C', 'G']
    fig, axes = plt.subplots(2, 2, figsize=(16, 14), sharey=True)
    axes = axes.flatten()

    modified_cols = sorted([col for col in df.columns if col.startswith('Modified_rate_')])
    untreated_cols = sorted([col for col in df.columns if col.startswith('Untreated_rate_')])
    
    if not (modified_cols or untreated_cols):
        logging.warning(f"No rate columns found for {gene_name} to create nucleotide box plots.")
        plt.close(fig)
        return

    # Extract unique treatments robustly from column names
    treatments_from_mod = [c.replace('Modified_rate_', '') for c in modified_cols]
    treatments_from_unt = [c.replace('Untreated_rate_', '') for c in untreated_cols]
    treatments = sorted(list(set(treatments_from_mod + treatments_from_unt)))


    for i, nucleotide in enumerate(nucleotides):
        ax = axes[i]
        data_for_nuc = df[df['Sequence'] == nucleotide]
        
        data_to_plot = []
        labels = []
        colors = []
        
        # First, gather all modified rates
        for treatment in treatments:
            mod_col = f'Modified_rate_{treatment}'
            if mod_col in data_for_nuc.columns:
                data_points = data_for_nuc[mod_col].dropna()
                if not data_points.empty:
                    data_to_plot.append(data_points)
                    labels.append(f'Mod_{treatment}')
                    colors.append('#A8DADC')
        
        # Second, gather all untreated rates
        for treatment in treatments:
            unt_col = f'Untreated_rate_{treatment}'
            if unt_col in data_for_nuc.columns:
                data_points = data_for_nuc[unt_col].dropna()
                if not data_points.empty:
                    data_to_plot.append(data_points)
                    labels.append(f'Unt_{treatment}')
                    colors.append('#F1FAEE')
        
        if data_to_plot:
            bplot = ax.boxplot(data_to_plot, patch_artist=True, labels=labels)

            for patch, color in zip(bplot['boxes'], colors):
                patch.set_facecolor(color)
                patch.set_edgecolor('#1D3557')
            
            # Add median values on top of each box
            for j, line in enumerate(bplot['medians']):
                median_val = line.get_ydata()[0]
                box = bplot['boxes'][j]
                top_of_box = box.get_path().vertices[:, 1].max()
                ax.text(j + 1, top_of_box, f'{median_val:.3f}',
                         horizontalalignment='center', 
                         verticalalignment='bottom',
                         fontsize=8,
                         color='#1D3557',
                         bbox=dict(facecolor='white', alpha=0.5, edgecolor='none', boxstyle='round,pad=0.2'))


            ax.set_title(f'Probing Rates for Nucleotide {nucleotide}')
            ax.set_ylabel('probing rate')
            ax.set_ylim(0, 0.05)
            ax.tick_params(axis='x', rotation=45, labelsize=8)
            ax.grid(True, axis='y', linestyle='--', alpha=0.6)
        else:
            ax.text(0.5, 0.5, f'No data for nucleotide {nucleotide}', horizontalalignment='center', verticalalignment='center')
            ax.set_title(f'Probing Rates for Nucleotide {nucleotide}')

    fig.suptitle(f'Nucleotide-Specific Probing Rates for {gene_name}', fontsize=16)
    
    # Common legend for the figure
    legend_elements = [Patch(facecolor='#A8DADC', edgecolor='#1D3557', label='Modified Rate'),
                       Patch(facecolor='#F1FAEE', edgecolor='#1D3557', label='Untreated Rate')]
    fig.legend(handles=legend_elements, loc='upper right')

    fig.tight_layout(rect=[0, 0.03, 1, 0.95]) # Adjust layout to make room for suptitle

    plot_filename = f"{gene_name}_nucleotide_rates_boxplot.png"
    plot_filepath = os.path.join(output_dir, plot_filename)
    
    try:
        plt.savefig(plot_filepath)
        logging.info(f"Successfully saved nucleotide box plot to {plot_filepath}")
    except Exception as e:
        logging.error(f"Failed to save nucleotide box plot {plot_filepath}. Reason: {e}")
    finally:
        plt.close(fig)


def calculate_and_save_sliding_window(df, gene_name, output_dir):
    """
    Calculates a sliding window normalization for Norm_profile columns and saves to a TSV.
    For each value, it is divided by the median of a 51-row forward-looking window.

    Args:
        df (pd.DataFrame): The combined dataframe.
        gene_name (str): The gene name for labeling.
        output_dir (str): The directory to save the file in.
        
    Returns:
        pd.DataFrame or None: The dataframe with sliding window calculations, or None if no data.
    """
    norm_profile_cols = [col for col in df.columns if col.startswith('Norm_profile_')]
    
    if not norm_profile_cols:
        logging.info(f"No 'Norm_profile' columns found for {gene_name} to calculate sliding window.")
        return None

    logging.info(f"Calculating sliding window normalization for {gene_name}...")
    
    # Start with base columns for context
    sliding_window_df = df[['Nucleotide', 'Sequence']].copy()
    window_size = 51

    for col_name in norm_profile_cols:
        series = df[col_name]
        
        # Calculate a forward-looking rolling median.
        # .rolling() is backward-looking by default (places result at the end of the window).
        # We shift the result by -(window_size - 1) to align the median of 
        # window [i:i+window_size] with index i.
        rolling_median = series.rolling(window=window_size, min_periods=1).median().shift(-(window_size - 1))
        
        # Calculate the normalized value
        new_col_name = f"SlidingWindow_{col_name.replace('Norm_profile_', '')}"
        sliding_window_df[new_col_name] = series / rolling_median

    output_filename = f"{gene_name}_sliding_window.tsv"
    output_filepath = os.path.join(output_dir, output_filename)

    try:
        sliding_window_df.to_csv(output_filepath, sep='\t', index=False, float_format='%.5f')
        logging.info(f"Successfully saved sliding window data to {output_filepath}")
    except Exception as e:
        logging.error(f"Failed to save sliding window data to {output_filepath}. Reason: {e}")
        
    return sliding_window_df


def plot_sliding_window_data(df, gene_name, output_dir):
    """
    Generates and saves a plot from the sliding window dataframe.

    Args:
        df (pd.DataFrame): The dataframe containing sliding window calculations.
        gene_name (str): The gene name for labeling.
        output_dir (str): The directory to save the plot in.
    """
    if df is None or df.empty:
        logging.warning(f"Sliding window dataframe is empty for {gene_name}. Skipping plot.")
        return
        
    sliding_window_cols = [col for col in df.columns if col.startswith('SlidingWindow_')]

    if not sliding_window_cols:
        logging.warning(f"No 'SlidingWindow' columns found for {gene_name} to plot.")
        return

    plt.figure(figsize=(12, 7))
    
    for col_name in sliding_window_cols:
        # Extract treatment from column name for the legend
        treatment_name = col_name.replace('SlidingWindow_', '')
        plt.plot(df.index, df[col_name], label=treatment_name, alpha=0.8)

    plt.title(f'Sliding Window Normalized Profiles for {gene_name}')
    plt.xlabel('Nucleotide Position')
    plt.ylabel('Normalized Value (Sliding Window)')
    plt.legend(title='Treatments')
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.tight_layout()

    plot_filename = f"{gene_name}_sliding_window_plot.png"
    plot_filepath = os.path.join(output_dir, plot_filename)

    try:
        plt.savefig(plot_filepath)
        logging.info(f"Successfully saved sliding window plot to {plot_filepath}")
    except Exception as e:
        logging.error(f"Failed to save plot {plot_filepath}. Reason: {e}")
    finally:
        plt.close()


def calculate_and_plot_correlation(df, gene_name, output_dir):
    """
    Calculates, saves, and plots the correlation matrix for Norm_profile columns.

    Args:
        df (pd.DataFrame): The combined dataframe.
        gene_name (str): The gene name for labeling.
        output_dir (str): The directory to save the outputs.
    """
    norm_profile_cols = [col for col in df.columns if col.startswith('Norm_profile_')]

    if len(norm_profile_cols) < 2:
        logging.info(f"Skipping correlation analysis for {gene_name}: less than two 'Norm_profile' columns found.")
        return

    logging.info(f"Calculating correlation matrix for {gene_name}...")
    
    correlation_df = df[norm_profile_cols].corr(method='pearson')

    # Clean up column and index names for readability
    clean_names = [col.replace('Norm_profile_', '') for col in correlation_df.columns]
    correlation_df.columns = clean_names
    correlation_df.index = clean_names

    # Save the correlation matrix to a TSV file
    output_table_name = f"{gene_name}_norm_profile_correlation.tsv"
    output_table_path = os.path.join(output_dir, output_table_name)
    try:
        correlation_df.to_csv(output_table_path, sep='\t', float_format='%.4f')
        logging.info(f"Successfully saved correlation table to {output_table_path}")
    except Exception as e:
        logging.error(f"Failed to save correlation table to {output_table_path}. Reason: {e}")

    # --- Create and save the heatmap ---
    plt.figure(figsize=(10, 8))
    
    sns.heatmap(correlation_df, annot=True, cmap='viridis', fmt='.2f', linewidths=.5)
    
    plt.title(f'Correlation of Normalized Profiles for {gene_name}')
    plt.xticks(rotation=45, ha="right")
    plt.yticks(rotation=0)
    plt.tight_layout()

    plot_filename = f"{gene_name}_norm_profile_correlation_heatmap.png"
    plot_filepath = os.path.join(output_dir, plot_filename)
    
    try:
        plt.savefig(plot_filepath)
        logging.info(f"Successfully saved correlation heatmap to {plot_filepath}")
    except Exception as e:
        logging.error(f"Failed to save correlation heatmap to {plot_filepath}. Reason: {e}")
    finally:
        plt.close()


def process_gene_group(gene_name, dir_info_tuples, base_output_dir):
    """
    Processes a group of directories corresponding to a single gene.
    Combines specified columns from their profile files into one table and plot.
    """
    logging.info(f"--- Processing Gene: {gene_name} ---")
    
    output_gene_dir = os.path.join(base_output_dir, gene_name)
    os.makedirs(output_gene_dir, exist_ok=True)
    
    combined_df = None
    base_cols = ["Nucleotide", "Sequence"]
    data_cols_to_extract = ["Modified_rate", "Untreated_rate", "Norm_profile"]

    # dir_info_tuples contains (directory_path, info_directory_name)
    for dir_path, info_dir_name in sorted(dir_info_tuples):
        _, treatment_name = extract_info_from_name(info_dir_name)
        
        profile_file = find_profile_file(dir_path)
        if not profile_file:
            continue
            
        try:
            temp_df = pd.read_csv(profile_file, sep='\t')
            
            # Filter out non-positive values for rate columns
            temp_df['Modified_rate'] = temp_df['Modified_rate'].where(temp_df['Modified_rate'] > 0)
            temp_df['Untreated_rate'] = temp_df['Untreated_rate'].where(temp_df['Untreated_rate'] > 0)

            required_cols = base_cols + data_cols_to_extract
            if not all(col in temp_df.columns for col in required_cols):
                logging.warning(f"Skipping {profile_file}: missing one or more required columns.")
                continue

            if combined_df is None:
                combined_df = temp_df[base_cols].copy()

            data_to_add = temp_df[data_cols_to_extract].copy()
            data_to_add.columns = [f"{col}_{treatment_name}" for col in data_to_add.columns]
            
            combined_df = pd.concat([combined_df, data_to_add], axis=1)

        except Exception as e:
            logging.error(f"Failed to process file {profile_file}. Reason: {e}")

    if combined_df is not None and not combined_df.empty:
        # --- Reorder columns before saving ---
        modified_rate_cols = sorted([col for col in combined_df.columns if col.startswith('Modified_rate_')])
        untreated_rate_cols = sorted([col for col in combined_df.columns if col.startswith('Untreated_rate_')])
        norm_profile_cols = sorted([col for col in combined_df.columns if col.startswith('Norm_profile_')])
        
        new_column_order = base_cols + modified_rate_cols + untreated_rate_cols + norm_profile_cols
        combined_df = combined_df[new_column_order]

        output_table_name = f"{gene_name}_combined_data.tsv"
        output_table_path = os.path.join(output_gene_dir, output_table_name)
        combined_df.to_csv(output_table_path, sep='\t', index=False)
        logging.info(f"Successfully created combined data table: {output_table_path}")
        
        # --- Create and save nucleotide summary table ---
        summary_data = []
        for nucleotide in ['A', 'U', 'C', 'G']:
            nuc_df = combined_df[combined_df['Sequence'] == nucleotide]
            for col in modified_rate_cols + untreated_rate_cols:
                stats_series = nuc_df[col].dropna().describe()
                if '50%' in stats_series:
                    stats_series = stats_series.rename({'50%': 'median'})
                
                stats = stats_series.to_dict()
                # Correctly split treatment from rate type
                parts = col.split('_', 1)
                rate_type = parts[0]
                treatment = parts[1].replace('rate_', '')

                row = {'Nucleotide': nucleotide, 'Rate_Type': rate_type, 'Treatment': treatment}
                row.update(stats)
                summary_data.append(row)
        
        summary_df = pd.DataFrame(summary_data)
        if not summary_df.empty:
            summary_table_name = f"{gene_name}_nucleotide_summary.tsv"
            summary_table_path = os.path.join(output_gene_dir, summary_table_name)
            summary_df.to_csv(summary_table_path, sep='\t', index=False, float_format='%.5f')
            logging.info(f"Successfully created nucleotide summary table: {summary_table_path}")

        # --- Calculate and save sliding window normalization ---
        sliding_window_df = calculate_and_save_sliding_window(combined_df, gene_name, output_gene_dir)

        # --- Plotting and Analysis Steps ---
        calculate_and_plot_correlation(combined_df, gene_name, output_gene_dir)
        plot_combined_data(combined_df, gene_name, output_gene_dir)
        plot_box_plots(combined_df, gene_name, output_gene_dir)
        plot_nucleotide_boxplots(combined_df, gene_name, output_gene_dir)
        plot_sliding_window_data(sliding_window_df, gene_name, output_gene_dir)
    else:
        logging.warning(f"No data was combined for gene {gene_name}. No output file will be created.")


def main():
    """
    Main function to find, group, and process data for each gene.
    It searches for profile files in immediate child and grandchild directories.
    """
    start_directory = os.getcwd()
    logging.info(f"Starting script in directory: {start_directory}")
    
    main_output_dir = os.path.join(start_directory, OUTPUT_DIR_NAME)
    os.makedirs(main_output_dir, exist_ok=True)
    
    gene_groups = defaultdict(list)
    
    # This dictionary maps the info_dir_name to the data_dir_path
    # to prevent processing the same dataset twice.
    data_dirs_to_process = {}

    # Search for profile files in child and grandchild directories
    for item in os.listdir(start_directory):
        child_path = os.path.join(start_directory, item)

        # Skip the output directory and any files
        if not os.path.isdir(child_path) or item == OUTPUT_DIR_NAME:
            continue

        # Case 1: Look for profile file directly in the child directory
        if glob.glob(os.path.join(child_path, f'*{PROFILE_FILE_SUFFIX}')):
            info_dir_name = item
            # Only add if we haven't found this treatment directory before
            if info_dir_name not in data_dirs_to_process:
                data_dirs_to_process[info_dir_name] = child_path
        
        # Case 2: Look for profile file in grandchild directories
        for sub_item in os.listdir(child_path):
            grandchild_path = os.path.join(child_path, sub_item)
            if os.path.isdir(grandchild_path):
                if glob.glob(os.path.join(grandchild_path, f'*{PROFILE_FILE_SUFFIX}')):
                    info_dir_name = sub_item
                    # Only add if we haven't found this treatment directory before
                    if info_dir_name not in data_dirs_to_process:
                        data_dirs_to_process[info_dir_name] = grandchild_path

    # Group the found directories by gene name
    for info_dir_name, data_dir in data_dirs_to_process.items():
        gene_name, _ = extract_info_from_name(info_dir_name)
        if gene_name:
            # Store a tuple of (data_directory_path, info_directory_name)
            gene_groups[gene_name].append((data_dir, info_dir_name))

    if not gene_groups:
        logging.warning("No subdirectories containing profile files matching the naming convention were found.")
        return

    for gene_name, dir_info_tuples in gene_groups.items():
        process_gene_group(gene_name, dir_info_tuples, main_output_dir)

    logging.info("Script finished processing all directories.")


if __name__ == "__main__":
    main()

