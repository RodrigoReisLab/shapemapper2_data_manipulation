# **Manipulation of probing data from shapemapper2**

## **Description**

Python pipeline that automates the extraction, processing, analysis, and visualization of probing data output from [shapemapper2](https://github.com/Weeks-UNC/shapemapper2). 

The pipeline extracts data from "_profile.txt" found in nested directory structure by searching for a simple naming convention (currently, based on Arabidopsis thaliana gene name, i.e., TAIR ID)

## **Features**

* **Flexible data discovery:** Automatically locates _profile.txt files in both immediate subdirectories (./data_dir/) and nested subdirectories (./parent_dir/data_dir/).  
* **Automates data aggregation:** For each gene, combines data from multiple treatment conditions into a single .tsv file.  
* **Basic probing data visualization:** Generates multiple plots for each gene:  
  * Comparative line plots for normalized profiles.  
  * Box plots showing the distribution of **probing rates**.  
  * A four-panel plot comparing **probing rates for each nucleotide (A, U, C, G)**.  
  * A heatmap of the Pearson correlation between normalized profiles.  
* **Statistical analysis:** Calculates summary statistics for each nucleotide and computes correlation matrices between different experimental conditions (also plotted as a heatmap).  
* **Organized output:** Creates a structured output directory (gene_data_output/) with a dedicated subfolder for each gene, containing all generated tables and plots.

## **Expected directory structure**

The script is designed to **work with either one of the following two directory structures**. You should run the script from the root directory that contains your experimental folders.

**Structure 1**: Child Directory  
The _profile.txt file is located in an immediate subdirectory.  

/your_project_root/  
|-- treatment_GeneID_1_.../  
|   |-- ...\_profile.txt  
|-- treatment_GeneID_2_.../  
|   |-- ...\_profile.txt  
|-- ProbingData_manipulation.py  

**Structure 2**: Grandchild Directory  
The \_profile.txt file is located in a nested subdirectory.  

/your_project_root/  
|-- directory_1/  
|   |-- treatment_GeneID_1_.../  
|   |   |-- ...\_profile.txt  
|-- directory_2/  
|   |-- treatment_GeneID_2_.../  
|   |   |-- ...\_profile.txt  
|-- ProbingData_manipulation.py  

**CRITICAL**: The **name** of the directory containing _profile.txt file must have a valid gene ID (currently, TAIR ID only) with a prefix (typically, a short name for the treatment). Suffix to gene ID is ignored.

### **Prerequisites**

* Python 3.8 or newer  

## ** Dependencies**

* [pandas](https://pandas.pydata.org/)
* [matplotlib](https://matplotlib.org/)
* [seaborn](https://seaborn.pydata.org/)

### **Usage**

1. Place the ProbingData_manipulation.py script in your root project directory, alongside your data folders.  
2. Go to the root directory in your terminal.  
3. Run the script:  
   python ProbingData_manipulation.py

The script will automatically find the data, process it, and create a gene_data_output directory with all the results.

## **Output files**

For each gene (e.g., AT3G23000), the script will generate a dedicated folder containing the following files:

* AT3G23000_combined_data.tsv: A table combining the Modified_rate, Untreated_rate, and Norm_profile data from all conditions.  
* AT3G23000_nucleotide_summary.tsv: A table with summary statistics for each nucleotide across all conditions.  
* AT3G23000_sliding_window.tsv: The results of the sliding window normalization.  
* AT3G23000_norm_profile_correlation.tsv: A table of Pearson correlation coefficients between the Norm_profile of each condition.  
* AT3G23000_combined_plot.png: A line plot comparing the Norm_profile across treatments.  
* AT3G23000_rates_boxplot.png: A box plot comparing the distribution of probing rates.  
* AT3G23000_nucleotide_rates_boxplot.png: A 4-panel figure with box plots for each nucleotide type.  
* AT3G23000_sliding_window_plot.png: A line plot of the sliding window data.  
* AT3G23000_norm_profile_correlation_heatmap.png: A heatmap visualizing the correlation matrix.
