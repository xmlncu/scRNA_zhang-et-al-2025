# -*- coding: utf-8 -*-
"""
Created on Fri Dec 20 18:41:52 2024

@author: Xiaoming
"""

import os
print("Working Directory:", os.getcwd())

# Set your working directory
os.chdir(r"D:\Your\Working\Directory")

# Confirm the working directory
print("Current working directory:", os.getcwd())

# convert all .10x format to .h5ad 
import os
import scanpy as sc

# Define the folder paths
input_folder = r"D:\AI work\IBD8sample\10x format"
output_folder = r"D:\AI work\IBD8sample\h5ad"

# Create the output folder if it doesn't exist
os.makedirs(output_folder, exist_ok=True)

# Iterate over all files and directories in the input folder
for dir_name in os.listdir(input_folder):
    dir_path = os.path.join(input_folder, dir_name)
    
    # Check if it is a directory (10X matrix folders are typically directories)
    if os.path.isdir(dir_path):
        try:
            # Read the 10X matrix
            print(f"Reading 10X matrix from {dir_path}...")
            adata = sc.read_10x_mtx(dir_path, var_names='gene_symbols', cache=True)
            
            # Ensure unique variable names
            adata.var_names_make_unique()
            
            # Define the output file name
            output_file = os.path.join(output_folder, f"{dir_name}.h5ad")
            
            # Save as .h5ad
            print(f"Saving to {output_file}...")
            adata.write(output_file)
        except Exception as e:
            print(f"Error processing {dir_name}: {e}")
print("Conversion completed!")


# 1. necessary setting
# potential code for advanced analysis
# new anlysis and cell_label
# cell trajectory analysis
import omicverse as ov
#print(f"omicverse version: {ov.__version__}")
import scanpy as sc
#print(f"scanpy version: {sc.__version__}")
ov.utils.ov_plot_set()
sc.settings.set_figure_params(dpi=200, facecolor='white')
import os
import tempfile
import scvi
import seaborn as sns
import torch
scvi.settings.seed = 0
torch.set_float32_matmul_precision("high")
import scvelo as scv

# Core scverse libraries
import scanpy as sc
import anndata as ad

# Data retrieval
import pooch


# 2. Parameter settings

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

# Define a custom colormap
custom_cmap = LinearSegmentedColormap.from_list(
    'custom_blue_green', ['cyan', 'red'], N=256
)



# 2. Potential data filtration preprocessing part
import scanpy as sc
import os
from pathlib import Path

def process_h5ad(input_file, output_dir):
    """
    Process a single h5ad file with QC metrics and filtering
    
    Parameters:
    -----------
    input_file : str or Path
        Path to input h5ad file
    output_dir : str or Path
        Directory to save processed file
    """
    print(f"Processing {input_file}")
    
    # Read the data
    adata = sc.read_h5ad(input_file)
    
    # Calculate QC metrics
    # Mitochondrial genes
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    # Ribosomal genes
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
    # Hemoglobin genes
    adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")
    
    # Calculate QC metrics
    sc.pp.calculate_qc_metrics(
        adata, 
        qc_vars=["mt", "ribo", "hb"], 
        inplace=True, 
        log1p=True
    )
    
    # Plot QC metrics before filtering
    fig_dir = Path(output_dir) / "figures"
    fig_dir.mkdir(exist_ok=True)
    
    sc.pl.violin(
        adata,
        ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
        jitter=0.4,
        multi_panel=True,
        save=f"{Path(input_file).stem}_qc_before_filtering.pdf"
    )
    
    # Filter cells and genes
    sc.pp.filter_cells(adata, min_genes=100)
    sc.pp.filter_genes(adata, min_cells=3)
    
    # Save count data
    adata.layers["counts"] = adata.X.copy()
    
    # Filter based on mitochondrial percentage and number of genes
    adata = adata[adata.obs['pct_counts_mt'] < 20].copy()
    adata = adata[adata.obs['n_genes_by_counts'] < 10000].copy()
    adata = adata[adata.obs['n_genes_by_counts'] > 100].copy()
    
    # Save processed file
    output_file = Path(output_dir) / f"{Path(input_file).stem}_filtered.h5ad"
    adata.write_h5ad(output_file)
    
    print(f"Saved processed file to {output_file}")
    return

def batch_process_folder(input_dir, output_dir):
    """
    Process all h5ad files in a directory
    
    Parameters:
    -----------
    input_dir : str
        Input directory containing h5ad files
    output_dir : str
        Output directory for processed files
    """
    # Create output directory if it doesn't exist
    output_path = Path(output_dir)
    output_path.mkdir(exist_ok=True)
    
    # Process all h5ad files
    input_path = Path(input_dir)
    h5ad_files = list(input_path.glob("*.h5ad"))
    
    if not h5ad_files:
        print(f"No h5ad files found in {input_dir}")
        return
    
    print(f"Found {len(h5ad_files)} h5ad files to process")
    
    for file in h5ad_files:
        try:
            process_h5ad(file, output_path)
        except Exception as e:
            print(f"Error processing {file}: {str(e)}")
            continue
    
    print("Batch processing complete!")

# Run the batch processing
input_directory = r'I:\zhangtongmei\H5ad_Loom'
output_directory = r'I:\zhangtongmei\H5ad_Loom\filtered'
batch_process_folder(input_directory, output_directory)



# 3. SCVI-TOOLS for data integrations
import os
import scanpy as sc
import scvi

# Define the folder where your .h5ad files are located
data_folder = r"D:/AI work/HumanAD"

# Get a list of all .h5ad files in the directory
samples = {os.path.splitext(f)[0]: f for f in os.listdir(data_folder) if f.endswith(".h5ad")}

# Load all samples
adatas = {}
for sample_id, filename in samples.items():
    file_path = os.path.join(data_folder, filename)
    sample_adata = sc.read(file_path)
    sample_adata.var_names_make_unique()
    adatas[sample_id] = sample_adata

# Concatenate AnnData objects and ensure unique obs names
adata = sc.concat(adatas, label="sample")
adata.obs_names_make_unique()

# Print cell count per sample
print(adata.obs["sample"].value_counts())

# Prepare the AnnData object for scVI
scvi.model.SCVI.setup_anndata(adata, batch_key="sample")

# Initialize and train the SCVI model
model = scvi.model.SCVI(adata)
model.train()

# Get the latent representation
adata.obsm["X_scVI"] = model.get_latent_representation()

# Optional: Save the resulting AnnData object
adata.write_h5ad(r"D:/AI work/HumanAD/HumanADTau_scVI_adata.h5ad")

adata.layers["counts"] = adata.X.copy()  # preserve counts
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata  # freeze the state in `.raw`

# afterwards adata trainings
sc.pp.highly_variable_genes(
    adata_micro,
    n_top_genes=1200,
    subset=True,
    layer="counts",
    flavor="seurat_v3",
    batch_key="Group",
)
scvi.model.SCVI.setup_anndata(
    adata_micro,
    layer="counts",
    categorical_covariate_keys=["cell_type", "Group"],
    #continuous_covariate_keys=["percent_mito", "percent_ribo"],
)
model = scvi.model.SCVI(adata_micro)
model.train()
SCVI_LATENT_KEY = "X_scVI"

latent = model.get_latent_representation()
adata_micro.obsm[SCVI_LATENT_KEY] = latent
latent.shape
denoised = model.get_normalized_expression(adata_micro, library_size=1e4)
denoised.iloc[:5, :5]
SCVI_NORMALIZED_KEY = "scvi_normalized"

adata_micro.layers[SCVI_NORMALIZED_KEY] = model.get_normalized_expression(library_size=10e4)

# 3.2 optional for define name of group or cell_type/cell_label

adata.obs["Group"] = adata.obs["sample"].map(
    {
    'AD_Hig_1': 'High_Tau', 
    'AD_Hig_2': 'High_Tau',  
    'AD_Hig_3': 'High_Tau',
    'AD_Low_1': 'Non_AD', 
    'AD_Low_2': 'Non_AD',  
    'AD_Low_3': 'Non_AD',
    'AD_Low_4': 'Non_AD', 
    'AD_Mid_1': 'Low_Tau',
    'AD_Mid_2': 'Low_Tau', 
    'AD_Mid_3': 'Low_Tau',
    'AD_Mid_4': 'Low_Tau',     
    'AD_Mid_5': 'Low_Tau',
    'AD_Mid_6': 'Low_Tau',     
    'AD_Mid_7': 'Low_Tau',
    'AD_Mid_8': 'Low_Tau',     
    'AD_Mid_9': 'Low_Tau',
    'AD_Mid_10': 'Low_Tau',     
    'AD_Mid_11': 'Low_Tau',
    'AD_Mid_12': 'Low_Tau',          
    }
)




# 4. Post_SCVI-processing

import numpy as np
adata.X=adata.X.astype(np.int64)

adata.X.max()

adata.var_names_make_unique()
adata.obs_names_make_unique()


# data filitration %%time
adata=ov.pp.qc(adata,
              tresh={'mito_perc': 0.25, 'nUMIs': 500, 'detected_genes': 250},
               #doublets_method='sccomposite',
              batch_key='sample')
adata

# define hvg %%time
adata=ov.pp.preprocess(adata,mode='shiftlog|pearson',n_HVGs=3000,
                      batch_key=None)
adata

# store the raw data %%time
adata.raw = adata
adata = adata[:, adata.var.highly_variable_features]
adata

# scale the data for further processing %%time
ov.pp.scale(adata,max_value=12)
adata

# %%time define the pca dimension
ov.pp.pca(adata,layer='scaled',n_pcs=50)
adata
adata.obsm['X_pca']=adata.obsm['scaled|original|X_pca']
ov.pl.embedding(adata,
                  basis='X_pca',
                  color='GJA1',
                  frameon='small')

# Print the counts for each cell type
print("Cell numbers by 'cell_type':")
print(adata.obs['cell_label'].value_counts())

# Print the counts for each group
print("\nCell numbers by 'Group':")
print(adata.obs['Group'].value_counts())

# plot information 
marker_genes =['CDH5']
sc.pl.dotplot(adata, marker_genes, groupby='Sample',
             standard_scale='var');

# scvi UMAP PROCESS
sc.pp.highly_variable_genes(
    adata,
    n_top_genes=1200,
    subset=True,
    layer="counts",
    flavor="seurat_v3",
    batch_key="Group",
)
scvi.model.SCVI.setup_anndata(
    adata,
    layer="counts",
    categorical_covariate_keys=["Group","sample"],
    #continuous_covariate_keys=["percent_mito", "percent_ribo"],
)
model = scvi.model.SCVI(adata)
model.train()
SCVI_LATENT_KEY = "X_scVI"

latent = model.get_latent_representation()
adata.obsm[SCVI_LATENT_KEY] = latent
latent.shape
SCVI_NORMALIZED_KEY = "scvi_normalized"

adata.layers[SCVI_NORMALIZED_KEY] = model.get_normalized_expression(library_size=10e4)
# run PCA then generate UMAP plots
sc.tl.pca(adata)
sc.pp.neighbors(adata, n_pcs=30, n_neighbors=20)
sc.tl.umap(adata, min_dist=0.3)
sc.pl.umap(
    adata,
    color=["Group"],
    frameon=False,
)
sc.pl.umap(
    adata,
    color=["Group", "sample"],
    ncols=2,
    frameon=False,
)



# splite the group by leiden

# Convert columns to strings before concatenation
groupby_combined = adata.obs['Group'].astype(str) + "_" + adata.obs['cell_label'].astype(str)
adata.obs['group_split'] = groupby_combined

# Generate the dotplot
sc.pl.dotplot(
    adata, 
    marker_genes, 
    groupby='group_split', 
    standard_scale='var', 
    dendrogram=True
)

# Cleanup
del adata.obs['group_split']



# vln plot to show
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt

# Convert columns to strings before concatenation
groupby_combined = adata.obs['Group'].astype(str) + "_" + adata.obs['cell_label'].astype(str)
adata.obs['group_split'] = groupby_combined

# Define a custom color palette
custom_palette = sns.color_palette("Set2", len(adata.obs['group_split'].unique()))

# Set a style (you can adjust this to your liking)
sns.set_theme(style="whitegrid")

# Generate the violin plot with a specified figure size
plt.figure(figsize=(30, 10))
sc.pl.violin(
    adata,
    keys=['SORL1'],  # List of marker genes
    groupby='group_split',
    stripplot=False,  # Disable strip plot
    rotation=90,      # Rotate x-axis labels
    #figsize=(10, 6),  # Set figure size (width, height)
    show=False        # Don't show immediately to allow further customization
)

# Customize plot appearance with title and axis labels
plt.suptitle("Expression of SORL1 by cell_label", fontsize=16, y=1.02)
plt.xlabel("cell_label", fontsize=12)
plt.ylabel("SORL1 Expression", fontsize=12)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.tight_layout()

# Disable the grid
plt.grid(False)

# Remove the top and right spines
ax = plt.gca()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_color('black')
ax.spines['left'].set_color('black')

# Show the plot
plt.show()


# better bar plot visualization
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt

# Ensure proper concatenation of group information if needed
groupby_combined = adata.obs['Group'].astype(str) + "_" + adata.obs['cell_label'].astype(str)
adata.obs['group_split'] = groupby_combined

# Define a custom color palette for distinct cell types
custom_palette = sns.color_palette("Set2", len(adata.obs['cell_label'].unique()))

# Set a style (you can adjust this to your liking)
sns.set_theme(style="whitegrid")

# Generate the violin plot grouped by 'cell_label' with a specified figure size
plt.figure(figsize=(30, 10))  # Adjust figure size for better visualization
sc.pl.violin(
    adata,
    keys=['SORL1'],        # Gene to plot
    groupby='cell_label',  # Group by cell_label
    stripplot=False,       # Disable strip plot
    rotation=90,           # Rotate x-axis labels for readability
    show=False             # Don't show immediately to allow further customization
)

# Customize plot appearance with title and axis labels
plt.suptitle("Expression of SORL1 by Cell Type", fontsize=20, y=1.02)
plt.xlabel("cell_label", fontsize=16)
plt.ylabel("SORL1 Expression", fontsize=16)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.tight_layout()

# Disable the grid for a cleaner appearance
plt.grid(False)

# Customize axes by removing the top and right spines
ax = plt.gca()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_color('black')
ax.spines['left'].set_color('black')

# Show the plot
plt.show()


# cell_type visualization
import matplotlib.pyplot as plt
fig,ax=plt.subplots(figsize = (12,4))
ov.pl.cellproportion(adata=adata,celltype_clusters='cell_label',
                    groupby='sample',legend=True,ax=ax)


# leiden and batch group informations
sc.tl.dendrogram(adata,'cell_label',use_rep='scaled|original|X_pca')
sc.tl.rank_genes_groups(adata, 'cell_label', use_rep='scaled|original|X_pca',
                        method='t-test',use_raw=False,key_added='leiden_ttest')
sc.pl.rank_genes_groups_dotplot(adata,groupby='cell_label',
                                cmap='Spectral_r',key='leiden_ttest',
                                standard_scale='var',n_genes=5)

# pathway analysis 
ov.utils.download_pathway_database()
ov.utils.download_geneid_annotation_pair()

......Pathway Geneset download start: GO_Biological_Process_2021
......Downloading dataset save to genesets/GO_Biological_Process_2021.txt
......Creating directory genesets
......[GO_Biological_Process_2021 Size of file]: 0.15 MB
......[Downloader]: >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>100.00%
.......Finish！4.126541614532471.2f s
......Pathway Geneset download start: GO_Cellular_Component_2021
......Downloading dataset save to genesets/GO_Cellular_Component_2021.txt
......[GO_Cellular_Component_2021 Size of file]: 0.03 MB
......[Downloader]: >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>100.00%
.......Finish！2.420053482055664.2f s
......Pathway Geneset download start: GO_Molecular_Function_2021
......Downloading dataset save to genesets/GO_Molecular_Function_2021.txt
......[GO_Molecular_Function_2021 Size of file]: 0.03 MB
......[Downloader]: >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>100.00%
.......Finish！5.607283115386963.2f s
......Pathway Geneset download start: WikiPathway_2021_Human
......Downloading dataset save to genesets/WikiPathway_2021_Human.txt
......[WikiPathway_2021_Human Size of file]: 0.02 MB
......[Downloader]: >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>100.00%
.......Finish！2.965646743774414.2f s
......Pathway Geneset download start: WikiPathways_2019_Mouse
......Downloading dataset save to genesets/WikiPathways_2019_Mouse.txt
......[WikiPathways_2019_Mouse Size of file]: 0.01 MB
......[Downloader]: >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>100.00%
...
......[pair_danRer7 Size of file]: 0.09 MB
......[Downloader]: >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>100.00%
.......Finish！3.2669999599456787.2f s
......Geneid Annotation Pair download finished!


pathway_dict=ov.utils.geneset_prepare('genesets/GO_Biological_Process_2021.txt',organism='Human')

##Assest all pathways
adata_aucs=ov.single.pathway_aucell_enrichment(adata,
                                                  pathways_dict=pathway_dict,
                                                  num_workers=8)


adata_aucs.obs=adata[adata_aucs.obs.index].obs
adata_aucs.obsm=adata[adata_aucs.obs.index].obsm
adata_aucs.obsp=adata[adata_aucs.obs.index].obsp
adata_aucs

#adata_aucs.uns['log1p']['base']=None
sc.tl.rank_genes_groups(adata_aucs, 'cell_label', method='t-test',n_genes=100)
sc.pl.rank_genes_groups_dotplot(adata_aucs,groupby='cell_label',
                                cmap='Spectral_r',
                                standard_scale='var',n_genes=3)

adata.uns['log1p']['base']=None
sc.tl.rank_genes_groups(adata, 'cell_label', method='t-test',n_genes=100)

res=ov.single.pathway_enrichment(adata,pathways_dict=pathway_dict,organism='Human',
                                     group_by='cell_label',plot=True)

ax=ov.single.pathway_enrichment_plot(res,plot_title='Enrichment',cmap='Reds',
                                         xticklabels=True,cbar=False,square=True,vmax=10,
                                         yticklabels=True,cbar_kws={'label': '-log10(qvalue)','shrink': 0.5,})


# better visualization

# Use it in the plot
sc.pl.umap(
    adata,
    color=['leiden','CD14','CD3E','NKG7','COL1A1','EPHX2','IRAK3','Group','NFKB1','TNF','EPCAM'],
    cmap=custom_cmap,
    size=15,
    alpha=0.8,
    frameon=False,
    legend_loc = 'on data',
    show=True
)

# change it to the theme, thanks!
import matplotlib.pyplot as plt
fig,ax=plt.subplots(figsize = (1,4))
ov.pl.cellproportion(adata=adata,celltype_clusters='clusters',
                    groupby='age(days)',legend=True,ax=ax)



#subset the cell_type


# Ensure cell labels are strings
adata.obs['Group'] = adata.obs['Group'].astype(str)

# Remove NaN or missing values
adata = adata[~adata.obs['cell_label'].isnull()].copy()

# Run rank_genes_groups
sc.tl.rank_genes_groups(adata, 'cell_label', method='t-test', n_genes=100)

# Check valid cell types
valid_celltypes = adata.obs['cell_label'].unique()
print("Valid cell types:", valid_celltypes)

# Perform pathway enrichment
res = ov.single.pathway_enrichment(
    adata,
    pathways_dict=pathway_dict,
    organism='Human',
    group_by='cell_label',
    plot=True
)

# Plot pathway enrichment results
ax = ov.single.pathway_enrichment_plot(
    res,
    plot_title='Enrichment',
    cmap='Reds',
    xticklabels=True,
    cbar=False,
    square=True,
    vmax=10,
    yticklabels=True,
    cbar_kws={'label': '-log10(qvalue)', 'shrink': 0.5},
)


# convert human to mouse ID

prior_network = ov.single.convert_human_to_mouse_network(prior_network, server_name='ensembl.ensembl.org')  # Europe server




# Prepare pathway dictionaries
pathway_dict_BP = ov.utils.geneset_prepare(r'D:\AI work\hubeihospital\genesets\GO_Biological_Process_2021.txt', organism='Human')
pathway_dict_CC = ov.utils.geneset_prepare(r'D:\AI work\hubeihospital\genesets\GO_Cellular_Component_2021.txt', organism='Human')
pathway_dict_MF = ov.utils.geneset_prepare(r'D:\AI work\hubeihospital\genesets\GO_Molecular_Function_2021.txt', organism='Human')
pathway_dict_Rea = ov.utils.geneset_prepare(r'D:\AI work\hubeihospital\genesets\Reactome_2022.txt', organism='Human')
pathway_dict_wiki = ov.utils.geneset_prepare(r'D:\AI work\hubeihospital\genesets\WikiPathway_2021_Human.txt', organism='Human')

# Set log1p base to None
adata_micro.uns['log1p']['base'] = None

# Perform differential expression analysis
sc.tl.rank_genes_groups(adata_micro, 'leiden_scVI', method='t-test', n_genes=100)

# Perform pathway enrichment analysis for each gene set
res_BP = ov.single.pathway_enrichment(adata_micro, pathways_dict=pathway_dict_BP, organism='Human', group_by='leiden_scVI', plot=False)
res_CC = ov.single.pathway_enrichment(adata_micro, pathways_dict=pathway_dict_CC, organism='Human', group_by='leiden_scVI', plot=False)
res_MF = ov.single.pathway_enrichment(adata_micro, pathways_dict=pathway_dict_MF, organism='Human', group_by='leiden_scVI', plot=False)
res_Rea = ov.single.pathway_enrichment(adata_micro, pathways_dict=pathway_dict_Rea, organism='Human', group_by='leiden_scVI', plot=False)
res_wiki = ov.single.pathway_enrichment(adata_micro, pathways_dict=pathway_dict_wiki, organism='Human', group_by='leiden_scVI', plot=False)

# Define the output folder
output_folder = r'D:\AI work\zhang tongmei\pathway\adata_microglia_leiden'

# Create the folder if it doesn't exist
if not os.path.exists(output_folder):
    os.makedirs(output_folder)

# Save the results to CSV files
res_BP.to_csv(os.path.join(output_folder, 'res_BP.csv'), index=False)
res_CC.to_csv(os.path.join(output_folder, 'res_CC.csv'), index=False)
res_MF.to_csv(os.path.join(output_folder, 'res_MF.csv'), index=False)
res_Rea.to_csv(os.path.join(output_folder, 'res_Rea.csv'), index=False)
res_wiki.to_csv(os.path.join(output_folder, 'res_wiki.csv'), index=False)

print(f"All pathway enrichment results saved to {output_folder}")
# Combine all results into one DataFrame
combined_res = pd.concat([res_BP, res_CC, res_MF, res_Rea, res_wiki], axis=0, ignore_index=True)

# Save the combined results to a single CSV file
combined_res.to_csv(os.path.join(output_folder, 'combined_pathway_results.csv'), index=False)

print(f"All pathway enrichment results saved to {output_folder}")
print(f"Combined results saved as 'combined_pathway_results.csv'")



import sys
#!{sys.executable} -m pip -q install palantir fa2
import os
import fa2
os.environ['R_HOME'] = "C:/Program Files/R/R-4.2.2"
os.environ['PATH'] += os.pathsep + "C:/Program Files/R/R-4.2.2/bin/x64"
import warnings
warnings.simplefilter(action='ignore', category=Warning)
import scanpy as sc
import scFates as scf
import rpy2.robjects as robjects
robjects.r('''
  print("Hello from R")
''')
sc.settings.verbosity = 3
sc.settings.logfile = sys.stdout

















