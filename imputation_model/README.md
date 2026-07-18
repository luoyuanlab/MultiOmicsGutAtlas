# GCN-Based Gene Expression Imputation for CosMx Spatial Transcriptomics

A graph convolutional network (GCN) model that imputes full-transcriptome gene expression in lower-plex CosMx spatial transcriptomics panels using a whole-transcriptome (WTX) reference. The model constructs mutual nearest neighbor (MNN) graphs between reference and query cells and trains a two-layer graph convolutional autoencoder to predict held-out gene expression via weighted MSE loss.

**Primary use case:** imputing the CosMx 1K panel to full transcriptome scale, enabling cell type annotation and downstream spatial analyses that require broader gene coverage. The demo uses a 6K panel query as a faster example (fewer genes to impute).

---

## System Requirements

### Software Dependencies

- Python ≥ 3.9
- PyTorch ≥ 2.0 (with CUDA support recommended)
- scanpy ≥ 1.9
- anndata ≥ 0.9
- scipy ≥ 1.10
- numpy ≥ 1.24
- pandas ≥ 1.5
- scikit-learn ≥ 1.2
- tqdm
- matplotlib
- seaborn

Install all dependencies via pip:

```bash
pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu118
pip install scanpy anndata scipy numpy pandas scikit-learn tqdm matplotlib seaborn
```

### Tested On

- Ubuntu 20.04 / CentOS 8
- Python 3.10
- PyTorch 2.1 with CUDA 11.8
- NVIDIA A100 / V100 GPUs

### Hardware

- **GPU strongly recommended.** Training on CPU is possible but will be slow for large datasets.
- Minimum 16 GB GPU memory recommended for datasets with >30,000 reference cells.
- For the demo dataset (~30,000 reference + ~57,000 query cells): ≥24 GB GPU memory.

---

## Installation Guide

1. Clone or download this repository.

2. Install dependencies (see above). A conda environment is recommended:

```bash
conda create -n gcn_imputation python=3.10
conda activate gcn_imputation
pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu118
pip install scanpy anndata scipy numpy pandas scikit-learn tqdm matplotlib seaborn jupyter
```

3. No compilation or build steps are required. The model is imported directly from the `Imputation_model/` package.

**Typical install time:** 5–10 minutes on a standard desktop with internet access.

---

## Demo

### Data

Demo data (UC WTX reference and 6K query datasets) are available on Figshare:
https://figshare.com/articles/dataset/33023585

Download the following files and update the paths in `Tutorial.ipynb`:
- `ref_PATH`: WTX reference dataset (`.h5ad` or `.csv`, cells × genes)
- `query_PATH`: CosMx 6K panel query dataset (`.h5ad` or `.csv`, cells × genes)

### Running the Demo

Open and run `Tutorial.ipynb`:

```bash
jupyter notebook Tutorial.ipynb
```

The notebook walks through the full pipeline:
1. Load and preprocess reference (WTX) and query (panel) datasets
2. Compute MNN-based intra- and inter-dataset adjacency matrices
3. Train the GCN model
4. Save imputed gene expression to CSV

### Expected Output

Two CSV files saved to `save_PATH`:
- `Gene_imputation_best_log1p.csv` — imputed expression at best validation epoch (log1p-normalized)
- `Gene_imputation_final_log1p.csv` — imputed expression at final epoch (log1p-normalized)

Rows = query cells, columns = imputed genes (genes present in reference but not in query panel).

### Expected Run Time (Demo)

With 30,000 reference cells and 57,000 query cells (~6,000 input genes, ~13,000 predicted genes):

| Setting | Epochs | Time |
|---------|--------|------|
| GPU (A100) | 600 | ~80 minutes |
| GPU (A100) | 10,000 (recommended) | ~22 hours |
| CPU | 600 | several hours (not recommended) |

---

## Instructions for Use

### Running on Your Own Data

Edit the path and parameter block at the top of `Tutorial.ipynb`:

```python
ref_PATH          = "/path/to/your/wtx_reference.h5ad"   # WTX reference (.h5ad or .csv)
query_PATH        = "/path/to/your/panel_query.h5ad"      # Panel query (.h5ad or .csv)
save_PATH         = "/path/to/save/results"
ref_subsample_num = 30000   # Set to None to use all reference cells (increases memory/time)
scale_type        = 'gene'  # 'gene', 'cell', or False
log1p             = True

device  = torch.device('cuda:0')  # Change to 'cpu' if no GPU available
epochs  = 10000                   # Recommended ≥10,000
lr      = 0.00001
dropout = 0.0
```

### Input Format

- **Reference** (`ref_PATH`): whole-transcriptome scRNA-seq or CosMx WTX data, cells × genes, raw or normalized counts. Accepts `.h5ad` or `.csv`.
- **Query** (`query_PATH`): lower-plex CosMx panel data (1K or 6K), cells × genes, same normalization as reference. Accepts `.h5ad` or `.csv`.
- Genes are matched by name; overlapping genes are used for graph construction, unique reference genes are imputed.

### Tips

- Use `ref_subsample_num` to subsample the reference if GPU memory is limited.
- Set `scale_type='gene'` (default) for gene-wise max normalization, which generally gives the best performance.
- The best model weights (lowest validation loss) are saved automatically to `save_PATH/model_best_wts.pt` and can be reloaded:

```python
from Imputation_model.GAE import CosMxGAE
model = CosMxGAE(adj_RNA, input_dim, output_dim, hidden_dim, dropout)
model.load_state_dict(torch.load(f'{save_PATH}/model_best_wts.pt'))
model.eval()
```

### Output

Imputed expression values are log1p-normalized and clipped to non-negative values. To recover linear-scale counts: `np.expm1(imputed_log1p)`.
