from pathlib import Path
from urllib.request import Request, urlopen


def download_file(url: str, dest: Path | str, chunk_size: int = 8192):
    """Download a file using urllib with a User-Agent header."""
    req = Request(
        url, headers={"User-Agent": "Mozilla/5.0 (compatible; Python urllib)"}
    )
    with urlopen(req) as response, open(dest, "wb") as out_file:
        while True:
            chunk = response.read(chunk_size)
            if not chunk:
                break
            out_file.write(chunk)


def setup_dandelion_tutorial_bcr(path: Path | str | None = None) -> None:
    """Download example BCR datasets for Dandelion tutorial."""
    base = Path("./dandelion_tutorial") if path is None else Path(path)
    base.mkdir(parents=True, exist_ok=True)

    datasets = {
        "vdj_v1_hs_pbmc3": {
            "filtered_feature_bc_matrix.h5": "https://cf.10xgenomics.com/samples/cell-vdj/3.1.0/vdj_v1_hs_pbmc3/vdj_v1_hs_pbmc3_filtered_feature_bc_matrix.h5",
            "filtered_contig_annotations.csv": "https://cf.10xgenomics.com/samples/cell-vdj/3.1.0/vdj_v1_hs_pbmc3/vdj_v1_hs_pbmc3_b_filtered_contig_annotations.csv",
            "filtered_contig.fasta": "https://cf.10xgenomics.com/samples/cell-vdj/3.1.0/vdj_v1_hs_pbmc3/vdj_v1_hs_pbmc3_b_filtered_contig.fasta",
        },
        "vdj_nextgem_hs_pbmc3": {
            "filtered_feature_bc_matrix.h5": "https://cf.10xgenomics.com/samples/cell-vdj/3.1.0/vdj_nextgem_hs_pbmc3/vdj_nextgem_hs_pbmc3_filtered_feature_bc_matrix.h5",
            "filtered_contig_annotations.csv": "https://cf.10xgenomics.com/samples/cell-vdj/3.1.0/vdj_nextgem_hs_pbmc3/vdj_nextgem_hs_pbmc3_b_filtered_contig_annotations.csv",
            "filtered_contig.fasta": "https://cf.10xgenomics.com/samples/cell-vdj/3.1.0/vdj_nextgem_hs_pbmc3/vdj_nextgem_hs_pbmc3_b_filtered_contig.fasta",
        },
        "sc5p_v2_hs_PBMC_10k": {
            "filtered_feature_bc_matrix.h5": "https://cf.10xgenomics.com/samples/cell-vdj/4.0.0/sc5p_v2_hs_PBMC_10k/sc5p_v2_hs_PBMC_10k_filtered_feature_bc_matrix.h5",
            "filtered_contig_annotations.csv": "https://cf.10xgenomics.com/samples/cell-vdj/4.0.0/sc5p_v2_hs_PBMC_10k/sc5p_v2_hs_PBMC_10k_b_filtered_contig_annotations.csv",
            "filtered_contig.fasta": "https://cf.10xgenomics.com/samples/cell-vdj/4.0.0/sc5p_v2_hs_PBMC_10k/sc5p_v2_hs_PBMC_10k_b_filtered_contig.fasta",
        },
        "sc5p_v2_hs_PBMC_1k": {
            "filtered_feature_bc_matrix.h5": "https://cf.10xgenomics.com/samples/cell-vdj/4.0.0/sc5p_v2_hs_PBMC_1k/sc5p_v2_hs_PBMC_1k_filtered_feature_bc_matrix.h5",
            "filtered_contig_annotations.csv": "https://cf.10xgenomics.com/samples/cell-vdj/4.0.0/sc5p_v2_hs_PBMC_1k/sc5p_v2_hs_PBMC_1k_b_filtered_contig_annotations.csv",
            "filtered_contig.fasta": "https://cf.10xgenomics.com/samples/cell-vdj/4.0.0/sc5p_v2_hs_PBMC_1k/sc5p_v2_hs_PBMC_1k_b_filtered_contig.fasta",
        },
    }

    for dirname, files in datasets.items():
        dirpath = base / dirname
        dirpath.mkdir(parents=True, exist_ok=True)

        for filename, url in files.items():
            outfile = dirpath / filename
            if outfile.exists():
                continue
            print(f"Downloading {filename} → {outfile}")
            download_file(url, outfile)


def setup_dandelion_tutorial_tcr(path: Path | str | None = None) -> None:
    """Download example TCR datasets for Dandelion tutorial."""
    base = Path("./dandelion_tutorial") if path is None else Path(path)
    base.mkdir(parents=True, exist_ok=True)

    datasets = {
        "sc5p_v2_hs_PBMC_10k": {
            "filtered_feature_bc_matrix.h5": "https://cf.10xgenomics.com/samples/cell-vdj/4.0.0/sc5p_v2_hs_PBMC_10k/sc5p_v2_hs_PBMC_10k_filtered_feature_bc_matrix.h5",
            "airr_rearrangement.tsv": "https://cf.10xgenomics.com/samples/cell-vdj/4.0.0/sc5p_v2_hs_PBMC_10k/sc5p_v2_hs_PBMC_10k_t_airr_rearrangement.tsv",
            "filtered_contig_annotations.csv": "https://cf.10xgenomics.com/samples/cell-vdj/4.0.0/sc5p_v2_hs_PBMC_10k/sc5p_v2_hs_PBMC_10k_t_filtered_contig_annotations.csv",
            "filtered_contig.fasta": "https://cf.10xgenomics.com/samples/cell-vdj/4.0.0/sc5p_v2_hs_PBMC_10k/sc5p_v2_hs_PBMC_10k_t_filtered_contig.fasta",
        },
        "sc5p_v1p1_hs_melanoma_10k": {
            "filtered_feature_bc_matrix.h5": "https://cf.10xgenomics.com/samples/cell-vdj/4.0.0/sc5p_v1p1_hs_melanoma_10k/sc5p_v1p1_hs_melanoma_10k_filtered_feature_bc_matrix.h5",
            "airr_rearrangement.tsv": "https://cf.10xgenomics.com/samples/cell-vdj/4.0.0/sc5p_v1p1_hs_melanoma_10k/sc5p_v1p1_hs_melanoma_10k_t_airr_rearrangement.tsv",
            "filtered_contig_annotations.csv": "https://cf.10xgenomics.com/samples/cell-vdj/4.0.0/sc5p_v1p1_hs_melanoma_10k/sc5p_v1p1_hs_melanoma_10k_t_filtered_contig_annotations.csv",
            "filtered_contig.fasta": "https://cf.10xgenomics.com/samples/cell-vdj/4.0.0/sc5p_v1p1_hs_melanoma_10k/sc5p_v1p1_hs_melanoma_10k_t_filtered_contig.fasta",
        },
    }

    for dirname, files in datasets.items():
        dirpath = base / dirname
        dirpath.mkdir(parents=True, exist_ok=True)

        for filename, url in files.items():
            outfile = dirpath / filename
            if outfile.exists():
                continue
            print(f"Downloading {filename} → {outfile}")
            download_file(url, outfile)
