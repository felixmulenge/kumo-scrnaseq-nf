#!/usr/bin/env python3
import tempfile
import shutil
import os
import argparse
import scanpy as sc
import matplotlib.pyplot as plt
from scib_metrics.benchmark import Benchmarker, BioConservation, BatchCorrection


def main(args):
  
  #print(f"Loading AnnData: {args.input_h5ad}")
  adata = sc.read_h5ad(args.input_h5ad)

  embeddings = args.embeddings.split(",")
  print(f"Embeddings to benchmark: {embeddings}")

# -------------------------
# Initialize Benchmarker
# -------------------------
  bm = Benchmarker(
    adata,
    batch_key=args.batch_key,
    label_key=args.label_key,
    embedding_obsm_keys=embeddings,
    bio_conservation_metrics=BioConservation(),
    batch_correction_metrics=BatchCorrection(),
    n_jobs=args.n_jobs,
  )

  # -------------------------
  # Run benchmarking
  # -------------------------
  print("Running scIB benchmark...")
  bm.benchmark()

  # -------------------------
  # Save results tables
  # -------------------------
  raw_df = bm.get_results(min_max_scale=False)
  raw_df = raw_df.transpose()
  raw_out = f"{args.output_prefix}_scib_metrics_raw.tsv"
  raw_df.to_csv(raw_out, sep="\t")
  print(f"Saved raw metrics table: {raw_out}")

  scaled_df = bm.get_results(min_max_scale=True)
  scaled_df = scaled_df.transpose()
  scaled_out = f"{args.output_prefix}_scib_metrics_scaled.tsv"
  scaled_df.to_csv(scaled_out, sep="\t")
  print(f"Saved scaled metrics table: {scaled_out}")


  # -------------------------
  # Plot: scaled metrics
  # -------------------------
  with tempfile.TemporaryDirectory() as tmpdir:
      bm.plot_results_table(min_max_scale=False,show=False,save_dir=tmpdir)
      src = os.path.join(tmpdir,"scib_results.svg")
      dst =f"{args.output_prefix}_scib_metrics_raw.svg"
      shutil.move(src,dst)

  # -------------------------
  # Plot: unscaled metrics
  # -------------------------
  with tempfile.TemporaryDirectory() as tmpdir:
      bm.plot_results_table(min_max_scale=True,show=False,save_dir=tmpdir)
      src = os.path.join(tmpdir,"scib_results.svg")
      dst =f"{args.output_prefix}_scib_metrics_scaled.svg"
      shutil.move(src,dst)


  print("scIB benchmarking completed successfully!")


if __name__ == "__main__":
  
  parser = argparse.ArgumentParser(description="Run scIB integration benchmarking on an h5ad file")
  parser.add_argument("--input_h5ad",required=True,help="Input AnnData (.h5ad)")
  parser.add_argument("--output_prefix",required=True,help="Prefix for all output files")
  parser.add_argument("--batch_key",default="batch",help="Batch column in adata.obs (default: batch)")
  parser.add_argument("--label_key",default="cell_type",help="Cell type column in adata.obs (default: cell_type)")
  parser.add_argument("--embeddings",default=" ",help="Comma-separated embedding keys in adata.obsm")
  parser.add_argument("--n_jobs",type=int,default=1,help="Number of parallel jobs (default: 1)")
  
  args = parser.parse_args()
  main(args)




#python run_scib_benchmark.py \
#--input_h5ad lung_integrated.h5ad \
#--output_prefix lung \
#--batch_key batch \
#--label_key cell_type \
#--embeddings X_pca,X_scVI,X_harmony \
#--n_jobs 8
