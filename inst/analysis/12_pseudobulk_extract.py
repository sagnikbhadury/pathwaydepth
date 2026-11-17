"""Pseudo-bulk expression of the DeepSpot predictions: the mean of the predicted
expression over all spots of a slide, for each of the predicted genes.

The DeepSpot files are not distributed with the package (CC BY-NC-ND 4.0); they
are obtained from their authors. The output is a slides x genes matrix that is
kept locally.

Usage: python 12_pseudobulk_extract.py <folder with TCGA_LUAD/FFPE and TCGA_LUSC/FFPE> <output .csv.gz>
"""
import os, sys, time
import numpy as np, pandas as pd, h5py

root, out = sys.argv[1], sys.argv[2]
rows, names, genes, n_spots = [], [], None, []
t0 = time.time()
files = [(sub, f) for sub in ("TCGA_LUAD/FFPE", "TCGA_LUSC/FFPE")
         for f in sorted(os.listdir(os.path.join(root, sub))) if f.endswith(".h5ad")]
for k, (sub, f) in enumerate(files):
    with h5py.File(os.path.join(root, sub, f), "r") as h:
        X = h["X"]
        assert not isinstance(X, h5py.Group), f          # dense in the DeepSpot release
        g = [v.decode() if isinstance(v, bytes) else str(v) for v in h["var"]["_index"][:]]
        if genes is None:
            genes = g
        assert g == genes, f                             # same gene order in every file
        m = np.zeros(X.shape[1], dtype=np.float64)
        for i in range(0, X.shape[0], 2000):              # chunked, float64 accumulation
            m += X[i:i + 2000].astype(np.float64).sum(axis=0)
        rows.append(m / X.shape[0]); names.append(f[:-5]); n_spots.append(X.shape[0])
    if k % 50 == 0:
        print(f"{k + 1}/{len(files)} {f} spots={X.shape[0]} elapsed={time.time() - t0:.0f}s", flush=True)
df = pd.DataFrame(np.vstack(rows), index=names, columns=genes)
df.insert(0, "n_spots", n_spots)
df.to_csv(out, compression="gzip")
print("written", out, df.shape, f"in {time.time() - t0:.0f}s")
