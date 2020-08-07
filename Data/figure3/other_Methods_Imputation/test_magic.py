import magic
import scanpy as sc

adata = sc.read_10x_mtx("/home/math/hyl2016/Imputation/all_data_5W/")
magic_operator = magic.MAGIC(n_pca =20)
X_magic = magic_operator.fit_transform(adata, genes="all_genes")

print("MAGIC Finished")