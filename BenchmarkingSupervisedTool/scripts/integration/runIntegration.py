#!/usr/bin/env python
# coding: utf-8

import scanpy as sc
import scib
import warnings
import numpy as np

warnings.filterwarnings('ignore')

def scanvi(adata, batch, labels, hvg=None, max_epochs=None):
    """scANVI wrapper function
    Based on scvi-tools version >=0.16.0 (available through `conda <https://docs.scvi-tools.org/en/stable/installation.html>`_)
    .. note::
        Use non-normalized (count) data for scANVI!
    :param adata: preprocessed ``anndata`` object
    :param batch: batch key in ``adata.obs``
    :param labels: label key in ``adata.obs``
    :param hvg: list of highly variables to subset to. If ``None``, the full dataset will be used
    :return: ``anndata`` object containing the corrected feature matrix as well as an embedding representation of the
        corrected data
    """
    try:
        from scvi.model import SCANVI
    except ModuleNotFoundError as e:
        raise OptionalDependencyNotInstalled(e)

    # # Defaults from SCVI github tutorials scanpy_pbmc3k and harmonization
    # this n_epochs_scVI is now default in scvi-tools
    if max_epochs is None:
        n_epochs_scVI = int(np.min([round((20000 / adata.n_obs) * 400), 400]))  # 400
        n_epochs_scANVI = int(np.min([10, np.max([2, round(n_epochs_scVI / 3.0)])]))
    else:
        n_epochs_scVI = max_epochs
        n_epochs_scANVI = max_epochs

    vae = scib.integration.scvi(adata, batch, hvg, return_model=True, max_epochs=n_epochs_scVI)

    # STEP 2: RUN scVI to initialize scANVI
    scanvae = SCANVI.from_scvi_model(
        scvi_model=vae,
        labels_key=labels,
        unlabeled_category="unknown",  # pick anything definitely not in a dataset
    )
    scanvae.train(max_epochs=n_epochs_scANVI, train_size=1.0)
    adata.obsm["X_emb"] = scanvae.get_latent_representation()

    return adata
  

def scgen(adata, batch, cell_type, epochs=100, latent_dims = 30, hvg=None, **kwargs):
    """scGen wrapper function
    Based on `scgen`_ version 2.1.0 with parametrization taken from the tutorial `notebook`_.
    with adapted default latent space size
    .. _scgen: https://github.com/theislab/scgen
    .. _notebook: https://scgen.readthedocs.io/en/stable/tutorials/scgen_batch_removal.html
    :param adata: preprocessed ``anndata`` object
    :param batch: batch key in ``adata.obs``
    :param hvg: list of highly variables to subset to. If ``None``, the full dataset will be used
    :return: ``anndata`` object containing the corrected feature matrix
    """
    try:
        from scgen import SCGEN
    except ModuleNotFoundError as e:
        raise OptionalDependencyNotInstalled(e)

    scib.utils.check_sanity(adata, batch, hvg)

    net_adata = adata.copy()
    if hvg is not None:
        net_adata = net_adata[:, hvg].copy()

    SCGEN.setup_anndata(net_adata, batch_key=batch, labels_key=cell_type)
    model = SCGEN(net_adata,n_latent = latent_dims)  # adapted n_latent for consitency with other methods
    model.train(
        max_epochs=epochs,
        batch_size=32,
        early_stopping=True,
        early_stopping_patience=25,
    )
    corrected_adata = model.batch_removal(**kwargs)
    return corrected_adata
  

def scanorama(adata, batch, hvg = None,latent_dims = 30, **kwargs):
    """Scanorama wrapper function
    Based on `scanorama <https://github.com/brianhie/scanorama>`_ version 1.7.0
    with adapted default latent space size
    :param adata: preprocessed ``anndata`` object
    :param batch: batch key in ``adata.obs``
    :param hvg: list of highly variables to subset to. If ``None``, the full dataset will be used
    :return: ``anndata`` object containing the corrected feature matrix as well as an embedding representation of the
        corrected data
    """
    try:
        import scanorama
    except ModuleNotFoundError as e:
        raise OptionalDependencyNotInstalled(e)

    scib.utils.check_sanity(adata, batch, hvg)
    split, categories = scib.utils.split_batches(adata.copy(), batch, return_categories=True)
    corrected = scanorama.correct_scanpy(split, return_dimred=True, dimred = latent_dims, **kwargs) # adapted dimred for consitency with other methods
    corrected = scib.utils.merge_adata(
        *corrected, batch_key=batch, batch_categories=categories, index_unique=None
    )
    corrected.obsm["X_emb"] = corrected.obsm["X_scanorama"]
    # corrected.uns['emb']=True

    return corrected


def runIntegration(inPath, outPath, method, hvg, batch, celltype=None):
    """
    params:
        method: name of method
        batch: name of `adata.obs` column of the batch
        max_genes_hvg: maximum number of HVG
    """

    adata = sc.read(inPath)

    if celltype is not None:
        integrated = method(adata, batch, celltype)
    else:
        integrated = method(adata, batch)

    sc.write(outPath, integrated)


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Run the integration methods')

    parser.add_argument('-m', '--method', required=True)
    parser.add_argument('-i', '--input_file', required=True)
    parser.add_argument('-o', '--output_file', required=True)
    parser.add_argument('-b', '--batch', required=True, help='Batch variable')
    parser.add_argument('-v', '--hvgs', help='Number of highly variable genes', default=2000)
    parser.add_argument("-c", '--celltype', help='Cell type variable', default=None)

    args = parser.parse_args()
    file = args.input_file
    out = args.output_file
    batch = args.batch
    hvg = int(args.hvgs)
    celltype = args.celltype
    method = args.method
    methods = {
        'scanorama': scanorama,
        'trvae': scib.integration.trvae,
        'trvaep': scib.integration.trvaep,
        'scgen': scgen,
        'mnn': scib.integration.mnn,
        'bbknn': scib.integration.bbknn,
        'scvi': scib.integration.scvi,
        'scanvi': scanvi,
        'combat': scib.integration.combat,
        'saucie': scib.integration.saucie,
        'desc': scib.integration.desc
    }

    if method not in methods.keys():
        raise ValueError(f'Method "{method}" does not exist. Please use one of '
                         f'the following:\n{list(methods.keys())}')

    run = methods[method]
    runIntegration(file, out, run, hvg, batch, celltype)
