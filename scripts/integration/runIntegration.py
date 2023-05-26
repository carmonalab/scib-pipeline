#!/usr/bin/env python
# coding: utf-8

import scanpy as sc
import scib
import warnings

warnings.filterwarnings('ignore')

def scgen(adata, batch, dims, cell_type, epochs=100, hvg=None, **kwargs):
    """scGen wrapper function
    Adapted from scib package version 1.1.3, based on `scgen`_ version 2.1.0 with parametrization taken from the tutorial `notebook`_.
    .. _scgen: https://github.com/theislab/scgen
    .. _notebook: https://scgen.readthedocs.io/en/stable/tutorials/scgen_batch_removal.html
    :param adata: preprocessed ``anndata`` object
    :param batch: batch key in ``adata.obs``
    :param dims: number of dimensions to use for integration
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
    model = SCGEN(net_adata,n_latent = dims)  
    model.train(
        max_epochs=epochs,
        batch_size=32,
        early_stopping=True,
        early_stopping_patience=25,
    )
    corrected_adata = model.batch_removal(**kwargs)
    return corrected_adata
  

def scanorama(adata, batch, dims, hvg = None, **kwargs):
    """Scanorama wrapper function
    Adapted from scib package version 1.1.3, based on `scanorama <https://github.com/brianhie/scanorama>`_ version 1.7.0
    :param adata: preprocessed ``anndata`` object
    :param batch: batch key in ``adata.obs``
    :param dims: number of dimensions to use for integration
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
    corrected = scanorama.correct_scanpy(split, return_dimred=True, dimred = dims, **kwargs) 
    corrected = scib.utils.merge_adata(
        *corrected, batch_key=batch, batch_categories=categories, index_unique=None
    )
    corrected.obsm["X_emb"] = corrected.obsm["X_scanorama"]
    # corrected.uns['emb']=True

    return corrected
  

def scvi(adata, batch, dims, hvg=None, return_model=False, max_epochs=None):
    """scVI wrapper function

    Adapted from scib package version 1.1.3, based on scvi-tools version >=0.16.0 (available through `conda <https://docs.scvi-tools.org/en/stable/installation.html>`_)

    .. note::
        scVI expects only non-normalized (count) data on highly variable genes!

    :param adata: preprocessed ``anndata`` object
    :param batch: batch key in ``adata.obs``
    :param dims: number of dimensions to use for integration
    :param hvg: list of highly variables to subset to. If ``None``, the full dataset will be used
    :return: ``anndata`` object containing the corrected feature matrix as well as an embedding representation of the
        corrected data
    """
    import numpy as np
    try:
        from scvi.model import SCVI
    except ModuleNotFoundError as e:
        raise OptionalDependencyNotInstalled(e)

    scib.utils.check_sanity(adata, batch, hvg)

    # Check for counts data layer
    if "counts" not in adata.layers:
        raise TypeError(
            "Adata does not contain a `counts` layer in `adata.layers[`counts`]`"
        )
        
    n_latent = dims #Defaults from SCVI github tutorials scanpy_pbmc3k and harmonization is 30
    
    # Defaults from SCVI github tutorials scanpy_pbmc3k and harmonization
    n_hidden = 128
    n_layers = 2

    # copying to not return values added to adata during setup_anndata
    net_adata = adata.copy()
    if hvg is not None:
        net_adata = adata[:, hvg].copy()
    SCVI.setup_anndata(net_adata, layer="counts", batch_key=batch)

    vae = SCVI(
        net_adata,
        gene_likelihood="nb",
        n_layers=n_layers,
        n_latent=n_latent,
        n_hidden=n_hidden,
    )
    train_kwargs = {"train_size": 1.0}
    if max_epochs is not None:
        train_kwargs["max_epochs"] = max_epochs
    vae.train(**train_kwargs)
    adata.obsm["X_emb"] = vae.get_latent_representation()

    if not return_model:
        return adata
    else:
        return vae
      

def scanvi(adata, batch, dims, labels, hvg=None, max_epochs=None):
    """scANVI wrapper function

    Adapted from scib package version 1.1.3, based on scvi-tools version >=0.16.0 (available through `conda <https://docs.scvi-tools.org/en/stable/installation.html>`_)

    .. note::
        Use non-normalized (count) data for scANVI!

    :param adata: preprocessed ``anndata`` object
    :param batch: batch key in ``adata.obs``
    :param dims: number of dimensions to use for integration
    :param labels: label key in ``adata.obs``
    :param hvg: list of highly variables to subset to. If ``None``, the full dataset will be used
    :return: ``anndata`` object containing the corrected feature matrix as well as an embedding representation of the
        corrected data
    """
    import numpy as np
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

    vae = scvi(adata, batch, dims, hvg, return_model=True,max_epochs=n_epochs_scVI)

    # STEP 2: RUN scVI to initialize scANVI
    scanvae = SCANVI.from_scvi_model(
        scvi_model=vae,
        labels_key=labels,
        unlabeled_category="unknown",  # updated to match STACAS keyword
    )
    scanvae.train(max_epochs=n_epochs_scANVI, train_size=1.0)
    adata.obsm["X_emb"] = scanvae.get_latent_representation()

    return adata

def combat(adata, batch, dims):
    """ComBat wrapper function (``scanpy`` implementation)

    Adapted from scib package version 1.1.3, using scanpy implementation of `Combat <https://scanpy.readthedocs.io/en/stable/generated/scanpy.pp.combat.html>`_

    :param adata: preprocessed ``anndata`` object
    :param batch: batch key in ``adata.obs``
    :param dims: number of dimensions to use for integration (not used by combat)
    :return: ``anndata`` object containing the corrected feature matrix
    """
    adata_int = adata.copy()
    sc.pp.combat(adata_int, key=batch)
    return adata_int

  

def runIntegration(inPath, outPath, method, hvg, dims, batch, celltype=None):
    """
    params:
        method: name of method
        batch: name of `adata.obs` column of the batch
        max_genes_hvg: maximum number of HVG
    """

    adata = sc.read(inPath)

    if celltype is not None:
        integrated = method(adata, batch, dims, celltype)
    else:
        integrated = method(adata, batch, dims)

    sc.write(outPath, integrated)


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Run the integration methods')

    parser.add_argument('-m', '--method', required=True)
    parser.add_argument('-i', '--input_file', required=True)
    parser.add_argument('-o', '--output_file', required=True)
    parser.add_argument('-b', '--batch', required=True, help='Batch variable')
    parser.add_argument('-n', '--latent_dims', default=30)
    parser.add_argument('-v', '--hvgs', help='Number of highly variable genes', default=2000)
    parser.add_argument("-c", '--celltype', help='Cell type variable', default=None)

    args = parser.parse_args()
    file = args.input_file
    out = args.output_file
    batch = args.batch
    hvg = int(args.hvgs)
    dims = int(args.latent_dims)
    celltype = args.celltype
    method = args.method
    methods = {
        'scanorama': scanorama,
        #'trvae': scib.integration.trvae, not included in our pipeline
        #'trvaep': scib.integration.trvaep, not included in our pipeline
        'scgen': scgen,
        #'mnn': scib.integration.mnn, not included in our pipeline
        #'bbknn': scib.integration.bbknn, not included in our pipeline
        'scvi': scvi,
        'scanvi': scanvi,
        'combat': combat,
        #'saucie': scib.integration.saucie, not included in our pipeline
        #'desc': scib.integration.desc, not included in our pipeline
    }

    if method not in methods.keys():
        raise ValueError(f'Method "{method}" does not exist. Please use one of '
                         f'the following:\n{list(methods.keys())}')

    run = methods[method]
    runIntegration(file, out, run, hvg, dims, batch, celltype)
