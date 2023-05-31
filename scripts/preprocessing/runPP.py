#!/usr/bin/env python
# coding: utf-8

import scanpy as sc
import scib
import warnings

warnings.filterwarnings('ignore')

import numpy as np

warnings.filterwarnings('ignore')

np.random.seed(2023)

def shufflePortion(arr, percentage, leaveUnknown = True): 
    x = np.array(arr)
    if leaveUnknown:
        annoCellsIdx = np.where(x != "unknown")[0]
    else:
        annoCellsIdx = x        
    
    shuf = np.random.choice(annoCellsIdx,  
                            round(arr.shape[0]*percentage/100), 
                            replace=False) 
    arr[np.sort(shuf)] = arr[shuf]
    return arr
  
def unassignPortion(arr, percentage): 
    shuf = np.random.choice(np.arange(arr.shape[0]),  
                            round(arr.shape[0]*percentage/100), 
                            replace=False) 
    arr[np.sort(shuf)] = 'unknown' 
    return arr





def runPP(inPath, outPath, hvg, batch, rout, scale, seurat, label, inputAnno):
    """
    params:
        inPath: path of the anndata object
        outPath: path of the preprocessed file to be written
        hvg: number of highly variable genes to use
        rout: set to true to save a Seurat object
        scale: set to true to activate scaling
        seurat: set to true to produce hvg list
    """

    adata = sc.read(inPath)
    hvgs = adata.var.index

    # remove HVG if already precomputed
    if 'highly_variable' in adata.var:
        del adata.var['highly_variable']

    if hvg > 500:
        print("Computing HVGs ...")
        if seurat:
            hvgs = scib.preprocessing.hvg_batch(
                adata,
                batch_key=batch,
                target_genes=hvg,
                adataOut=False
            )
        else:
            adata = scib.preprocessing.hvg_batch(
                adata,
                batch_key=batch,
                target_genes=hvg,
                adataOut=True
            )
    if scale:
        print("Scaling data ...")
        adata = scib.preprocessing.scale_batch(adata, batch)
        
    for i in inputAnno:
        print(i)
        if i == "original":
            print("making a copy of original anno to guide integration")
            adata.obs['original_'+label] = adata.obs[label].copy()
        elif (i.startswith("unknown_")) & ("shuffled_" not in i):
            unknown = i.split("_")[1]
            print("Hiding "+ str(unknown) +"% of annotations")
            adata.obs["unknown_"+str(unknown)+"_"+ label] = adata.obs[label].copy()
            adata.obs["unknown_"+str(unknown)+"_"+ label] = adata.obs["unknown_"+str(unknown)+"_"+ label].cat.add_categories(['unknown'])
            adata.obs["unknown_"+str(unknown)+"_"+ label] = unassignPortion(adata.obs["unknown_"+str(unknown)+"_"+ label].values, float(unknown))        
        elif (i.startswith("shuffled_")) & ("unknown" not in i):
            shuffled = i.split("_")[1]
            print("Suffling "+ str(shuffled) +"% of annotations")
            adata.obs["shuffled_"+str(shuffled)+"_"+ label] = adata.obs[label].copy()
            adata.obs["shuffled_"+str(shuffled)+"_"+ label] = shufflePortion(adata.obs["shuffled_"+str(shuffled)+"_"+ label].values, float(shuffled))
        elif (i.startswith("unknown_")) & ("shuffled_" in i):
            unknown = i.split("_")[1]
            shuffled = i.split("_")[3]
            print("Hiding "+ str(unknown) +"% of annotations and shuffling " + str(shuffled) + "% of remaining annotations")
            adata.obs["unknown_"+str(unknown)+"_"+ label] = adata.obs[label].copy()
            adata.obs["unknown_"+str(unknown)+"_"+ label] = adata.obs["unknown_"+str(unknown)+"_"+ label].cat.add_categories(['unknown'])
            adata.obs["unknown_"+str(unknown)+"_"+ label] = unassignPortion(adata.obs["unknown_"+str(unknown)+"_"+ label].values, float(unknown))
            #Shuffling
            alteredAnnoCol = "unknown_"+str(unknown)+"_shuffled_"+str(shuffled)+"_"+ label
            adata.obs[alteredAnnoCol] = adata.obs["unknown_"+str(unknown)+"_"+ label].copy()
            adata.obs[alteredAnnoCol] = shufflePortion(adata.obs[alteredAnnoCol].values, float(shuffled))
        else:
            otherAnno =  i+'_'+label
            if not any([j.startswith(i) for j in adata.obs.columns]):
                print(i +" annotations not available.")
                print("Adding "+i+" annotations equivalent to original annotations...")
                adata.obs['broad_'+label] = adata.obs[label].copy()

    if rout:
        print("Save as RDS")
        scib.preprocessing.saveSeurat(adata, outPath, batch, hvgs)

    else:
        print("Save as HDF5")
        sc.write(outPath, adata)


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Run the integration methods')

    parser.add_argument('-i', '--input_file', required=True)
    parser.add_argument('-o', '--output_file', required=True)
    parser.add_argument('-b', '--batch', required=True, help='Batch variable')
    parser.add_argument('-t', '--label', required=True, help='label variable for making noise')
    parser.add_argument('-v', '--hvgs', help='Number of highly variable genes', default=2000)
    parser.add_argument('-r', '--rout', help='Save output for R methods', action='store_true')
    parser.add_argument('-s', '--scale', action='store_true', help='Scale the data per batch')
    parser.add_argument('-a', '--input_annotations', help='all different annotations that will be use to guide supervised tools (separated by ,', default="original")
    parser.add_argument('-l', '--seurat', help='Generate output for seurat including hvg list', action='store_true')

    args = parser.parse_args()
    file = args.input_file
    out = args.output_file
    batch = args.batch
    label = args.label
    hvg = int(args.hvgs)
    rout = args.rout
    seurat = args.seurat
    scale = args.scale
    inputAnno = args.input_annotations.split(",")

    runPP(file, out, hvg, batch, rout, scale, seurat,label,inputAnno)
