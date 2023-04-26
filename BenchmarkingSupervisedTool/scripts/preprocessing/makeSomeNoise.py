#!/usr/bin/env python
# coding: utf-8

import scanpy as sc
import scib
import warnings
import numpy as np

warnings.filterwarnings('ignore')

np.random.seed(2023)

def shuffle_portion(arr, percentage): 
    shuf = np.random.choice(np.arange(arr.shape[0]),  
                            round(arr.shape[0]*percentage/100), 
                            replace=False) 
    arr[np.sort(shuf)] = arr[shuf] 
    return arr
  



def makeNoise(inPath, outPath, label, batch, percentNoise):
    """
    params:
        inPath: path of the anndata object
        outPath: path of the processed file to be written
        label: cell type label column
        percentNoise: percentage of noise
    """

    adata = sc.read(inPath)
    for (n in [10,20,50]):
      adata.obs["noisy_"+n+"_"+ label] = adata.obs[label].copy()
      adata.obs["noisy_"+n+"_"+ label] = shuffle_portion(adata.obs["noisy_"+n+"_"+ label].values, n)
    
    
    #adata.obs['correct'] = adata.obs["noisy_" + label] == adata.obs[label]
    #adata.obs["labelling"] = ["correct" if c else "wrong" for c in list(adata.obs['correct'])]
    sc.write(outPath, adata)
    scib.preprocessing.saveSeurat(adata, outPath, batch)



if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Run the integration methods')

    parser.add_argument('-i', '--input_file', required=True)
    parser.add_argument('-o', '--output_file', required=True)
    parser.add_argument('-s', '--output_file_seurat', required=True)
    parser.add_argument('-l', '--label', required=True, help='label variable')
    parser.add_argument('-b', '--batch', required = True, help='batch variable')
    #parser.add_argument('-p', '--percent', help='Noise (shuffling) percentage', default=20)
  
    args = parser.parse_args()
    file = args.input_file
    out = args.output_file
    outSeurat = = args.output_file_seurat
    label = args.label
    batch = args.batch
    #percent = float(args.percent)

    makeNoise(file, out, label, batch, percent)
