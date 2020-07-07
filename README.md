# CSN
## @AUTHOR Gialitsis and Voutsadaki
### Implemented as project for the course Machine Learning in Computational Biology 2020

File Descriptions:   

        Already Implemented:
               csndm.m: this function performs the transformation from gene expression matrix to network degree matrix (ndm).
               csnedge.m: this function calculates the normalized statistic of edge x-y corresponding to two genes
               csnet.m: The function performs the transformation from gene expression matrix to cell-specific network (csn)
        Implemented By Us:
                buettner_analysis.r : this script performs the differential expression analysis on the gene expression data 
                all_in_one_CSN.ipynb : this jupyter notebook constructs and visualizes the representative cell specific network according to various centrality measures 
      
               
               


Libraries and Dependencies
python3.6
Install GeoParse(2.0.1): !pip install GEOparse:

        Requirement pandas>=0.17 (1.0.5)
        Requirement requests>=2.21.0 (2.23.0)
        Requirement numpy>=1.7 (1.18.5)
        Requirement tqdm>=4.31.1 (4.41.1)
        Requirement python-dateutil>=2.6.1 (2.8.1)
        Requirement pytz>=2017.2 (2018.9)
        Requirement idna<3,>=2.5 (2.9)
        Requirement certifi>=2017.4.17(2020.6.20)
        Requirement  chardet<4,>=3.0.2 (3.0.4)
        Requirement urllib3!=1.25.0,!=1.25.1,<1.26,>=1.21.1 (1.24.3)
        Requirement six>=1.5 (1.12.0)

Import Libraries:

    import pandas as pd
    import copy
    import matplotlib
    import matplotlib.pyplot as plt
    import matplotlib.colors as colors
    import matplotlib.cm as cm
    import numpy as np
    import networkx as nx
    import sys
    from scipy.optimize import curve_fit
    import numpy as np
    from scipy.stats import norm
    import GEOparse
    import pandas

Execution:

1. Unzip Data folder
2. Open notebook all_in_one_CSN.ipynb (we used Google Collaboratory)
3. Import files ReprCSN.csv and gea_norm.csv in the working directory
4. Runtime -> Run All

  






' Forked Readme Below (ownership belongs to the owners of the respective authors)'
Cell-specific Network Constructed by Single-cell RNA Sequencing Data

    function ndm = csndm(data,alpha,boxsize,normalize)

 Construction of network degree matrix
 
 The function performs the transformation from gene expression matrix to network degree matrix (ndm).
 
 data: Gene expression matrix (TPM/FPKM/RPKM/count), rows = genes, columns = cells
 
 alpha: Significant level (eg. 0.001, 0.01, 0.05 ...), Default = 0.01
 
 boxsize: Size of neighborhood, Default = 0.1
 
 normalize: 1: result is normalized (Default); 0: result is not normalized
 
    
    
    
 
    function csn = csnet(data,c,alpha,boxsize,weighted)
 Construction of cell-specific network
 
 The function performs the transformation from gene expression matrix to cell-specific network (csn).
 
 data: Gene expression matrix, rows = genes, columns = cells
 
 c: Construct the CSNs for all cells, set c = [] (Default); Construct the CSN for cell k, set c = k
 
 alpha: Significant level (eg. 0.001, 0.01, 0.05 ...), larger alpha leads to more edges, Default = 0.01
 
 boxsize: Size of neighborhood, Default = 0.1
 
 weighted: 1: edge is weighted; 0: edge is not weighted (Default)
 
 csn: Cell-specific network, the kth CSN is in csn{k}, rows = genes, columns = genes
 
 Note that too many cells or genes may lead to out of memory.
 
    
    
    
 
    function edge = csnedge(gx,gy,boxsize)

 The normalized statistic of edge x-y
 
 gx gy: Gene expression values of gene x and gene y. If there are n cells, gx and gy are 1-by-n vectors.
 
 boxsize: Size of neighborhood, Default = 0.1
 
 edge: 1-by-n vector, the normalized statistic of edge x-y in all cells
