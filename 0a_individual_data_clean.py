# Process and Assemble the Selected Microarray Datasets into Single Matrix
# Copyright (C) 2023-4 Y. David Chen

import pandas as pd
import numpy as np
from os import path
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler

## Reusable Constants:
DIR = "*************PATH MASKED*************"
DIR_MATS = path.join(DIR, "geo_serial_matrices/expr_mats/")
PATH_META = path.join(DIR, "series_master.xlsx")
AFFY_CTRLS = ["Control sequence", "Exemplar sequence"]
STRNA = ["", "---", "NA", "nan", "N/A", "n/a", "NaN", "NAN"]
CHROMS = ["chr"+str(s) for s in np.arange(1,22+1, 1)] #autosomes ONLY
PERC_VAR_FILTER = 0.005 #per-dataset variance filter

## Reusable Functions:
def read_txt_matrix(accession, prefix=DIR_MATS, suffix=".txt"):
    mat = None
    if suffix == ".csv":
        mat = pd.read_csv(path.join(prefix, accession+suffix), index_col=0, low_memory=False)
    else:
        compression = "gzip" if suffix == ".txt.gz" else None
        mat = pd.read_table(path.join(prefix, accession+suffix), index_col=0, 
                            compression=compression, low_memory=False)
    
    mat = mat.loc[mat.index != "!series_matrix_table_end"]
    mat.dropna(axis=1, how="all", inplace=True)
    return mat

handle_ensembl = lambda names: [x.split('.')[0] for x in names]

glanced = lambda d: pd.DataFrame(d, index=[0]).T.head()

def make_naming_d(dataset_name, metaFile, sampCol="GEO_title", accessionCol="GEO_accession"):
    """ Match-Switch Sample Title and GEO Accession """
    metaFile = metaFile[metaFile["Dataset"]==dataset_name].copy()
    return dict(zip(metaFile[sampCol], metaFile[accessionCol]))

def make_feature_d(annotFile, probeCol="ID_REF", featureCol="Gene Symbol"):
    """ Match-switch platform ID (rowname) with Gene Symbol """
    annotFile = annotFile.loc[~annotFile[featureCol].isna()]
    return dict(zip(annotFile[probeCol], annotFile[featureCol]))

def preview_distribution(df, dropCol=None):
    tmp = df.copy()
    
    if dropCol is not None:
        tmp.drop(dropCol, inplace=True)
    
    plt.figure(figsize=(5,5))
    plt.hist(tmp.to_numpy().flatten())
    plt.xlabel("Binned Values")
    plt.ylabel("Frequency")
    plt.show()
    
def filter_invariant(raw_df, dropFlat0=True, lowestPerc=PERC_VAR_FILTER):
    df = raw_df.copy()
    print("Input size: %d genes x %d samples" % df.shape)
    
    zeroCond = (df==0).all(axis=1)
    print("%d (%.2f percent) genes NOT expressed" % (np.sum(zeroCond),100.0*np.mean(zeroCond)))
    
    sampSds = df.std(axis=1)
    novarCond = (sampSds == 0)
    print("%d (%.2f percent) genes INVARIANT" % (np.mean(novarCond),100.0*np.mean(novarCond)))
    
    if dropFlat0:
        df1 = df.loc[~ zeroCond].copy()
        print("After filtering all zeros: %d" % df.shape[0])
    
    if lowestPerc is None:
        return df
    elif lowestPerc == 0:
        df = df.loc[~ novarCond].copy()
        print("After filtering all zero-variant: %d" % df.shape[0])
    elif lowestPerc > 0:
        num_filter = int(lowestPerc * df.shape[0])
        print("Filtering additional %.2f percent (%d) least variable genes" % (100.0*lowestPerc,num_filter))
        
        thresh = sampSds.sort_values(ascending=True)[num_filter]
        df = df.loc[sampSds >= thresh].copy()
        print("After filtering invariant: %d" % df.shape[0])
        
    return df

def standardize_intrasamp(df, useRows=True):
    """ Intra-sample (column) standard-scaling of expression matrix """
    scaler = StandardScaler()
    newX = scaler.fit_transform(df.T)
    newX = pd.DataFrame(newX).T
    newX.set_index(df.index, inplace=True)
    newX.columns = df.columns
    return newX

## Manually curated sample sheet:
samp_meta = pd.read_excel(PATH_META, sheet_name="array_samps", na_values=["","NA","NaN"])
samp_meta = samp_meta[~samp_meta["Gender"].isna()]
samp_meta["Sex"] = samp_meta["Gender"].str.capitalize()
samp_meta.drop(["Exclude","Gender"], axis=1, inplace=True)

pd.crosstab(samp_meta["isHF"], samp_meta["Sex"])
pd.crosstab(samp_meta["isHF"], samp_meta["Age"]>50)

# ## I. Prepare individual ProbeID-to-GeneSymbol maps
gpl570 = pd.read_table(DIR+"platforms/GPL570-55999.txt", low_memory=False, na_values=STRNA)
gpl570 = gpl570.loc[~ gpl570["Sequence Type"].isin(AFFY_CTRLS)] 
gpl570 = gpl570.loc[~ gpl570["Gene Symbol"].isna()]
gpl570 = gpl570.loc[~ gpl570["Gene Symbol"].str.contains(" /// ", regex=False)]

gpl5175 = pd.read_table(DIR+"platforms/GPL5175-3188.txt", low_memory=False, na_values=STRNA)
gpl5175["ID"] = gpl5175["ID"].astype(str)
gpl5175 = gpl5175.loc[gpl5175["category"]=="main"]
gpl5175 = gpl5175.loc[gpl5175["seqname"].isin(CHROMS)]
gpl5175 = gpl5175.loc[~ gpl5175["gene_assignment"].isna()]
gpl5175["gene_assign_symbol"] = [g.split(" // ")[1].strip() for g in gpl5175["gene_assignment"]]

gpl6244 = pd.read_table(DIR+"platforms/GPL6244-17930.txt", skiprows=12, low_memory=False, na_values=STRNA)
gpl6244 = gpl6244.loc[gpl6244["category"]=="main"]
gpl6244 = gpl6244.loc[gpl6244["seqname"].isin(CHROMS)]
gpl6244 = gpl6244.loc[~ gpl6244["gene_assignment"].isna()]
gpl6244["gene_assign_symbol"] = [g.split(" // ")[1].strip() for g in gpl6244["gene_assignment"].values]
gpl6244["ID"] = gpl6244["ID"].astype(str)

gpl10558 = pd.read_table(DIR+"platforms/GPL10558-50081.txt", skiprows=30, dtype=str, low_memory=False, na_values=STRNA) #Illumina
gpl10558["seqname"] = ["chr"+str(s) for s in gpl10558["Chromosome"].values]
gpl10558 = gpl10558.loc[gpl10558["seqname"].isin(CHROMS)]
gpl10558 = gpl10558.loc[~ gpl10558["Symbol"].isna()]
gpl10558 = gpl10558.loc[~gpl10558["Symbol"].str.contains("Mar|Dec",regex=True)]

gpl11532 = pd.read_table(DIR+"platforms/GPL11532-32230.txt", skiprows=12, low_memory=False, na_values=STRNA)
gpl11532 = gpl11532.loc[gpl11532["category"]=="main"]
gpl11532 = gpl11532.loc[gpl11532["seqname"].isin(CHROMS)]
gpl11532 = gpl11532.loc[~ gpl11532["gene_assignment"].isna()]
gpl11532["gene_assign_id"] = [str(s).split(" // ")[1].strip() for s in gpl11532["gene_assignment"].values]
gpl11532["ID"] = gpl11532["ID"].astype(str)

gpl17586 = pd.read_table(DIR+"platforms/GPL17586-45144.txt", skiprows=15, low_memory=False, na_values=STRNA)
gpl17586 = gpl17586.loc[gpl17586["category"] == "main"]
gpl17586 = gpl17586.loc[gpl17586["seqname"].isin(CHROMS)]
gpl17586 = gpl17586.loc[~ gpl17586["gene_assignment"].isna()]
gpl17586["GeneSymbol"] = [g.split(" // ")[1] for g in gpl17586["gene_assignment"]]
gpl17586 = gpl17586.loc[~ gpl17586["GeneSymbol"].isna()]

## Reusable gene-feature maps as dictionaries:
f570 = make_feature_d(gpl570, "ID", "Gene Symbol")
f5175 = make_feature_d(gpl5175, "ID", "gene_assign_symbol")
f6244 = make_feature_d(gpl6244, "ID", "gene_assign_symbol")
f10558 = make_feature_d(gpl10558, "ID", "Symbol")
f11532 = make_feature_d(gpl11532, "ID", "gene_assign_id")
f17586 = make_feature_d(gpl17586, "ID", "GeneSymbol")

## II. Load & Process Individual Expression Matrices:
def wrapper_intradat(geo_accession, feature_dict, sampList=samp_meta["GEO_accession"], feature_colname="ID_REF"):
    """ Mega-wrapper to load & perform intra-dataset further processing """
    ## Load data & subset to Sample Metadata:
    gse = read_txt_matrix(geo_accession)
    if sampList is not None:    
        cond_samps = gse.columns.isin(sampList)
        if np.mean(cond_samps) < 1:
            print("Including %d of %d samples in UPDATED Sample Metadata" % (np.sum(cond_samps),gse.shape[1]))
            gse = gse.loc[:, cond_samps]
        else:
            print("All %d samples included!" % gse.shape[1])
    
    ## Subset to genes available on platform annotation
    cond_avail = gse.index.isin(feature_dict.keys()) #bool
    print("Proportion of available features %.3f" % np.mean(cond_avail))
    gse = gse.loc[cond_avail].copy()
    gse.rename(index=feature_dict, inplace=True)
    
    ## Perform per-gene aggregation, if multiple probes per gene (yes mostly)
    prop_dup = np.mean(gse.index.duplicated())
    print("Proportion of probes with 1+ genes %.3f" % prop_dup)
    if prop_dup > 0:
        print("Aggregating by gene symbol...")
        gse = gse.groupby(feature_colname).mean()
    
    ## Custom processing with visualization before & after:
    preview_distribution(gse)
    gse = filter_invariant(gse)
    gse = standardize_intrasamp(gse)
    preview_distribution(gse)
    return gse

datA = wrapper_intradat("datA", f570)
datB = wrapper_intradat("datB", f570)
datC = wrapper_intradat("datC", f570)
datD = wrapper_intradat("datD", f5175)
datE = wrapper_intradat("datE", f6244)
datF = wrapper_intradat("datF", f10558)
datG = wrapper_intradat("datG", f11532)
datH = wrapper_intradat("datH", f10558)
datI = wrapper_intradat("datI", f17586)

## III. Inner-join & Export
datArray = datA.merge(datB, left_index=True, right_index=True)\
                   .merge(datC, left_index=True, right_index=True)\
                   .merge(datD, left_index=True, right_index=True)\
                   .merge(datE, left_index=True, right_index=True)\
                   .merge(datF, left_index=True, right_index=True)\
                   .merge(datG, left_index=True, right_index=True)\
                   .merge(datH, left_index=True, right_index=True)\
                   .merge(datI, left_index=True, right_index=True)

datArray = datArray.iloc[datArray.index.argsort(), datArray.columns.argsort()]
datArray

samp_meta = samp_meta.loc[samp_meta["GEO_accession"].isin(datArray.columns)]
samp_meta = samp_meta.iloc[samp_meta["GEO_accession"].argsort()]
samp_meta["GEO_title"] = samp_meta["GEO_title"].str.replace(",| ","_", regex=True)

if all(datArray.columns==samp_meta["GEO_accession"]) and np.mean(datArray.index.duplicated())==0:
    print("Checkpoints PASSED! Proceed to exporting as CSVs...")
    datArray.to_csv(DIR+"calculated_profiles/microarray_only/merged_microarrays_for_dge.csv", index=True)
    samp_meta.to_csv(DIR+"calculated_profiles/microarray_only/microarray_samples.csv", index=False)
else:
    raise ValueError("One or more checkpoints FAILED!")
