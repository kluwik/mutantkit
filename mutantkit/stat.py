#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
import os

import matplotlib.pyplot as plt
import seaborn as sns

import math
from sklearn.metrics import mean_squared_error
from scipy.stats import mannwhitneyu, pearsonr, spearmanr

################################# DDG DATAFRAME PROCESSING AND STATISTICS #################################

def splitMutation(dataframe, mutation_column):
    df = dataframe.copy()
    df['wt_aa'] = df[mutation_column].apply(lambda x: x[0])
    df['position'] = df[mutation_column].apply(lambda x: x[1:-1]).astype('int64')
    df['mut_aa'] = df[mutation_column].apply(lambda x: x[-1])
    return df

def removeOutliers(df, y_ddg):
    Q1 = df[y_ddg].quantile(0.25)
    Q3 = df[y_ddg].quantile(0.75)
    IQR = Q3-Q1
    df_trim = df[~((df[y_ddg]<(Q1-1.5*IQR)) | (df[y_ddg]>(Q3+1.5*IQR)))]
    return df_trim

def getDDG(data, y_ddg, outliers=False):
    df = data[data[y_ddg].notna()] # remove NaN ddg
    if outliers == True:
        df = removeOutliers(df, y_ddg)
    return df

def getStat(x,y, prnt=False):
    pcc = pearsonr(x,y)[0]
    sp = spearmanr(x,y)[0]
    p_pcc = pearsonr(x,y)[1]
    se = math.sqrt((1 - pcc**2)/(len(x)-2))
    sp_se = math.sqrt((1 - sp**2)/(len(x)-2))
    mse = mean_squared_error(x,y, squared=False)
    mwu = mannwhitneyu(x,y)[1]
    if prnt == True:
        print(f'PCC = {pcc:.2f}±{se:.2f}, SCC = {sp:.2f}±{sp_se:.2f}, MSE = {mse:.2f}')
    return round(pcc,2), round(se,2), round(sp,2), round(sp_se,2), round(mse,2)
             
ddg_class_2 = lambda x: 'D' if x >= 0 else 'S'    
ddg_class_3 = lambda x: 'D' if x >= 0.5 else ('S' if x <= -0.5 else 'N') 
    
def getAntisymmetry(f, r):
    pcc =  getStat(f, r)[0]
    b = np.sum((f + r)) / (2*len(f))
    return round(pcc,2), round(b,2)

