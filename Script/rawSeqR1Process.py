import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt
import re
import sys

sample = sys.argv[1]
regionX = sys.argv[2]


def synonymousConstraint(variant,regionX):
    """
    to identify whether the variant is synonymous
    """
    L = ['CTA','CTT','CTC','CTG']
    T = ['ACA','ACT','ACC','ACG']
    F = ['TTC','TTT']
    I = ['ATA','ATT','ATC']
    G = ['GGA','GGT','GGC','GGG']
    K = ['AAA','AAG']
    C = ['TGC','TGT']
    Q = ['CAA','CAG']
    N = ['AAC','AAT']
    V = ['GTA','GTT','GTC','GTG']
    R = ['CGA','CGT','CGC','CGG']
    
    codon1,codon2,codon3,codon4,codon5,codon6,codon7,codon8,codon9,codon10,codon11,codon12 = re.findall('.{'+str(3)+'}', variant)
    
    if str(regionX) == 'region2': 
        keep = True
        if (codon1 not in L)|(codon2 not in T)|(codon3 not in L)|(codon4 not in K)|(codon5 not in F)|(codon6 not in I):
            keep = False
        if (codon7 not in C)|(codon8 not in T)|(codon9 not in T)|(codon10 not in G)|(codon11 not in K)|(codon12 not in L):
            keep = False
        
    if str(regionX) == 'region5': 
        keep = True
        if (codon1 not in Q)|(codon2 not in K)|(codon3 not in N)|(codon4 not in G)|(codon5 not in I)|(codon6 not in K):
            keep = False
        if (codon7 not in V)|(codon8 not in N)|(codon9 not in F)|(codon10 not in K)|(codon11 not in I)|(codon12 not in R):
            keep = False
            
    return keep

# filter reads with wrong barcode
def baseRevCom(seq):
    """
    reverse and compliment
    """
    baseRev = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    seqcom = []
    # com
    for i in range(len(seq)):
        com = baseRev[seq[i]]
        seqcom.append(com)
    # rev
    seqcom.reverse()
    revcom = ''.join(seqcom)
    return (revcom)

def barcodeSearch(barcode):
    """input: barcode sequence, return fraction name"""
    
    validBarcode = {'GTAGAG':'bin1','GTCCGC':'bin2','GTGAAA':'bin3','GTGGCC':'bin4','GTTTCG':'bin5','CGTACG':'bin6','GAGTGG':'bin7','GGTAGC':'bin8'}
    barcode = baseRevCom(barcode)
    if barcode in validBarcode:
        fraction = validBarcode[barcode]
        return fraction
    else:
        return "WrongBarcode"

def cleanOutProcessRead1(sample,regionX):
    print("Your sample: "+sample + " Analysis start ...")
    sampleFile = sample + '.clean.out.gz'
    df = pd.read_csv(sampleFile,sep='\t',compression='gzip')
    df_read1 = df[['seqid','r1_insert','r1_barcode','r1_UMI']]
    print(sample + " total reads: "+ str(df_read1.shape[0]))
    
    # filter reads could not find constant region
    df_read1_lociRight = df_read1[(~df_read1['r1_insert'].str.contains('o'))&(~df_read1['r1_barcode'].str.contains('o'))&(~df_read1['r1_UMI'].str.contains('o'))]
    print(sample + "Constant region found reads: "+ str(df_read1_lociRight.shape[0]))
    # barcode
    df_read1_lociRight['fraction'] = df_read1_lociRight.apply(lambda x: barcodeSearch(x['r1_barcode']),axis = 1)
    df_read1_lociRight_barcode = df_read1_lociRight[df_read1_lociRight['fraction'] != 'WrongBarcode']
    print(sample + "Right barcode: "+ str(df_read1_lociRight_barcode.shape[0]))
    
    # 保留同义突变序列
    df_read1_lociRight_barcode['syn'] = df_read1_lociRight_barcode.apply(lambda x: synonymousConstraint(x['r1_insert'],regionX),axis = 1)
    df_read1_lociRight_barcode_syn = df_read1_lociRight_barcode[df_read1_lociRight_barcode['syn'] == True]
    print(sample + "Synonymous mutation: "+ str(df_read1_lociRight_barcode_syn.shape[0]))
    
    del df_read1_lociRight_barcode_syn['syn']
    # 去除UMI重复的reads
    df_read1_lociRight_barcode_syn_dropDup = df_read1_lociRight_barcode_syn[['r1_insert','fraction','r1_UMI']].drop_duplicates()
    # 将不同的UMI进行累加
    df_read1_lociRight_barcode_syn_dropDup_count = df_read1_lociRight_barcode_syn_dropDup.groupby(['r1_insert','fraction'])['r1_UMI'].size().reset_index(name='count')
    
    df_read1_lociRight_barcode_syn_dropDup_count.columns = ['variable','barcode','Count']
    df = df_read1_lociRight_barcode_syn_dropDup_count
    
    for bin in np.unique(df[['barcode']].values):
        tmp = df[df.barcode == bin][['variable','Count']]
        tmp.columns = ['variable',bin]
        if bin == 'bin1':
            df_tmp = tmp
        else:
            df_tmp = pd.merge(df_tmp,tmp,on='variable',how='outer')
            
    outFile = sample + '.Read1.bin.gz'
    df_tmp.to_csv(outFile, sep="\t", index=False,compression='gzip')

cleanOutProcessRead1(sample,regionX)

