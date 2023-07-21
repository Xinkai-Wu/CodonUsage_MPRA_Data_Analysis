import pandas as pd
import sys

"""
parameters needed:
1. fasta file for R1
2. fasta file for R2
3. output table name
4. tolerance mismatch
example usage: ./CodonUsageBias_01dataSplit.py R1.fa R2.fa out 2
"""

f_r1 = sys.argv[1]
f_r2 = sys.argv[2]
outname = sys.argv[3]
mismatch = int(sys.argv[4])

def readFasta(fastafile):
    seq = {}
    for line in fastafile:
        if line.startswith('>'):
            name = line.replace('>', '').split()[0]
            seq[name] = ''
        else:
            seq[name] += line.replace('\n', '').strip()
    return seq

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


def MismatchCount(seq1, seq2):
    """
    return the matched base count for aligned two sequence
    """
    mutation = [i for i in range(len(seq1)) if seq1[i] != seq2[i]]
    return len(mutation)


def slidingMatchSearch(query, target, mismatch):
    """
    search for matched sequence using sliding window, return the target start locus
    """
    querylen = len(query)
    start = 0
    match = False
    while start + querylen < len(target):
        seqMismatch = MismatchCount(query, target[start:start + querylen])
        if seqMismatch <= mismatch:
            match = True
            break
        start += 1
    if match:
        return start
    else:
        return 'nomatch'  # indicating no mismatch


def MainInfoFeatch(constant5prime,constant3prime,read1_,read2_rc_):
    """
    fetch variants,barcode and UMI
    """
    # read1 info
    read1_constant5 = slidingMatchSearch(constant5prime,read1_,mismatch)
    read1_constant3 = slidingMatchSearch(constant3prime,read1_,mismatch)
    
    if str(read1_constant5) != 'nomatch':
        if read1_constant5 + 63 <= len(read1_): # the minmum lenth to fetch insert seq
            r1_insert_start = read1_constant5 + 27 # the length of constant5
            r1_insert = read1_[r1_insert_start:r1_insert_start+36]
        else: 
            r1_insert = 'outRange'
    else:
        r1_insert = 'no_constant5'
    
    if str(read1_constant3) != 'nomatch':
        if read1_constant3 + 26 <= len(read1_):
            r1_barcode_start = read1_constant3 + 20
            r1_barcode = read1_[r1_barcode_start:r1_barcode_start+6]
            if read1_constant3 + 36 <= len(read1_):
                r1_UMI = read1_[r1_barcode_start+6:r1_barcode_start+16]
            else:
                r1_UMI = 'outRange'
        else:
            r1_barcode = 'outRange'
            r1_UMI = 'outRange'
    else:
        r1_barcode = 'no_constant3'
        r1_UMI = 'no_constant3'
    
    # read2_rc info
    read2_constant5 = slidingMatchSearch(constant5prime,read2_rc_,mismatch)
    read2_constant3 = slidingMatchSearch(constant3prime,read2_rc_,mismatch)
    
    if str(read2_constant5) != 'nomatch':
        if read2_constant5 + 63 <= len(read2_rc_):
            r2_insert_start = read2_constant5 + 27 
            r2_insert = read2_rc_[r2_insert_start:r2_insert_start+36]
        else: 
            r2_insert = 'outRange'
    else:
        r2_insert = 'no_constant5'
    
    if str(read2_constant3) != 'nomatch':
        if read2_constant3 + 26 <= len(read2_rc_):
            r2_barcode_start = read2_constant3 + 20
            r2_barcode = read2_rc_[r2_barcode_start:r2_barcode_start+6]
            if read2_constant3 + 36 <= len(read2_rc_):
                r2_UMI = read2_rc_[r2_barcode_start+6:r2_barcode_start+16]
            else:
                r2_UMI = 'outRange'
        else:
            r2_barcode = 'outRange'
            r2_UMI = 'outRange'
    else:
        r2_barcode = 'no_constant3'
        r2_UMI = 'no_constant3'
    
    return r1_insert,r1_barcode,r1_UMI,r2_insert,r2_barcode,r2_UMI


# Open R1 file
f_r1=open(f_r1)
read1 = readFasta(f_r1)
f_r1.close()
# Open R2 file
f_r2=open(f_r2)
read2 = readFasta(f_r2)
f_r2.close()

# merge paired reads
df1 = pd.DataFrame([read1])
df1 = pd.DataFrame(df1.values.T, index=df1.columns, columns=df1.index)
df1 = df1.reset_index()
df1.columns = ['seqid','read1']
df2 = pd.DataFrame([read2])
df2 = pd.DataFrame(df2.values.T, index=df2.columns, columns=df2.index)
df2 = df2.reset_index()
df2.columns = ['seqid','read2']
read_df = pd.merge(df1,df2,on='seqid')

# to identify whether the constant sequence in paired reads
constant5prime = 'CACAACGTCTATATCATGGCCGACAAG'
constant3prime = 'TGCAGCTCGCCGACCACTAC' # reverse compliment to match

read_df['read2_rc'] = read_df.apply(lambda x: baseRevCom(x['read2']),axis = 1)
del read_df['read2']

fetchFuc = lambda x: pd.Series(MainInfoFeatch(constant5prime,constant3prime,x.read1,x.read2_rc))
newcols = read_df.apply(fetchFuc, axis=1)
newcols.columns = ['r1_insert','r1_barcode','r1_UMI','r2_insert','r2_barcode','r2_UMI']
read_df = read_df.join(newcols)

read_df.to_csv(outname + 'all.out' , sep="\t", index=False)
del read_df['read1']
del read_df['read2_rc']
read_df.to_csv(outname + 'clean.out', sep="\t", index=False)
