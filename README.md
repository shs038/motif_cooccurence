# motif_cooccurence
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
import scipy
from scipy import stats
%matplotlib inline
def co_occur(n):
    '''
    input is a matix contains all data of the genomic position of motif in chr1.
    output is a matix contains the frequency of co-occurance of every two motifs. 
    '''
    (row,col)=n.shape #find the dimension of the data
    count_matrix=np.zeros((col,col),dtype=np.int) #a matrix to store future data
    col_vector=np.zeros((row,1),dtype=np.int) #a column vector use in loop
    for i in range (col):#use the column vector to store each column one by one
        col_vector=n[:,i]!= -1# Find if the motifs exist
        for j in range (col): 
            log_v=n[:,j]!= -1#Find if the motifs exist
            vector_sum = 1*col_vector + 1*log_v
            log_input = vector_sum == 2
            input=np.sum(1*log_input)#convert logical ou
            count_matrix[i,j]=count_matrix[i,j]+input
            np.fill_diagonal(count_matrix, 0)# change diagonal back to zero
    return count_matrix
    # read in data
df = pd.read_excel('/Users/sun/Desktop/Lab/chr1new.xlsx',parse_cols="B:GO")#load part of file contain data
data_list=df.as_matrix()#conver data into matrix form
# find co-occurence
co=co_occur(data_list)
motifs = df.columns.values
co_frame = pd.DataFrame(co, columns = motifs)
co_frame.index = motifs
#plot z socre for each row of motif(one motif co-occur with all the other)
def motif_co(motif_id):
    #plot the z score for each motif, i in range z_score.shape[0]
    sns.distplot(co_frame[motif_id])
    plt.ylabel('Frequency')
    plt.xlabel('co_occur')
    plt.title('co_occurance') 
motif_co('ap-1')
def Find_Z(n):
    '''
    For all pairs of motifs - is there a pair that co-occurs more or less often than you would expect.
    input: a matirc contains frequency of co-occurence of every pair of motifs.
    output:an array contains z score of each pair of motifs.
    '''
    z_matrix = np.zeros((n.shape[0],n.shape[1]-1),dtype=np.float)
    for i in range (n.shape[0]):
        co_motif = co[i,:]
        co_motif = np.delete(co_motif,i)#remove data of the motif co-occur with itself
        z_score=stats.zscore(co_motif)# find z socre 
        z_matrix[i,:]=z_score
    return z_matrix
# find z score of co-occurence data for each row of motif
z_score=Find_Z(co)
#assign co-occurance z socre of one motif with itself as 100 to make dataframe 
zscore_all=np.zeros((z_score.shape[0],z_score.shape[0]),dtype=np.float)
for i in range (z_score.shape[0]):
        z_motif_self=z_score[i,:]
        z_motif_self=np.insert(z_motif_self,i,100)
        zscore_all[i,:]=z_motif_self
print zscore_all
zscore_frame = pd.DataFrame(zscore_all, columns=motifs)
zscore_frame.index = motifs
def motif_z(i):
    #plot the z score for each motif, i in range z_score.shape[0]
    sns.distplot(z_score[i,:])
    plt.ylabel('Frequency')
    plt.xlabel('z')
    plt.title('z_score') 
    print motifs[i]
motif_z(0)
# find significant pair of motif with certain z score threshold
def Find_sigpair(siglevel):
    threshold = siglevel # threshold of z score
    significantPairCount = 0
    for i in range(co.shape[0]):
        for j in range(co.shape[0]):
            current_zscore = zscore_all[i][j]
            if current_zscore >= threshold and i != j:#ignore the data of one motif with itself
                # display motif pair
                motif1 = motifs[i]
                motif2 = motifs[j]
                print(motif1, motif2)
                significantPairCount += 1
    return(significantPairCount)
sigpaircount=Find_sigpair(2.25)
least_fq_z=zscore_frame.min(axis=1)
least_fq_id=zscore_frame.idxmin(axis=1)
for i in range(196):
    zscore_frame.ix[i,i]=-100
most_fq_z=zscore_frame.max(axis=1)
most_fq_id=zscore_frame.idxmax(axis=1)
df1=most_fq_id.to_frame(name='most_fq_id')
df2=most_fq_z.to_frame(name='most_fq_z')
df3=least_fq_id.to_frame(name='least_fq_id')
df4=least_fq_z.to_frame(name='least_fq_z')
result= pd.merge(pd.merge(pd.merge(df1,df2,left_index=True, right_index=True),df3,left_index=True, right_index=True),df4,left_index=True, right_index=True)
result['new col'] = range(result.shape[0])
del result['new col'] 
result.to_csv('/Users/sun/Desktop/results.tsv', sep='\t')
from collections import Counter
most_frequent_motifs = result['most_fq_id'].values # retreive motifs that co-occur the most often
most_frequent_motifs = list(most_frequent_motifs) # convert numpy array to Python list
c = Counter(most_frequent_motifs) # create a counter to count ocurrences of each motif
c = dict(c) # convert counter to a dictionary
motif_names = []
counts = []
for key in c:
    motif_names.append(key)
    counts.append(c[key])
frame = pd.DataFrame({'Motif Name':motif_names,
                       'Count': counts}, 
                      )

sns.factorplot(data=frame.sort('Count'), 
               x='Motif Name', 
               y='Count',
               size=8,
               kind='bar')
plt.xticks(rotation=90)
most_frequent_motifs = result['least_fq_id'].values # retreive motifs that co-occur the most often
most_frequent_motifs = list(most_frequent_motifs) # convert numpy array to Python list
c = Counter(most_frequent_motifs) # create a counter to count ocurrences of each motif
c = dict(c) # convert counter to a dictionary
motif_names = []
counts = []
for key in c:
    motif_names.append(key)
    counts.append(c[key])
frame = pd.DataFrame({'Motif Name':motif_names,
                       'Count': counts}, 
                      )

sns.factorplot(data=frame.sort('Count'), 
               x='Motif Name', 
               y='Count',
               size=8,
               kind='bar')
plt.xticks(rotation=90)
plt.title('Most common least co-ocurring motifs')
