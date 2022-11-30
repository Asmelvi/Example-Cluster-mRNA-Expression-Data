###### PART OF IMPORTS
import types

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pylab as plt
from scipy.cluster.hierarchy import dendrogram
import seaborn as sns
import scipy
import sys

###### PART OF GENERAL CONFIGURATIONS

np.set_printoptions(threshold=np.inf, precision=2)

infile = open("C:/Users/Avalon/Desktop/Cluster mRNA/lymphoma.res", 'r')

#################### PART 1

counter = 0
list_of_values_per_row = []
list_of_lists = []

for line in infile:
    if counter <=2:
        counter+=1
    elif counter >2:
        words = line.split('\t')
        for word in words:
            try:
                value = float(word)
                list_of_values_per_row.append(float(value))
            except:
                pass
        list_of_lists.append(list_of_values_per_row)
        list_of_values_per_row = []

lymp_numpy_array = []
lymp_numpy_array2 = []

for list in list_of_lists:
    lymp_numpy_array.append(np.array(list))

lymp_numpy_array2 = np.array(lymp_numpy_array)

#print(lymp_numpy_array2) ##PART 1 OUTPUT


#################### PART 2

Z = scipy.cluster.hierarchy.average(lymp_numpy_array2)
#print(Z) ##PART 2 OUTPUT

#################### PART 3

leaves = scipy.cluster.hierarchy.leaves_list(Z)
ordered_row_lymp = lymp_numpy_array2[leaves]
#print(ordered_row_lymp) ##PART 3 OUTPUT

#plt.matshow(lymp_numpy_array2) ###GRAPH 1
#plt.show() ###GRAPH 1

#plt.matshow(ordered_row_lymp) ###GRAPH 2
#plt.show() ###GRAPH 2

#################### PART 4

column_lymp = ordered_row_lymp.transpose()
Z_c = scipy.cluster.hierarchy.average(column_lymp)

#ordering the column and transposing to get the complete row and column ordered matrix
leaves_c = scipy.cluster.hierarchy.leaves_list(Z_c)
ordered_column_lymp = column_lymp[leaves_c]
ordered_lymp = ordered_column_lymp.transpose()
#print(ordered_lymp) #PART 4 OUTPUT

##Now we print the clusters genes/strains

#Genes

fcluster_row = scipy.cluster.hierarchy.fcluster(Z,t=20, criterion='maxclust')
count=0
lst=[]
description = ""
for i in fcluster_row:
    if i==1:
        lst.append(count)
    count+=1
counter_row=-3

#WE HAVE TO OPEN THE FILE AGAIN SINCE WE ALREADY REACHED THE BOTTON OF THE FILE AND IT WILL NOT RE-READ IT AGAIN
infile = open("C:/Users/Avalon/Desktop/Cluster mRNA/lymphoma.res",'r')

list_of_descriptions = []
for line2 in infile:
    if counter_row in lst:
        description = line2.split('\t')[0]
        list_of_descriptions.append(description)
    counter_row+=1

list_of_descriptions.sort()
#print(list_of_descriptions) #PART 4 OUTPUT

#Strains

fcluster_col = scipy.cluster.hierarchy.fcluster(Z_c,t=20, criterion='maxclust')
count2 = -2
lst2 = []
strain = ""
for i in fcluster_col:
    if i==1:
        lst2.append(count2)
    count2+=1

#WE HAVE TO OPEN THE FILE AGAIN SINCE WE ALREADY REACHED THE BOTTON OF THE FILE AND IT WILL NOT RE-READ IT AGAIN
infile = open("C:/Users/Avalon/Desktop/Cluster mRNA/lymphoma.res", 'r')
counter_col = 0
list_of_strains = []
strain_names = infile.readlines()[0].split('\t')

for word in strain_names:
    if counter_col in lst2:
        list_of_strains.append(word)
    counter_col+=1

list_of_strains.sort()
#print(list_of_strains) #PART 4 OUTPUT

#################### PART 5

#Genes

dict_row = {}

for value in fcluster_row:
    if value not in dict_row.keys():
        dict_row[value]=1
    else:
        dict_row[value]+=1

sum = 0
for key in dict_row.keys():
    if dict_row[key] > sum:
        sum = dict_row[key]
        id = key

count=0
lst=[]
for i in fcluster_row:
    if i==id:
        lst.append(count)
    count+=1

#WE HAVE TO OPEN THE FILE AGAIN SINCE WE ALREADY REACHED THE BOTTON OF THE FILE AND IT WILL NOT RE-READ IT AGAIN
infile = open("C:/Users/Avalon/Desktop/Cluster mRNA/lymphoma.res",'r')

list_of_descriptions = []
counter_row=-3
for line in infile:
    if counter_row in lst:
        description = line.split('\t')[0]
        list_of_descriptions.append(description)
    counter_row += 1
list_of_descriptions.sort()
#print(str(sum) + str(list_of_descriptions)) #PART 5 OUTPUT

#Strains

dict_col = {}

for value in fcluster_col:
    if value not in dict_col.keys():
        dict_col[value]=1
    else:
        dict_col[value]+=1

sum2 = 0
for key in dict_col.keys():
    if dict_col[key] > sum:
        sum2 = dict_col[key]
        id2 = key

count=0
lst2=[]
for i in fcluster_col:
    if i==id2:
        lst2.append(count)
    count+=1

#WE HAVE TO OPEN THE FILE AGAIN SINCE WE ALREADY REACHED THE BOTTON OF THE FILE AND IT WILL NOT RE-READ IT AGAIN
infile = open("C:/Users/Avalon/Desktop/Cluster mRNA/lymphoma.res", 'r')
counter_col = -2
list_of_strains = []
strain_names = infile.readlines()[0].split('\t')

for word in strain_names:
    if counter_col in lst2:
        list_of_strains.append(word)
    counter_col+=1

list_of_strains.sort()
#print(str(sum2) + str(list_of_strains)) #PART 5 OUTPUT

###################################### PLOTTING #########################################################

#WE HAVE TO OPEN THE FILE AGAIN SINCE WE ALREADY REACHED THE BOTTON OF THE FILE AND IT WILL NOT RE-READ IT AGAIN

fig = plt.figure(figsize=(8,8))
ax1 = fig.add_axes([0.09,0.1,0.2,0.6])
z1 = scipy.cluster.hierarchy.dendrogram(Z,orientation='left')
ax1.set_xticks([])
ax1.set_yticks([])

# Compute and plot second dendrogram.
ax2 = fig.add_axes([0.3,0.71,0.6,0.15])
z2 = scipy.cluster.hierarchy.dendrogram(Z_c)
ax2.set_xticks([])
ax2.set_yticks([])

# Plot distance matrix.
axmatrix = fig.add_axes([0.3,0.1,0.6,0.6])
im = axmatrix.matshow(ordered_lymp, aspect='auto', origin='lower')
axmatrix.set_xticks([])
axmatrix.xaxis.set_tick_params(pad=95)

axmatrix.yaxis.tick_right()

fig.subplots_adjust(bottom=0.1)
axcolor = fig.add_axes([0.3,0.07,0.6,0.02])
fig.colorbar(im,cax=axcolor,orientation='horizontal')
#plt.show()
#plt.savefig(sys.argv[3],dpi=600)  #LAST GRAPH