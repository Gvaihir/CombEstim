""" Function to simulate overlapping pools of constructs in the wells
and count number of combinations covered in the entire screen


Variables:
v1 - number of genes
v2 - number of wells
v3 - proportion of a library as number of genes/well
v4 - number of genes in pathway

Constants:

c1 - standard deviation (0.15 * v2) for normal distr of guides/well
c2 - well exclusion (0.025 * v1)

"""
v1=150
v2=100
v3=0.4
v4=6



import pandas as pd
import numpy as np
import itertools as it


def combCovEstim (v1, v2, v3, v4, c1 = 0.15, c2 = 0.025):

    # make sensing matrix

    A = pd.DataFrame(np.zeros(shape=(v2, v1)))

    for index_pos, row in A.iterrows():

        # model number of sgRNA/well
        sizeSamp = round(np.random.normal(loc = round(v3 * v1), scale = round(v1 * v3 * c1)))

        colsFill = np.random.choice(a = v1, size = abs(sizeSamp), replace = True)

        # fill with 1
        A.iloc[index_pos, colsFill] = 1


    # mark all possible combinations in each well
    res = pd.DataFrame()

    for index_comb, row in A.iterrows():

        # define gene "IDs" with in the well
        geneIdWell = A.columns[A.iloc[index_comb,] == 1]

        combDF = pd.DataFrame(list(it.combinations(geneIdWell,v4)))

        res = res.append(combDF, ignore_index = True)


    # sort every row in place
    res.values.sort(axis=1)

    # exclude duplicated combinations
    res.drop_duplicates(inplace=True)

    # return number of unique combinations

    return res.shape[0]

