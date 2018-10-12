""" Function to simulate overlapping pools of constructs in the wells
and count number of combinations covered in the entire screen


Variables:
v1 - number of genes
v2 - number of wells
v3 - proportion of a library as number of genes/well
v4 - number of genes in pathway
v5 - number of iterations in bootstrap


Constants:
c1 - standard deviation (0.15 * v2) for normal distr of guides/well
c2 - well exclusion (0.025 * v1)

"""

import pandas as pd
import numpy as np
import itertools as it
import math

# Help message
import argparse
from argparse import RawTextHelpFormatter

#from tqdm import tqdm

# calculate number of combinations
def nCr(n,r):
    f = math.factorial
    return f(n) // f(r) // f(n-r)


def combCovEstim (v1, v2, v3, v4, c1=0.15, c2=0.025):

    # make sensing matrix

    A = pd.DataFrame(np.zeros(shape=(v2, v1)))

    for index_pos, row in A.iterrows():

        # model number of sgRNA/well
        sizeSamp = round(np.random.normal(loc=round(v3 * v1), scale=round(v1 * v3 * c1)))

        colsFill = np.random.choice(a=v1, size=abs(sizeSamp), replace=True)

        # fill with 1
        A.iloc[index_pos, colsFill] = 1

    # exclude some "wells" of the array
    wellNumbEx = np.random.choice(a=v2, size=round(v2 * c2).astype(int), replace=False)
    A = A.drop(index=wellNumbEx)


    # mark all possible combinations in each well
    res = pd.DataFrame()

    for index_comb, row in A.iterrows():

        # define gene "IDs" with in the well
        geneIdWell = A.columns[A.loc[index_comb, ] == 1]

        combDF = pd.DataFrame(list(it.combinations(geneIdWell, v4)))

        res = res.append(combDF, ignore_index=True)

        # print index
        print(index_comb)

    # sort every row in place
    res.values.sort(axis=1)

    # exclude duplicated combinations
    res.drop_duplicates(inplace=True)

    # return number of unique combinations

    return res.shape[0]


### Simulation for array of variables

if __name__ == "__main__":

    v1_vect = np.linspace(50, 250, 5).astype(int)
    v2_vect = np.linspace(48, 192, 4).astype(int)
    v3_vect = np.linspace(0.25, 0.75, 3)
    v4_vect = np.linspace(3, 7, 5).astype(int)
    v1=250
    v5 = 10

    for v4 in v4_vect:

        # create output dataFrame
        out_df = pd.DataFrame(np.zeros(shape=(len(v3_vect), len(v2_vect))))

        # total number
        total_Nr = nCr(n=v1, r=v4)

        for v3 in v3_vect:
            for v2 in v2_vect:

                # list for results from iterations
                list_for_mean = [0] * v5

                # repeat v5 times
                for i in range(v5):
                    coverage = combCovEstim(v1=250, v2=v2, v3=v3, v4=v4, c1=0.15, c2=0.025)

                    # export proportion
                    list_for_mean[i] = coverage / total_Nr

                    #

                out_df.iloc[np.where(v3_vect == v3)[0], np.where(v2_vect == v2)[0]] = np.mean(list_for_mean[i])

        out_df.to_csv(path_or_buf="/home/anton/02_PhenoSudoku/oak180823_sudoku0014_PathwaySudokuModel/results"
                                  + str(v4) + ".txt", sep='\t', header=False, index=False)






