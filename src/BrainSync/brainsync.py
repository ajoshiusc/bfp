import scipy as sp
from numpy.random import random
from scipy.stats import special_ortho_group
from tqdm import tqdm
"""
Created on Tue Jul 11 22:42:56 2017
Author Anand A Joshi (ajoshi@usc.edu)
"""


def normalizeData(pre_signal):
    """
     normed_signal, mean_vector, std_vector = normalizeData(pre_signal)
     This function normalizes the input signal to have 0 mean and unit
     norm in time.
     pre_signal: Time x Original Vertices data
     normed_signal: Normalized (Time x Vertices) signal
     mean_vector: 1 x Vertices mean for each time series
     norm_vector : 1 x Vertices norm for each time series
    """

#    if sp.any(sp.isnan(pre_signal)):
#        print('there are NaNs in the data matrix, making them zero')

    pre_signal[sp.isnan(pre_signal)] = 0
    mean_vector = sp.mean(pre_signal, axis=0, keepdims=True)
    normed_signal = pre_signal - mean_vector
    norm_vector = sp.linalg.norm(normed_signal, axis=0, keepdims=True)
    norm_vector[norm_vector == 0] = 1e-116
    normed_signal = normed_signal / norm_vector

    return normed_signal, mean_vector, norm_vector


def brainSync(X, Y):
    """
   Input:
       X - Time series of the reference data (Time x Vertex) \n
       Y - Time series of the subject data (Time x Vertex)

   Output:
       Y2 - Synced subject data (Time x Vertex)\n
       R - The orthogonal rotation matrix (Time x Time)

   Please cite the following publication:
       AA Joshi, M Chong, RM Leahy, BrainSync: An Orthogonal Transformation
       for Synchronization of fMRI Data Across Subjects, Proc. MICCAI 2017,
       in press.
       """
    if X.shape[0] > X.shape[1]:
        print('The input is possibly transposed. Please check to make sure \
that the input is time x vertices!')

    C = sp.dot(X, Y.T)
    U, _, V = sp.linalg.svd(C)
    R = sp.dot(U, V)
    Y2 = sp.dot(R, Y)
    return Y2, R


def groupBrainSync(S):
    # Group BrainSync algorithm developed by Haleh Akrami

    numT = S.shape[0]
    numV = S.shape[1]
    SubNum = S.shape[2]

    # init random matrix for Os
    Os = sp.zeros((numT, numT, SubNum))
    for i in range(SubNum):  #initializeing O
        #        R = 2 * rnd.random(size=(numT, numT)) - 1; #define a random matrix with unity distributian from -1 to 1
        Os[:, :, i] = special_ortho_group.rvs(
            numT)  #(sp.dot(R , R.T)^(-1/2) , R;  #orthogonal rows of matrix

    Error = 1
    PreError = 1
    relcost = 1

    alpha = 1e-6
    var = 0
    Costdif = sp.zeros(10000)

    print('init done')

    # Initialize PreError from gloal average
    X = sp.zeros((numT, numV))
    for j in range(SubNum):  #calculate X
        X = sp.dot(Os[:, :, j], S[:, :, j]) + X

    X = X / SubNum
    InError = 0

    for j in range(SubNum):
        etemp = sp.dot(Os[:, :, j], S[:, :, j]) - X
        InError = InError + sp.trace(sp.dot(etemp,
                                            etemp.T))  #calculating error

    # Find best Orthogognal map, by minimizing error (distance) from average
    while relcost > alpha:
        var = var + 1

        print('subject iteration')
        for i in tqdm(range(SubNum)):
            X = sp.zeros((numT, numV))
            for j in range(SubNum):  #calculate X average excluded subject i
                if j != i:
                    X = sp.dot(Os[:, :, j], S[:, :, j]) + X
            # Y is i excluded average
            Y = X / (SubNum - 1)

            # Update Orthogonal matrix with BrainSync projection technique
            U, _, V = sp.linalg.svd(sp.dot(Y, S[:, :, i].T))
            Os[:, :, i] = sp.dot(U, V.T)

        print('calculate error')
        Error = 0
        # New Average with all subject updated orthogonal matrix
        # update last subject outside loop
        X2 = (X + sp.dot(Os[:, :, i], S[:, :, i])) / SubNum

        # Calculate error of all subjects from average map
        for j in range(SubNum):
            etemp = sp.dot(Os[:, :, j], S[:, :, j]) - X2
            Error = Error + sp.trace(sp.dot(etemp,
                                            etemp.T))  #calculating error

        relcost = sp.abs(Error - PreError) / sp.abs(InError)
        Costdif[var] = PreError - Error
        PreError = Error

        var
        relcost

    Costdif[var:] = []
    Costdif = Costdif[1:]
    TotalError = Error

    return X2, Os, Costdif, TotalError