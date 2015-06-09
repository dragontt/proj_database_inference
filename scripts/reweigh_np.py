import sys
import os
import numpy as nmp
import reweigh_pwm_matrix as rw_pwm

""" 
A function mostly containing importables for reweighing netprophet edges

This module is closely related to the reweigh_pwm_matrix module, which indicates
a need for refactoring TODO
"""


def reweigh_by_tf_dot_corr(M, W):
    """ Reweigh each row in M based on how well it matches with its 
    corresponding row in W.  This function uses the dot product and the
    pearson's correlation as similarity measures. """
    weights = rw_pwm.calculate_corr_dot_pval_weights(W, M)
    retval = nmp.ones(nmp.shape(M))
    for ix in xrange(nmp.shape(M)[0]):
        retval[ix, :] *= weights[ix]
    return retval
#end function

def reweigh_rows(M, row_weights):
    """ Given a matrix M and a list of weights, return a copy of M with each
    row multiplied by its corresponding weight. """
    retval = nmp.copy(M)
    for n, val in enumerate(row_weights):
        retval[n, :] *= val
    return retval
#end function

def reweigh_matrix(M, W):
    """ Given a target matrix M and a weight matrix W, 
    return a copy of M with each item multiplied by the corresponding
    item in W. """
    retval = nmp.copy(M)
    retval *= W
    return retval # there, that was easy!
#end function

def main():
    pass
#end function

if __name__ == "__main__":
    main(sys.argv)
