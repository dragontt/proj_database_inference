import sys
import numpy as nmp
import scipy
from scipy import stats

"""
"""

def safelog10(x):
    return nmp.log10(nmp.array(x)+.0001)
#end function

def just_corr(x, y):
    # first element is corr, which is contained in an array
    return stats.pearsonr(x, y)[0]
#end function

def covariance(v1, v2):
    '''  '''
    centered_copy_v1 = nmp.copy(v1) - nmp.mean(v1)
    centered_copy_v2 = nmp.copy(v2) - nmp.mean(v2)
    return nmp.dot(centered_copy_v1, centered_copy_v2)
#end function

def cosine_similarity(v1, v2):
    ''' Given two vectors v1 and v2, return the cosine of the angle between them. '''
    import numpy as nmp
    import operator as op
    numer = nmp.dot(v1, v2)
    def l2_norm(vec):
        return nmp.sqrt(nmp.sum(nmp.power(vec,2)))
    denom = reduce(op.mul, (map(l2_norm, [v1, v2])))
    if denom == 0:
        return 0
    return numer/denom
#end function

def reweigh_np_both(pwm_adjmtr, np_adjmtr):
    """ Given a pwm and a np adjmtr, reweighs np by using its correlation
    with the pwm. """
    pass
#end function

def reweigh_pwm_all(pwm_adjmtr, np_adjmtr):
    weights = calculate_

def reweigh_pwm_both(pwm_adjmtr, np_adjmtr):
    return reweigh_pwm_adjmtr(pwm_adjmtr,
                              calculate_corr_dot_pval_weights(pwm_adjmtr,
                                                              np_adjmtr))
#end function

def reweigh_pwm_dot(pwm_adjmtr, np_adjmtr):
    """ returns a copy of pwm_adjmtr where each row is reweighed
    by the -log10 of its empirical pval of similarity to its corresponding
    np_adjmtr row. """
    return reweigh_pwm_adjmtr(pwm_adjmtr, 
                              calculate_dotprod_pval_weights(pwm_adjmtr, 
                                                             np_adjmtr))
#end function

def reweigh_pwm_corr(pwm_adjmtr, np_adjmtr):
    return reweigh_pwm_adjmtr(pwm_adjmtr,
                              calculate_corr_pval_weights(pwm_adjmtr,
                                                          np_adjmtr))
#end function

def reweigh_pwm_cov(pwm_adjmtr, np_adjmtr):
    return reweigh_pwm_adjmtr(pwm_adjmtr,
                              calculate_metric_pval(pwm_adjmtr,
                                                    np_adjmtr,
                                                    covariance))
#end function

def calculate_all_pval_weights(pwm_adjmtr, np_adjmtr, metrics):
    pvals = map(lambda x: calculate_metric_pval(pwm_adjmtr, np_adjmtr, x),
                metrics)
    avged = nmp.mean(nmp.array(pvals), 0)
    return -1*safelog10(avged)
#end function

def calculate_corr_dot_pval_weights(pwm_adjmtr, np_adjmtr):
    return calculate_all_pval_weights(pwm_adjmtr, np_adjmtr, [just_corr, nmp.dot])
    # corr_pvals = calculate_metric_pval(pwm_adjmtr, np_adjmtr, just_corr)
    # dot_pvals = calculate_metric_pval(pwm_adjmtr, np_adjmtr, nmp.dot)
    # total_pvals = nmp.mean(nmp.vstack([corr_pvals, dot_pvals]), 0)
    # return -1*safelog10(total_pvals)
#end function

def calculate_corr_pval_weights(pwm_adjmtr, np_adjmtr):
    return -1*safelog10(calculate_metric_pval(pwm_adjmtr, np_adjmtr, just_corr))
#end function

def calculate_dotprod_pval_weights(pwm_adjmtr, np_adjmtr):
    return -1*safelog10(calculate_metric_pval(pwm_adjmtr, np_adjmtr, nmp.dot))
#end function

def calculate_metric_pval(pwm_adjmtr, np_adjmtr, metric_fn):
    pvals = []
    num_tfs = nmp.shape(np_adjmtr)[0]
    # loop through each row of the pwm_adjmtr
    for jx in xrange(nmp.shape(pwm_adjmtr)[0]):
        num_less_than = 0;
        my_metric = metric_fn(pwm_adjmtr[jx, :], np_adjmtr[jx, :])
        # loop through each row of the np_adjmtr
        for kx in xrange(nmp.shape(np_adjmtr)[0]):
            if kx == jx:
                continue
            other_metric = metric_fn(pwm_adjmtr[jx, :], np_adjmtr[kx, :])
            if my_metric < other_metric:
                num_less_than += 1
        pvals.append(float(num_less_than)/num_tfs)
    return pvals
#end function

def reweigh_pwm_adjmtr(pwm_adjmtr, weights):
    """ return a copy of pwm_adjmtr where each row has been multiplied
    by the corresponding value in weights. """
    if len(weights) != nmp.shape(pwm_adjmtr)[0]:
        raise ValueError("must have as many weights as rows in pwm_adjmtr!")
    reweighed = nmp.copy(pwm_adjmtr)
    for n, i in enumerate(weights):
        reweighed[n, :] *= i
    return reweighed
#end function

def main(argv):
    np_filename = argv[1]
    pwm_filename = argv[2]
    
    np_adjmtr = nmp.abs(nmp.loadtxt(argv[1]))
    pwm_adjmtr = nmp.loadtxt(argv[2])
    # new_adjmtr = reweigh_pwm_dot(pwm_adjmtr, np_adjmtr)
    #new_adjmtr = reweigh_pwm_corr(pwm_adjmtr, np_adjmtr)
    new_adjmtr = reweigh_pwm_cov(pwm_adjmtr, np_adjmtr)
    nmp.savetxt(sys.stdout, new_adjmtr)
#end function

if __name__ == "__main__":
    main(sys.argv)
