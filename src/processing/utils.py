import numpy as np

alpha_svp_tzvp = 10.390
alpha_tzvpp_qzvpp = 7.880

beta_svp_tzvp = 2.400
beta_tzvpp_qzvpp = 2.970

def correct_hf(e1, e2, n1, n2, alpha):
    return e2 + (e2 - e1) / (np.exp(alpha * (np.sqrt(n2) - np.sqrt(n1))) - 1)

def correct_wft(e1, e2, n1, n2, beta):
    return (n1 ** beta * e1 - n2 ** beta * e2) / (n1 ** beta - n2 ** beta)