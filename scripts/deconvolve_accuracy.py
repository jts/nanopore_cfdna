import numpy as np

def accuracy(deconv_output):


    X = np.fromtxt(deconv_output)

    X_corr = np.corrcoef(X, rowvar=False)

    
