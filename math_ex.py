import numpy as np

def tsum(x,y,cumu=0):
    npts = x.size
    xdif = x[1:]-x[:-1]
    yavg = (y[:-1] + y[1:])/2.
    if cumu==0: 
        sum = np.sum(np.append(np.array([0.]),xdif*yavg))
    else:
        sum = np.cumsum(np.append(np.array([0]),xdif*yavg))
    return sum
