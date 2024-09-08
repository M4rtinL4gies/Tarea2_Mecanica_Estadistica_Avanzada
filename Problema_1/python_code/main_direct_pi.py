import numpy as np
from functions import make_direct_data

if __name__ == "__main__":
    runs = np.array([10**1, 10**2, 10**3, 10**4, 10**5, 10**6, 10**7, 10**8])
    n = 20
    make_direct_data(runs=runs[0:3], n=n)