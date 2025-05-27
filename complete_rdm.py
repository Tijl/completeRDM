#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created May 2025

@author: tijl grootswagers
"""

import numpy as np

def complete_rdm(X, verbose=False):
    """
    Estimate missing values in a representational dissimilarity matrix (RDM).

    Parameters:
    X : np.ndarray
        Incomplete symmetric distance matrix containing np.nan for missing values.
    verbose : bool, optional
        If True, prints progress and estimation details during reconstruction.

    Returns:
    Y : np.ndarray
        Completed symmetric distance matrix with all entries filled (if the input has no isolated clusters)
    """
    X = np.array(X, dtype=np.float64)
    np.fill_diagonal(X, 0)  # Set diagonal to zero
    Y = X.copy()

    prev_missing = np.inf

    while True:
        missing_i, missing_j = np.where(np.triu(np.isnan(Y),1))
        nmissing = len(missing_i)

        if verbose:
            print(f"Reconstructing {nmissing} missing entries")

        # keep going until all missing entries are filled
        if nmissing == 0 or nmissing >= prev_missing:
            if verbose and nmissing > 0:
                print("No progress in reducing missing entries; stopping.")
            return Y

        prev_missing = nmissing

        for n in range(nmissing):
            i = missing_i[n]
            j = missing_j[n]

            ai = Y[i, :]
            bj = Y[j, :]

            # find reference entries
            known = (~np.isnan(ai)) & (~np.isnan(bj)) & (np.arange(len(ai)) != i) & (np.arange(len(ai)) != j)
            if np.any(known):
                a = ai[known]
                b = bj[known]

                # we estimate the missing distance using Pythagorean theorem, 
                # either as the square root of the sum of squares (right angle)
                # or the difference of squares of the known distances
                d1 = np.sqrt(a**2 + b**2)
                d2 = np.sqrt(np.abs(a**2 - b**2))
                # we take the median of all estimates
                d_est = np.median(np.concatenate([d1, d2]))

                # assign estimate to matrix
                Y[i, j] = d_est
                Y[j, i] = d_est

                if verbose:
                    print(f"{n+1}/{nmissing} Estimated ({i},{j}) and ({j},{i}) with {d_est:.4f} using {len(a)} references")


if __name__ == "__main__":
    # test case:
    # Create a symmetric matrix with a few missing values
    X = np.array([
        [0.0, 1.0, np.nan, 2.0],
        [1.0, 0.0, 1.5,   np.nan],
        [np.nan, 1.5, 0.0, 1.2],
        [2.0,   np.nan, 1.2, 0.0]
    ])

    print("Original matrix with NaNs:")
    print(X)

    # Complete the matrix
    Y = complete_rdm(X, verbose=True)

    print("\nCompleted matrix:")
    print(Y)

    # test 2, generate increasingly missing data, and plot the result
    from scipy.spatial.distance import pdist, squareform
    import matplotlib.pyplot as plt
    np.random.seed(1)  # For reproducibility

    # True condensed RDM from 100 random 2D points
    points = np.random.randn(100, 2)
    Xt = pdist(points)  # True distances (condensed)

    r = []
    for n in range(1, 81):  # 1% to 80% missing
        X = Xt.copy()
        n_missing = int(np.ceil(n * len(Xt) / 100))
        missing_idx = np.random.choice(len(X), n_missing, replace=False)
        X[missing_idx] = np.nan

        # Convert to square form and reconstruct
        Xsq = squareform(X)
        Ysq = complete_rdm(Xsq, verbose=False)
        Yt = squareform(Ysq)

        # Correlate only on originally known entries
        corr = np.corrcoef(Yt, Xt)[0, 1]
        r.append(corr)

    # Plot the results
    plt.figure(figsize=(6, 4))
    plt.plot(r, marker='o')
    plt.xlabel('% Missing Values')
    plt.ylabel('Reconstruction accuracy')
    plt.title('RDM Completion Accuracy vs Missing Data')
    plt.grid(True)
    plt.tight_layout()
    plt.ylim(((0,1)))
    plt.show()
