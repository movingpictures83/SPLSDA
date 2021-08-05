# SPLSDA
# Language: R
# Input: TXT (key, value pairs)
# Output: CSV
# Tested with: PluMA 1.1, R 4.0.0
# Dependency: mixOmics_6.12.1, RCurl_1.98.1.2, bitops_1.0.6, DiscriMiner_0.1.29

PluMA plugin to run Sparse Partial Least Squares Discriminant Analysis (Sparse PLSDA, LeCao et al 2011)
to determine the degree of separation between dataset(s), given values of observables in each set.

The plugin accepts input in the form of a text file, with (key, value) pairs.  The user provides in this file:

KEY         MEANING
---         -------
samples     CSV file containing a table where rows are samples and columns are observables
categories  TXT file containing sample categories (one per sample)
observables TXT file containing the name of each observable (one per observables)
targets     Categories of observables on which PLS-DA should be run

A plot will be produced, and in addition a CSV file of the (X, Y) coordinate values corresponding to principal coordinates one and two for every sample will be generated.
