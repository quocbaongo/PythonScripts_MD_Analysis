# Python script to extract minima structure from conformational landscape generated from PCA analysis

Principal Component Analysis (PCA) is a valuable technique in molecular dynamics (MD) simulation analysis, which enables the identification of different motion types of researched molecules within their molecular simulation trajectory.  PCA effectively filters out all the local and fast motions, while retaining only the global and collective (often slow) motions.  This can be achieved by constructing a covariance matrix of atomic fluctuations. Diagonalization of this matrix yields eigenvectors and eigenvalues, which describe the collective modes of motion. The eigenvectors associated with the largest eigenvalues are called "principal components" and correspond to the largest-amplitude collective motions. Projecting the trajectory onto the two principal components with the largest eigenvalues generates a conformational landscape that distinguishes different conformations of the simulated molecule during the simulation. An example of such a landscape is shown below. The landscape depicts sampled conformations of leptin interacting with the CRH2 domain of its receptor (the leptin receptor).

<p align="center">
  <img src="Screenshot from 2024-11-05 23-02-58.png" />
</p>
