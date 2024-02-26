# Umeyama

The Kabsch-Umeyama algorithm is a method for finding the optimal translation, rotation, and scaling that aligns two sets of points with minimum root-mean-square deviation (RMSD). It is named after Wolfgang Kabsch and Shinji Umeyama, who independently developed the algorithm for different applications. It is useful for comparing molecular and protein structures, point-set registration, and physics simulation.

Some of the main steps of the algorithm are:

* Calculate the centroids of the two sets of points and translate them to the origin.
* Compute the covariance matrix of the translated points and perform singular value decomposition on it.
* Determine the optimal rotation matrix and scale factor using the singular values and the sign of the determinant of the orthogonal matrices.
* Apply the translation, rotation, and scaling to the second set of points to align it with the first set.
