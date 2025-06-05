## IV. Indoor HAR Based on Point Cloud Matching ###

### A. Theory in Simple ###

The proposed method first converts the level sets corresponding to the micro-Doppler signature into contour point clouds using the MATLAB embedded contour() function. Secondly, the point cloud and the template are subjected to a similarity metric using the Mapper algorithm to obtain the final recognition results.

![Pointcloud_Matching](https://github.com/user-attachments/assets/506bf81d-1cd3-4fa8-9acd-e57a935a88e3)

Fig. 4. Schematic diagram of the proposed indoor HAR method based on point cloud topological structure similarity using Mapper algorithm.

### B. Codes Explanation (Folder: PointCloud_Matching) ###


#### 1. Mapper_Similarity ####

This function computes the topological similarity between two 2D point clouds using a simplified Mapper algorithm, measuring similarity via the Jaccard index of their Mapper graph edge sets.

**Input:** Two 2D point clouds `PC` and `PC_Class` of 2xN and 2xM matrices; Optional: integers `nx`, `ny` of grid squares, scalar `overlap_factor`.

**Output:** Scalar `similarity` represents the Jaccard similarity between edge sets.


#### 2. Select_Points_for_Columns ####

This function identifies the top points with the largest values in each column of a matrix, recording their row indices.

**Input:** 2D matrix `phi_1_Normalized` with size Estimation_Resolution x Estimation_Resolution, integer `Points_Num_Per_Column`.

**Output:** Matrix `points` in 2x(Points_Num_Per_Column*Estimation_Resolution) of row and column indices of the top points.


#### 3. contour ####

This function generates contour plots of a 2D matrix with options for specifying coordinates, levels, and line styles, supporting automatic or user-defined contour levels.

**Input:** Variable inputs: matrix `Z`, optional coordinates `X`, `Y`, levels `N` or `V`, axes handle, line specifications, and name-value pairs.

**Output:** Contour matrix `cout`, graphics handle `hand`.


#### 4. Wasserstein_Similarity ####

This function computes the 2-Wasserstein distance between two 2D point clouds using the optimal transport formulation with squared Euclidean distance costs.

**Input:** Two 2D point clouds `pointCloud1` and `pointCloud2` of 2xN and 2xM matrices.

**Output:** Scalar `similarity` represents the 2-Wasserstein distance.

### C. Datafiles Explanation (Folder: None) ###

None.
