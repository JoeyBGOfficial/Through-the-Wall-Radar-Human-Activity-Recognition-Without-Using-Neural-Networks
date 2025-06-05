## III. MICRO-DOPPLER SIGNATURE EXTRACTION BASED ON ACM ##

### A. Theory in Simple ###

The micro-Doppler foreground and noise background centers using image corner detection are first estimated. Using the two centers as starting points, the micro-Doppler signature extraction is implemented based on multiphase Chan-Vese ACM model.

![ACM_Schematic](https://github.com/user-attachments/assets/0e152377-deaf-4eed-bc82-f599f9b8c97e)

Fig. 3. Schematic diagram of the proposed ACM-based micro-Doppler signature extraction method.

### B. Codes Explanation (Folder: ACM_Based_Micro-Doppler_Extraction) ###

#### 1. backward_gradient ####

This function computes the backward differences of a 2D matrix `f`, approximating partial derivatives along rows and columns.

**Input:** 2D matrix `f` of level set function or image.

**Output:** Two 2D matrices: `bdy` for column-wise backward differences, `bdx` for row-wise backward differences.



#### 2. Corner_Representation ####

This function detects SIFT corner points in a grayscale image, computes their centroid, and finds the pixels farthest and nearest to it.

**Input:** 2D matrix `I` of grayscale image, scalar `Corner_Threshold_Ratio`.

**Output:** Vectors `farPixel` and `nearPixel` (1x2 vector), normalized image `II_Normalized`, corner locations, and coordinates of farthest/nearest pixels.



#### 3. forward_gradient ####

This function computes the forward differences of a 2D matrix `f`, approximating partial derivatives along rows and columns.

**Input:** 2D matrix `f`.

**Output:** Two 2D matrices: `fdy` for column-wise forward differences, `fdx` for row-wise forward differences.



#### 4. EVOLUTION_4PHASE ####

This function evolves two level set functions for four-phase image segmentation using curvature and data fitting terms over multiple iterations.

**Input:** Image `I`, initial level sets `phi0`, parameters `nu`, `lambda_1`, `lambda_2`, `delta_t`, `epsilon`, and `numIter`.

**Output:** Evolved level set functions `phi` in 3D matrix.


#### 5. EVOLUTION_4PHASE_DR ####

This function extends `EVOLUTION_4PHASE` by adding a distance regularization term to maintain level set regularity.

**Input:** Same as `EVOLUTION_4PHASE`, plus `mu` for regularization weight.

**Output:** Evolved level set functions `phi` in 3D matrix.



#### 6. Delta ####

This function computes a smooth approximation of the Dirac delta function for level set methods, focused near the zero level set.

**Input:** 2D matrix `phi` represnets the level set function, scalar `epsilon` for width parameter.

**Output:** 2D matrix `Delta_h` approximating the Dirac delta function.



#### 7. get_contour ####

This function generates contour points along the perimeter of a rectangular region inset by `margin` pixels in a grid.

**Input:** Image `I` (unused), integers `nrow`, `ncol` of grid size, `margin` of inset distance.

**Output:** Arrays `xcontour` and `ycontour` represent the coordinates of contour points.



#### 8. Heaviside ####

This function computes a smooth approximation of the Heaviside step function using the arctangent function.

**Input:** 2D matrix or scalar `phi`, scalar `epsilon` for smoothness.

**Output:** Smooth Heaviside function values `H`.


#### 9. initial_sdf2circle ####

This function initializes two level set functions as circles centered at specified pixels in a grid.

**Input:** Grid size `nrow`, `ncol`, unused `ic`, `jc`, `fun_n`, radius factor `r`, center coordinates `far_pixel`, `near_pixel`.

**Output:** 3D matrix `f` containing two level set functions.


#### 10. quadrifit ####

This function computes constants `C` that best fit an image `U` across four regions defined by two level set functions, using a smooth Heaviside approximation.

**Input:** Level sets `phi` in 3D matrix, image `U`, `epsilon`, number of functions `fun_n` .

**Output:** Vector `C` regional constants, intermediate factors `mult` in 4D matrix.


#### 11. CURVATURE_CV ####

This function computes the curvature of a 2D matrix `f` using different finite difference schemes specified by `diff_scheme`.

**Input:** 2D matrix `f` of level set function or image, Integer `diff_scheme` (0, 1, or 2) selecting the difference scheme.

**Output:** 2D matrix `K` representing the approximated curvature of `f`.

### C. Dataset Explanation (Folder: None) ###

None.
