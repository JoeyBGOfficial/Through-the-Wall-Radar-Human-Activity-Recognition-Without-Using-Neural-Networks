## III. MICRO-DOPPLER SIGNATURE EXTRACTION BASED ON ACM ##

### A. Theory in Simple ###

The micro-Doppler foreground and noise background centers using image corner detection are first estimated. Using the two centers as starting points, the micro-Doppler signature extraction is implemented based on multiphase Chan-Vese ACM model.

![ACM_Schematic](https://github.com/user-attachments/assets/0e152377-deaf-4eed-bc82-f599f9b8c97e)

Fig. 3. Schematic diagram of the proposed ACM-based micro-Doppler signature extraction method.

### B. Codes Explanation (Folder: ACM_Based_Micro-Doppler_Extraction) ###

#### 1. CURVATURE_CV  ####

This function computes the curvature of a 2D matrix f using different finite difference schemes specified by the input diff_scheme.

**Input:** 2D matrix representing the level set function or image; Integer (0, 1, or 2) selecting the difference scheme to use.

**Output:** 2D matrix representing the approximated curvature of the input function.

#### 2. Corner_Representation  ####

This function detects Scale-Invariant Feature Transform (SIFT) corner points in a grayscale image, computes the centroid of these points, and identifies the pixels farthest and nearest to this centroid.

**Input:** 2D matrix representing the grayscale image.

**Output:** Farest pixel; Nearest pixel.

#### 3. Delta  ####

This function computes an approximation to the Dirac delta function.

**Input:** 2D matrix representing the level set function; Scalar value determining the width of the delta function approximation.

**Output:** 2D matrix approximating the Dirac delta function applied to phi.


