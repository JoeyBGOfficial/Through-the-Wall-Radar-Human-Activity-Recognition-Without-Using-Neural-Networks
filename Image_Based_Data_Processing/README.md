## II. TWR HUMAN ECHO MODEL ##

### A. Theory in Simple ###

The proposed method first extracts the baseband signal of TWR echo by pulse compression, and then concatenates it along the slow time dimension, and the resulting image is a range-time map (RTM). The Doppler-time map (DTM) is obtained by summing all range bins of the RTM and doing the short time fourier transform (STFT) along the slow time dimension. The target image after clutter and noise suppression is obtained by doing Moving Target Indication Filtering (MTI) and Empirical Modal Decomposition (EMD) on both RTM and DTM, respectively.

![TWR_Echo_Model](https://github.com/user-attachments/assets/b2033710-afad-4ded-a4c6-5d4502d2a875)

Fig. 2. TWR human echo model and data processing.

### B. Codes Explanation (Folder: Image_Based_Data_Processing) ###

#### 1. MTI  ####

This function implements static clutter cancellation on RTM.

**Input:** 2D matrix representing the RTM, with rows as range and columns as time.

**Output:** 2D matrix representing the MTI filtered image.

#### 2. DTM_Generator  ####

This function implements DTM generation using STFT method implemented in MATLAB.

**Input:** 2D matrix representing the RTM; Sampling frequency; Window function; Number of overlapping samples between windows; Number of Fast Fourier Transform (FFT) points.

**Output:** 2D matrix representing the DTM, with rows as Doppler frequencies and columns as time segments.

#### 3. EMD ####

This function denoises radar RTM using Empirical Mode Decomposition (EMD). 

**Input:** 2D matrix representing the RTM; Number of initial intrinsic mode functions (IMFs) to discard as noise.

**Output:** 2D matrix representing the denoised RTM.

### C. Datafiles Explanation (Folder: Example_Datas) ###

Here includes four images as examples: one simulated RTM, one simulated DTM, one measured RTM, and one measured DTM.
