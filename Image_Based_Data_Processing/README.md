## II. TWR HUMAN ECHO MODEL ##

### A. Theory in Simple ###

The proposed method first extracts the baseband signal of TWR echo by pulse compression, and then concatenates it along the slow time dimension, and the resulting image is a range-time map (RTM). The Doppler-time map (DTM) is obtained by summing all range bins of the RTM and doing the short time fourier transform (STFT) along the slow time dimension. The target image after clutter and noise suppression is obtained by doing Moving Target Indication Filtering (MTI) and Empirical Modal Decomposition (EMD) on both RTM and DTM, respectively.

![TWR_Echo_Model](https://github.com/user-attachments/assets/b2033710-afad-4ded-a4c6-5d4502d2a875)

Fig. 2. TWR human echo model and data processing.


### B. Codes Explanation (Folder: Image_Based_Data_Processing) ###


#### 1. EMD ####

This function denoises a radar RTM using EMD by discarding initial intrinsic mode functions (IMFs).

**Input:** 2D matrix `RTM`; Integer `num_discard` denotes the number of IMFs to discard.

**Output:** 2D matrix `denoised_RTM`.



#### 2. MTI ####

This function applies a simple MTI filter to a RTM by subtracting adjacent columns.

**Input:** 2D matrix `I`, here we usually use RTM for MTI processing.

**Output:** 2D matrix `filtered_I`.



#### 3. DTM_Generator ####

This function transforms a radar RTM into a DTM using the Short-Time Fourier Transform (STFT).

**Input:** 2D matrix `I`, RTM is used here; Optional parameters: `fs` for sampling frequency, `window` for window function, `noverlap` for overlap between windows, `nfft` for FFT points.

**Output:** 2D matrix `DTM`.


### C. Datafiles Explanation (Folder: Example_Datas) ###

Here includes four images as examples: one simulated RTM, one simulated DTM, one measured RTM, and one measured DTM.
