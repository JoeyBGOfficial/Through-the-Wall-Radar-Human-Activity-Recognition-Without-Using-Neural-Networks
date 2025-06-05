## III. MICRO-DOPPLER SIGNATURE EXTRACTION BASED ON ACM ##

### A. Theory in Simple ###

The micro-Doppler foreground and noise background centers using image corner detection are first estimated. Using the two centers as starting points, the micro-Doppler signature extraction is implemented based on multiphase Chan-Vese ACM model.

![ACM_Schematic](https://github.com/user-attachments/assets/0e152377-deaf-4eed-bc82-f599f9b8c97e)

Fig. 3. Schematic diagram of the proposed ACM-based micro-Doppler signature extraction method.

### B. Codes Explanation (Folder: ACM_Based_Micro-Doppler_Extraction) ###

#### 1. MTI  ####

This function implements static clutter cancellation on RTM.

**Input:** 2D matrix representing the RTM, with rows as range and columns as time.

**Output:** 2D matrix representing the MTI filtered image.
