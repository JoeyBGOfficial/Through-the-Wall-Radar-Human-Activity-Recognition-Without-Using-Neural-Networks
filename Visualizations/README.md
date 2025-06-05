## V. Main Branch & Visualization ##

### A. Theory in Simple ###
For simulated RTM, simulated DTM, measured RTM, and measured DTM, we wrote different function scripts for the whole process of feature extraction and recognition. Main.m is used to achieve the inference prediction. The generation code of the template point cloud library and the code used to give all the visualized images in the paper are also open-sourced together.

![Simulated_Visualization](https://github.com/user-attachments/assets/38567527-71c4-42eb-8d46-34d9a8be968b)

Fig. 5. Simulated visualization results of the proposed method.

![Measured_Visualization](https://github.com/user-attachments/assets/d5a8d2b3-b9d6-4d47-9fb7-cea343d9052d)

Fig. 6. Measured visualization results of the proposed method.

### B. Codes Explanation (Folder: Root) ###


#### 1. Feature_Extraction_SimHRTM ####

This function deals with simulated RTM data, using a four-phase level set method for segmentation and creating a sparse point cloud from the extracted contour features.

**Input:** Path to the image data.

**Output:** Point cloud data of the extracted features.


#### 2. Feature_Extraction_SimHDTM ####

This function processes simulated DTM data, employing a four-phase level set method for image segmentation and generating a sparse point cloud from the extracted contour features.

**Input:** Path to the image data.

**Output:** Point cloud data representing the extracted features.


#### 3. Feature_Extraction_RWRTM ####

This function handles measured RTM data, utilizing a four-phase level set method for segmentation and producing a sparse point cloud from the contour of the extracted features.

**Input:** Path to the image data.

**Output:** Point cloud data of the extracted features.


#### 4. Feature_Extraction_RWDTM ####

This function processes measured DTM data, applying a four-phase level set method for image segmentation and generating a sparse point cloud from the contour features.

**Input:** Path to the image data.

**Output:** Point cloud data representing the extracted features.

#### 5. Main ####

This script serves as the primary interface for TWR HAR. It allows users to select data types (simulated or measured, RTM or DTM), performs feature extraction, and classifies the input based on similarity to template point clouds.

**Input:** User selection from a menu.

**Output:** Displays the predicted class name.


#### 6. Templates_Generator ####

This script generates template point clouds for simulated and measured RTM and DTM data. It processes images from predefined directories, extracts features, and saves the point clouds for use in classification.

**Input:** None, just use predefined paths.

**Output:** .mat files containing point cloud templates.


#### 7. Visualization_12_Activities ####

This script visualizes the feature extraction and recognition results for 12 activities across simulated and measured RTM and DTM data, displaying images without axes or labels.

**Input:** None, also use predefined paths.

**Output:** Visualized images saved to specified directories.


### C. Datafiles Explanation (Folder: Root, Visualizations) ###

#### 1. JoeyBG_CList.mat ####

My favorite colormap file used for generating figures in the paper.

#### 2. Class_Names.mat ####

Name strings of $12$ predefined classes of activities.

#### 3. Visualization Sub-Figures ####

High-resolution coordinate-free files for each subplot of the visualization experiments in the paper can be found in the Visualizations folder.
