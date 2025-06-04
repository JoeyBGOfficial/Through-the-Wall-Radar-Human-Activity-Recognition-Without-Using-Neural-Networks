%% Script for Generating Point Cloud Template Sets
% Author: JoeyBG.
% Improved By: JoeyBG.
% Affiliation: Beijing Institute of Technology.
% Date: 2025/5/28.
%
% Information:
% This script generates point cloud templates for four types of data: Simulated RTM, Simulated DTM,
% Measured RTM, and Measured DTM. For each data type, it processes images from predefined directories,
% extracts point cloud features using specific feature extraction functions, and saves the resulting
% point clouds as .mat files in corresponding output directories. The script iterates over a set of
% class names and processes a subset of images (1 to 300 with a step of 15) for each class. The output
% point clouds are organized by class and data type for use in subsequent classification tasks.
%
% The script assumes the existence of image files in the specified directories and relies on external
% feature extraction functions (e.g., Feature_Extraction_SimHRTM). It also loads a predefined colormap
% and class names for consistency with related scripts.
%
% Inputs:
%   - None (uses predefined paths and class names).
%
% Outputs:
%   - .mat files containing point cloud data, saved in the specified pointcloud_root directories.
%
% Dependencies:
%   - JoeyBG_CList.mat: Custom colormap preset.
%   - Class_Names.mat: Array of class names for organizing templates.
%   - Feature_Extraction_SimHRTM: Function for Simulated RTM point cloud extraction.
%   - Feature_Extraction_SimHDTM: Function for Simulated DTM point cloud extraction.
%   - Feature_Extraction_RWRTM: Function for Measured RTM point cloud extraction.
%   - Feature_Extraction_RWDTM: Function for Measured DTM point cloud extraction.
%
% See also: Feature_Extraction_SimHRTM, Feature_Extraction_SimHDTM, Feature_Extraction_RWRTM,
%           Feature_Extraction_RWDTM.

%% Initialization of Matlab Script
clear all;
close all;
clc;
disp('© Author: JoeyBG ©');
CList = load('JoeyBG_CList.mat').CList; % My favorite colormap preset.
CList_Flip = flipud(CList);
Class_Names = load("Class_Names.mat").Class_Names; 

%% Pointcloud Template Sets Generator for Simulated RTM
% Define the root paths for images and point clouds.
image_root = 'Image_Templates\SimH_RTM_Images\';
pointcloud_root = 'Pointcloud_Templates\SimH_RTM_Pointcloud\';

% Loop through each class in Class_Names.
for i = 1:length(Class_Names)
    class_name = Class_Names(i);
    
    % Construct the full paths for the image and pointcloud folders.
    image_folder = fullfile(image_root, class_name);
    pointcloud_folder = fullfile(pointcloud_root, class_name);
    
    % Create the pointcloud folder if it does not exist.
    if ~exist(pointcloud_folder, 'dir')
        mkdir(pointcloud_folder);
    end
    
    % Loop through each image from 1 to 300 with step 15.
    for j = 1:15:300
        image_file = fullfile(image_folder, [num2str(j) '.png']);
        pointcloud_file = fullfile(pointcloud_folder, [num2str(j) '.mat']);

        % Main process of pointcloud feature extraction.
        PointCloud = Feature_Extraction_SimHRTM(image_file); 
        
        % Save the point cloud to a .mat file.
        save(pointcloud_file, 'PointCloud');
    end
end

%% Pointcloud Template Sets Generator for Simulated DTM
% Define the root paths for images and point clouds.
image_root = 'Image_Templates\SimH_DTM_Images\';
pointcloud_root = 'Pointcloud_Templates\SimH_DTM_Pointcloud\';

% Loop through each class in Class_Names.
for i = 1:length(Class_Names)
    class_name = Class_Names(i);
    
    % Construct the full paths for the image and pointcloud folders.
    image_folder = fullfile(image_root, class_name);
    pointcloud_folder = fullfile(pointcloud_root, class_name);
    
    % Create the pointcloud folder if it does not exist.
    if ~exist(pointcloud_folder, 'dir')
        mkdir(pointcloud_folder);
    end
    
    % Loop through each image from 1 to 300 with step 15.
    for j = 1:15:300
        image_file = fullfile(image_folder, [num2str(j) '.png']);
        pointcloud_file = fullfile(pointcloud_folder, [num2str(j) '.mat']);

        % Main process of pointcloud feature extraction.
        PointCloud = Feature_Extraction_SimHDTM(image_file);
        
        % Save the point cloud to a .mat file.
        save(pointcloud_file, 'PointCloud');
    end
end

%% Pointcloud Template Sets Generator for Measured RTM
% Define the root paths for images and point clouds.
image_root = 'Image_Templates\RW_RTM_Images\';
pointcloud_root = 'Pointcloud_Templates\RW_RTM_Pointcloud\';

% Loop through each class in Class_Names.
for i = 1:length(Class_Names)
    class_name = Class_Names(i);
    
    % Construct the full paths for the image and pointcloud folders.
    image_folder = fullfile(image_root, class_name);
    pointcloud_folder = fullfile(pointcloud_root, class_name);
    
    % Create the pointcloud folder if it does not exist.
    if ~exist(pointcloud_folder, 'dir')
        mkdir(pointcloud_folder);
    end
    
    % Loop through each image from 1 to 300 with step 15.
    for j = 1:15:300
        image_file = fullfile(image_folder, [num2str(j) '.png']);
        pointcloud_file = fullfile(pointcloud_folder, [num2str(j) '.mat']);

        % Main process of pointcloud feature extraction.
        PointCloud = Feature_Extraction_RWRTM(image_file); 
        
        % Save the point cloud to a .mat file.
        save(pointcloud_file, 'PointCloud');
    end
end

%% Pointcloud Template Sets Generator for Simulated DTM
% Define the root paths for images and point clouds.
image_root = 'Image_Templates\RW_DTM_Images\';
pointcloud_root = 'Pointcloud_Templates\RW_DTM_Pointcloud\';

% Loop through each class in Class_Names.
for i = 1:length(Class_Names)
    class_name = Class_Names(i);
    
    % Construct the full paths for the image and pointcloud folders.
    image_folder = fullfile(image_root, class_name);
    pointcloud_folder = fullfile(pointcloud_root, class_name);
    
    % Create the pointcloud folder if it does not exist.
    if ~exist(pointcloud_folder, 'dir')
        mkdir(pointcloud_folder);
    end
    
    % Loop through each image from 1 to 300 with step 15.
    for j = 1:15:300
        image_file = fullfile(image_folder, [num2str(j) '.png']);
        pointcloud_file = fullfile(pointcloud_folder, [num2str(j) '.mat']);

        % Main process of pointcloud feature extraction.
        PointCloud = Feature_Extraction_RWDTM(image_file);
        
        % Save the point cloud to a .mat file.
        save(pointcloud_file, 'PointCloud');
    end
end