%% Script for Through-the-Wall Radar Human Activity Recognition Without Using Neural Networks
% Author: JoeyBG.
% Improved By: JoeyBG.
% Affiliation: Beijing Institute of Technology.
% Date: 2025/5/30.
%
% Information:
% This script performs feature extraction on an input image and classifies it by comparing the
% extracted point cloud features to predefined template point clouds. The user selects an option
% from a menu (Simulated RTM, Simulated DTM, Measured RTM, or Measured DTM), which determines
% the feature extraction method and the template directory used for classification. The script
% calculates similarity scores between the input point cloud and templates across 12 classes,
% predicting the class with the highest similarity.
%
% The script relies on external functions for feature extraction (e.g., Feature_Extraction_SimHRTM)
% and similarity computation (Mapper_Similarity), which are not defined here. It assumes the
% existence of predefined data files for class names and templates.
%
% Inputs:
%   - User selection from the menu (integer: 1 to 4).
%
% Outputs:
%   - Displays the predicted class name based on similarity to template point clouds.
%
% Dependencies:
%   - JoeyBG_CList.mat: Custom colormap preset.
%   - Class_Names.mat: List of class names for classification.
%   - Feature_Extraction_* functions: External functions for point cloud generation.
%   - Mapper_Similarity: External function for similarity computation.
%   - Template directories: Containing .mat files with point cloud templates.
%
% See also: Feature_Extraction_SimHRTM, Feature_Extraction_SimHDTM, Feature_Extraction_RWRTM,
%           Feature_Extraction_RWDTM, Mapper_Similarity.

%% Initialization of Matlab Script
clear all;
close all;
clc;
disp('© Author: JoeyBG ©');
CList = load('JoeyBG_CList.mat').CList; % My favorite colormap preset.
CList_Flip = flipud(CList);
Class_Names = load("Class_Names.mat").Class_Names; 

%% Definitions
File_Path = 'Example_Datas/';                               % The folder to store the input images for recognition inference.
File_Name = '10.png';                                       % Name of the input data.
Data_Path = strcat(File_Path,File_Name);                    % Complete path of the pre-defined input data.
nx = 100;                                                   % Integer, number of squares along the x-axis in the covering grid (default: 100).
ny = 100;                                                   % Integer, number of squares along the y-axis in the covering grid (default: 100).
overlap_factor = 1.5;                                       % Scalar, factor determining square size relative to step size for overlap (default: 1.5).

%% Create A Menu
options = {'Simulated RTM', 'Simulated DTM', 'Measured RTM', 'Measured DTM'};
choice = menu('Select an option:', options);

switch choice
    %% Feautre Extraction Based on Simulated RTM
    case 1
        % Generate Feature Pointcloud.
        PC = Feature_Extraction_SimHRTM(Data_Path);

        % Readin template pointclouds.
        root_dir = 'PointCloud_Templates/SimH_RTM_Pointcloud/';
        similarity_sums = zeros(12, 1);
        
        % Calculate similarity for all classes.
        for i = 1:12
            class_name = Class_Names(i);
            class_dir = fullfile(root_dir, class_name);
            mat_files = dir(fullfile(class_dir, '*.mat'));
            total_similarity = 0;
            for j = 1:length(mat_files)
                mat_file_path = fullfile(class_dir, mat_files(j).name);
                data = load(mat_file_path);
                PC_Template = data.PointCloud;
                sim = Mapper_Similarity(PC, PC_Template,nx,ny,overlap_factor);
                total_similarity = total_similarity + sim;
            end
            similarity_sums(i) = total_similarity;
        end

        % Find the best-match class between the input pointcloud and the templates.
        [~, max_index] = max(similarity_sums);
        max_class_name = Class_Names(max_index);
        disp(strcat('Prediction:',{32},max_class_name));
        
    case 2
        % Generate Feature Pointcloud.
        PC = Feature_Extraction_SimHDTM(Data_Path);

        % Readin template pointclouds.
        root_dir = 'PointCloud_Templates/SimH_DTM_Pointcloud/';
        similarity_sums = zeros(12, 1);
        
        % Calculate similarity for all classes.
        for i = 1:12
            class_name = Class_Names(i);
            class_dir = fullfile(root_dir, class_name);
            mat_files = dir(fullfile(class_dir, '*.mat'));
            total_similarity = 0;
            for j = 1:length(mat_files)
                mat_file_path = fullfile(class_dir, mat_files(j).name);
                data = load(mat_file_path);
                PC_Template = data.PointCloud;
                sim = Mapper_Similarity(PC, PC_Template,nx,ny,overlap_factor);
                total_similarity = total_similarity + sim;
            end
            similarity_sums(i) = total_similarity;
        end

        % Find the best-match class between the input pointcloud and the templates.
        [~, max_index] = max(similarity_sums);
        max_class_name = Class_Names(max_index);
        disp(strcat('Prediction:',{32},max_class_name));
        
    case 3
        % Generate Feature Pointcloud.
        PC = Feature_Extraction_RWRTM(Data_Path);

        % Readin template pointclouds.
        root_dir = 'PointCloud_Templates/RW_RTM_Pointcloud/';
        similarity_sums = zeros(12, 1);
        
        % Calculate similarity for all classes.
        for i = 1:12
            class_name = Class_Names(i);
            class_dir = fullfile(root_dir, class_name);
            mat_files = dir(fullfile(class_dir, '*.mat'));
            total_similarity = 0;
            for j = 1:length(mat_files)
                mat_file_path = fullfile(class_dir, mat_files(j).name);
                data = load(mat_file_path);
                PC_Template = data.PointCloud;
                sim = Mapper_Similarity(PC, PC_Template,nx,ny,overlap_factor);
                total_similarity = total_similarity + sim;
            end
            similarity_sums(i) = total_similarity;
        end

        % Find the best-match class between the input pointcloud and the templates.
        [~, max_index] = max(similarity_sums);
        max_class_name = Class_Names(max_index);
        disp(strcat('Prediction:',{32},max_class_name));

    case 4
        % Generate Feature Pointcloud.
        PC = Feature_Extraction_RWDTM(Data_Path);

        % Readin template pointclouds.
        root_dir = 'PointCloud_Templates/RW_DTM_Pointcloud/';
        similarity_sums = zeros(12, 1);
        
        % Calculate similarity for all classes.
        for i = 1:12
            class_name = Class_Names(i);
            class_dir = fullfile(root_dir, class_name);
            mat_files = dir(fullfile(class_dir, '*.mat'));
            total_similarity = 0;
            for j = 1:length(mat_files)
                mat_file_path = fullfile(class_dir, mat_files(j).name);
                data = load(mat_file_path);
                PC_Template = data.PointCloud;
                sim = Mapper_Similarity(PC, PC_Template,nx,ny,overlap_factor);
                total_similarity = total_similarity + sim;
            end
            similarity_sums(i) = total_similarity;
        end

        % Find the best-match class between the input pointcloud and the templates.
        [~, max_index] = max(similarity_sums);
        max_class_name = Class_Names(max_index);
        disp(strcat('Prediction:',{32},max_class_name));
        
    otherwise
        % Handle unexpected cases.
        disp('Invalid selection! Exiting.');
        return;

end