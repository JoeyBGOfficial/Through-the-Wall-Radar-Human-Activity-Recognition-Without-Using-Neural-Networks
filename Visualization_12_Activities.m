%% Script for Visualization Experiments in the Paper
% Author: JoeyBG.
% Improved By: JoeyBG.
% Affiliation: Beijing Institute of Technology.
% Date: 2025/6/4.
%
% Information:
% The script's function is to display images for RTM, DTM, feature extraction, and recognition. 
% Simulated and measured images are displayed separately. 
% The first row of images shows the results of RTM for 12 activities after corner detection and nearest-farthest point estimation, 
% the second row shows feature extraction, and the third row shows the micro-Doppler level set. 
% The fourth row shows the results of DTM for 12 activities after corner detection and nearest-farthest point calculation, 
% the fifth row shows feature extraction, and the sixth row shows the micro-Doppler level set. 
% The format for measured and simulated images is consistent.
%
% Inputs:
%   - None (uses predefined paths and class names).
%
% Outputs:
%   - Images without any axis and labels mentioned in the information above.
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

%% The Copying of Template Images into Visualization Folder
Number_of_Images = 87; % Pre-defined number of image in all folders used for visualization.
main_folders = {'SimH_RTM_Images', 'SimH_DTM_Images', 'RW_RTM_Images', 'RW_DTM_Images'};
base_input_path = 'Image_Templates\';
base_output_path = 'Visualizations\';

% Loop through each main folder.
for i = 1:length(main_folders)
    main_folder = main_folders{i};
    
    % Construct input and output paths for the current main folder.
    input_main_path = fullfile(base_input_path, main_folder);
    output_main_path = fullfile(base_output_path, main_folder);
    
    % Create the output main folder if it doesn't exist.
    if ~exist(output_main_path, 'dir')
        mkdir(output_main_path);
    end
    
    % Loop through each subfolder in Class_Names.
    for j = 1:length(Class_Names)
        sub_folder = Class_Names(j);
        input_sub_path = fullfile(input_main_path, sub_folder);
        img_path = fullfile(input_sub_path, strcat(num2str(Number_of_Images),'.png'));
        new_img_name = strcat(sub_folder,'.png');
        output_img_path = fullfile(output_main_path, new_img_name);
        copyfile(img_path, output_img_path);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generating the Visualizations of Simulated RTM Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for Image_Name = 1:length(Class_Names)
    Current_Image_Name = strcat(Class_Names(Image_Name),'.png'); % Define the name of the current image.
    Data_Path = fullfile('Visualizations','SimH_RTM_Images',Current_Image_Name); % Define the path of the current image.
    Storing_Path = fullfile('Visualizations','SimH_RTM_12_Activities');

    %% Re-Initialization of Matlab Script
    % clear all;
    close all;
    % clc;
    % CList = load('JoeyBG_CList.mat').CList; % My favorite colormap preset.
    % CList_Flip = flipud(CList);
    % Class_Names = load("Class_Names.mat").Class_Names; 

    %% Parameter Definition
    delta_t = .1;                                   % Time step for evolution of level set functions.
    lambda_1 = 1;                                   % Weight for data fitting term.
    lambda_2 = 1;                                   % Weight for data fitting term.
    h = 1;                                          % Spatial step size.
    epsilon = 1;                                    % Parameter for Heaviside function approximation
    nu = .001 * 255 * 255;                          % Weight for curvature term, scaled by image intensity range.
    fun_n = 2;                                      % Number of level sets (**Fixed in this script with value 2 for four-phase segmentation**).
    mu = 0.5;                                       % Coefficient for the distance regularizing term.
    numIter = 20;                                   % The number of iterations per evolution step.
    evolveStep = 70;                                % Evolve the level set functions over 'evolveStep' steps.
    r = 64;                                         % Radius of initial circles for level set functions.
    Corner_Threshold_Ratio = 0.3;                   % Threshold cutting ratio of image preprocessing in corner representation.
    Estimation_Resolution = 256;                    % Resize scale of images in this work.
    Max_Range = 4;                                  % Max range bin of RTM, unit: m.
    Max_Doppler = 20;                               % Max Doppler bin of DTM, unit: Hz.
    Max_Time = 4;                                   % Max slow time of both RTM and DTM, unit: s.
    
    %% Readin the Datas for Recognition
    % Load and preprocess the input image.
    Img = imread(Data_Path);  
    U = double(Img(:,:,2));      
    U = imresize(U, [Estimation_Resolution Estimation_Resolution]); 
    % Detect corner points for nearest / farest initialization.
    [far_pixel, near_pixel, II_Normalized, locations, farRow, farCol, nearRow, nearCol] = Corner_Representation(U,Corner_Threshold_Ratio);

    % Plot the results of corner representation.
    figure(1);
    set(gcf, 'Position', [50, 450, 400, 500]);  % Set window size to 4:5 ratio.
    imagesc(1-II_Normalized);
    
    % Set up views and labels of ACM numerical solving visualization.
    colormap(CList);
    xticks([0 Estimation_Resolution/4 Estimation_Resolution/2 Estimation_Resolution*3/4 Estimation_Resolution]);
    xticklabels({'','','','',''});
    yticks([0 Estimation_Resolution/4 Estimation_Resolution/2 Estimation_Resolution*3/4 Estimation_Resolution]);
    yticklabels({'','','','',''});
    set(gca,'FontName','TsangerYuMo W03','FontSize',14,'Fontweight','normal');
    % xlabel('Time (s)','FontSize',15,'FontName',"TsangerYuMo W03","FontWeight","normal");
    % ylabel('Range (m)','FontSize',15,'FontName',"TsangerYuMo W03","FontWeight","normal");
    % title('Corner Representation',"FontSize",16,"FontName","TsangerYuMo W03","FontWeight","bold");    
    hold on;
    plot(locations(:,1), locations(:,2), 'o', 'Color', '#BEE9FF', 'MarkerSize', 10, 'LineWidth', 2);
    % plot(centroid(1), centroid(2), 'co', 'MarkerSize', 10, 'LineWidth', 2);
    % text(centroid(1) + 10, centroid(2), 'Centroid', 'FontName', "TsangerYuMo W03", 'Color', 'b', 'FontSize', 14);
    plot(farCol, farRow, '*', 'Color', [0.25 0.25 0.25], 'MarkerSize', 18, 'LineWidth', 2);
    % text(farCol + 10, farRow, 'Farthest', 'FontName', "TsangerYuMo W03", 'Color', 'g', 'FontSize', 14);
    plot(nearCol, nearRow, '*', 'Color', [1 1 1], 'MarkerSize', 18, 'LineWidth', 2);
    % text(nearCol + 10, nearRow, 'Nearest', 'FontName', "TsangerYuMo W03", 'Color', 'y', 'FontSize', 14);
    % title('Corner Representation',"FontSize",16,"FontName","TsangerYuMo W03","FontWeight","bold");
    hold off;
    pause(0.1);
    exportgraphics(gcf,strcat(Storing_Path,'\Corner_SimH_RTM_',Current_Image_Name),'Resolution', 800);
    
    % Determine image dimensions.
    [nrow, ncol] = size(U);  
    
    %% ACM-Based Micro-Doppler Signature Extraction
    % Initialize the level set functions.
    ic = nrow / 2;  % Center row index.
    jc = ncol / 2;  % Center column index.
    phi = initial_sdf2circle(nrow, ncol, ic, jc, r, fun_n, far_pixel, near_pixel);  % Initialize two level set functions.
    
    % Prepare image data for computations.
    I = double(U);  
    I_max = max(max(I));
    I_min = min(min(I));
    I_Normalized = (I-I_min)/(I_max-I_min);
    
    % Visualize the initial level set functions.
    figure(2);       
    set(gcf, 'Position', [500, 450, 400, 500]);  % Set window size to 4:5 ratio.
    imagesc(1-I_Normalized);
    
    % Set up views and labels of ACM numerical solving visualization.
    colormap(CList);
    xticks([0 Estimation_Resolution/4 Estimation_Resolution/2 Estimation_Resolution*3/4 Estimation_Resolution]);
    xticklabels({'','','','',''});
    yticks([0 Estimation_Resolution/4 Estimation_Resolution/2 Estimation_Resolution*3/4 Estimation_Resolution]);
    yticklabels({'','','','',''});
    set(gca,'FontName','TsangerYuMo W03','FontSize',14,'Fontweight','normal');
    % xlabel('Time (s)','FontSize',15,'FontName',"TsangerYuMo W03","FontWeight","normal");
    % ylabel('Range (m)','FontSize',15,'FontName',"TsangerYuMo W03","FontWeight","normal");
    % title('RTM for ACM Solving',"FontSize",16,"FontName","TsangerYuMo W03","FontWeight","bold"); 

    % Contour visualization of the level set phi_1 input to the function.
    [Contour_Edge,Contour_Obj] = contour(phi(:,:,1),[0 0],'Color',[1 1 1],'LineWidth',2);
    % Set up views and labels of ACM numerical solving visualization.
    colormap(CList);
    xticks([0 Estimation_Resolution/4 Estimation_Resolution/2 Estimation_Resolution*3/4 Estimation_Resolution]);
    xticklabels({'','','','',''});
    yticks([0 Estimation_Resolution/4 Estimation_Resolution/2 Estimation_Resolution*3/4 Estimation_Resolution]);
    yticklabels({'','','','',''});
    set(gca,'FontName','TsangerYuMo W03','FontSize',14,'Fontweight','normal');
    % xlabel('Time (s)','FontSize',15,'FontName',"TsangerYuMo W03","FontWeight","normal");
    % ylabel('Range (m)','FontSize',15,'FontName',"TsangerYuMo W03","FontWeight","normal");
    % title('ACM Feature Extraction',"FontSize",16,"FontName","TsangerYuMo W03","FontWeight","bold");

    % Contour visualization of the level set phi_2 input to the function.
    contour(phi(:,:,2),[0 0],'Color',[0.25 0.25 0.25],'LineWidth',2);
    % Set up views and labels of ACM numerical solving visualization.
    colormap(CList);
    xticks([0 Estimation_Resolution/4 Estimation_Resolution/2 Estimation_Resolution*3/4 Estimation_Resolution]);
    xticklabels({'','','','',''});
    yticks([0 Estimation_Resolution/4 Estimation_Resolution/2 Estimation_Resolution*3/4 Estimation_Resolution]);
    yticklabels({'','','','',''});
    set(gca,'FontName','TsangerYuMo W03','FontSize',14,'Fontweight','normal');
    % xlabel('Time (s)','FontSize',15,'FontName',"TsangerYuMo W03","FontWeight","normal");
    % ylabel('Range (m)','FontSize',15,'FontName',"TsangerYuMo W03","FontWeight","normal");
    % title('ACM Feature Extraction',"FontSize",16,"FontName","TsangerYuMo W03","FontWeight","bold");
    
    % Main loop for evolving the level set functions.
    for k = 1:evolveStep
        phi = EVOLUTION_4PHASE_DR(I, phi, nu, lambda_1, lambda_2, mu, delta_t, epsilon, numIter);  % Update the level set functions.
        if mod(k, 2) == 0  
            pause(.01);     
            imagesc(1-I_Normalized);  
            colormap(gray);  

            hold on;            
            phi_1 = phi(:,:,1); % Extract the first level set function.
            phi_2 = phi(:,:,2); % Extract the second level set function.

            % Contour visualization of the level set phi_1 input to the function.
            [Contour_Edge,Contour_Obj] = contour(phi_1,[0 0],'Color',[1 1 1],'LineWidth',2);
            % Set up views and labels of ACM numerical solving visualization.
            colormap(CList);
            xticks([0 Estimation_Resolution/4 Estimation_Resolution/2 Estimation_Resolution*3/4 Estimation_Resolution]);
            xticklabels({'','','','',''});
            yticks([0 Estimation_Resolution/4 Estimation_Resolution/2 Estimation_Resolution*3/4 Estimation_Resolution]);
            yticklabels({'','','','',''});
            set(gca,'FontName','TsangerYuMo W03','FontSize',14,'Fontweight','normal');
            % xlabel('Time (s)','FontSize',15,'FontName',"TsangerYuMo W03","FontWeight","normal");
            % ylabel('Range (m)','FontSize',15,'FontName',"TsangerYuMo W03","FontWeight","normal");
            % title('ACM Feature Extraction',"FontSize",16,"FontName","TsangerYuMo W03","FontWeight","bold");
        
            % Contour visualization of the level set phi_2 input to the function.
            contour(phi_2,[0 0],'Color',[0.25 0.25 0.25],'LineWidth',2);
            % Set up views and labels of ACM numerical solving visualization.
            colormap(CList);
            xticks([0 Estimation_Resolution/4 Estimation_Resolution/2 Estimation_Resolution*3/4 Estimation_Resolution]);
            xticklabels({'','','','',''});
            yticks([0 Estimation_Resolution/4 Estimation_Resolution/2 Estimation_Resolution*3/4 Estimation_Resolution]);
            yticklabels({'','','','',''});
            set(gca,'FontName','TsangerYuMo W03','FontSize',14,'Fontweight','normal');
            % xlabel('Time (s)','FontSize',15,'FontName',"TsangerYuMo W03","FontWeight","normal");
            % ylabel('Range (m)','FontSize',15,'FontName',"TsangerYuMo W03","FontWeight","normal");
            % title('ACM Feature Extraction',"FontSize",16,"FontName","TsangerYuMo W03","FontWeight","bold");

            hold off;          
        end
    end
    pause(0.1);
    exportgraphics(gcf,strcat(Storing_Path,'\ACM_SimH_RTM_',Current_Image_Name),'Resolution', 800);
    
    %% Visualization of Feature Extraction
    % Normalization of two level sets.
    phi_1_max = max(max(phi_1));
    phi_1_min = min(min(phi_1));
    phi_1_Normalized = (phi_1-phi_1_min)/(phi_1_max-phi_1_min);
    phi_2_max = max(max(phi_2));
    phi_2_min = min(min(phi_2));
    phi_2_Normalized = (phi_2-phi_2_min)/(phi_2_max-phi_2_min);
    
    % Visualize the final level set functions as 3D meshes.
    figure(3);      
    set(gcf, 'Position', [950, 450, 400, 500]);  % Set window size to 4:5 ratio.
    mesh(flip(phi_1_Normalized));             
    
    % Set up views and labels of phi_1 visualization.
    colormap(CList_Flip);
    view(45,25); % Set to 3D view, 45 degree viewing azimuth angle, 45 degree viewing pitch angle.
    zticks([0 1]);
    set(gca,'zticklabel',[]);
    zticklabels({'',''});
    xticks([0 Estimation_Resolution/4 Estimation_Resolution/2 Estimation_Resolution*3/4 Estimation_Resolution]);
    xticklabels({'','','','',''});
    yticks([0 Estimation_Resolution/4 Estimation_Resolution/2 Estimation_Resolution*3/4 Estimation_Resolution]);
    yticklabels({'','','','',''});
    set(gca,'FontName','TsangerYuMo W03','FontSize',14,'Fontweight','normal');
    % xlabel('Time (s)','FontSize',15,'FontName',"TsangerYuMo W03","FontWeight","normal");
    % ylabel('Range (m)','FontSize',15,'FontName',"TsangerYuMo W03","FontWeight","normal");
    % zlabel('Amp','FontSize',15,'FontName',"TsangerYuMo W03","FontWeight","normal");
    % title('Micro-Doppler Signature',"FontSize",16,"FontName","TsangerYuMo W03","FontWeight","bold");
    pause(0.1);
    exportgraphics(gcf,strcat(Storing_Path,'\mD_SimH_RTM_',Current_Image_Name),'Resolution', 800);
    
    figure(4);   
    set(gcf, 'Position', [1400, 450, 400, 500]);  % Set window size to 4:5 ratio.
    mesh(flip(phi_2_Normalized)); 
    
    % Set up views and labels of phi_2 visualization.
    colormap(CList_Flip);
    view(45,25); % Set to 3D view, 45 degree viewing azimuth angle, 45 degree viewing pitch angle.
    zticks([0 1]);
    set(gca,'zticklabel',[]);
    zticklabels({'',''});
    xticks([0 Estimation_Resolution/4 Estimation_Resolution/2 Estimation_Resolution*3/4 Estimation_Resolution]);
    xticklabels({'','','','',''});
    yticks([0 Estimation_Resolution/4 Estimation_Resolution/2 Estimation_Resolution*3/4 Estimation_Resolution]);
    yticklabels({'','','','',''});
    set(gca,'FontName','TsangerYuMo W03','FontSize',14,'Fontweight','normal');
    % xlabel('Time (s)','FontSize',15,'FontName',"TsangerYuMo W03","FontWeight","normal");
    % ylabel('Range (m)','FontSize',15,'FontName',"TsangerYuMo W03","FontWeight","normal");
    % zlabel('Amp','FontSize',15,'FontName',"TsangerYuMo W03","FontWeight","normal");
    % title('Noise Background',"FontSize",16,"FontName","TsangerYuMo W03","FontWeight","bold");
    
    %% PointCloud Generation
    % PointCloud = Select_Points_for_Columns(phi_1_Normalized,1);
    PointCloud = Contour_Edge;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generating the Visualizations of Simulated DTM Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for Image_Name = 1:length(Class_Names)
    Current_Image_Name = strcat(Class_Names(Image_Name),'.png'); % Define the name of the current image.
    Data_Path = fullfile('Visualizations','SimH_DTM_Images',Current_Image_Name); % Define the path of the current image.
    Storing_Path = fullfile('Visualizations','SimH_DTM_12_Activities');
    
    %% Re-Initialization of Matlab Script
    % clear all;
    close all;
    % clc;
    % CList = load('JoeyBG_CList.mat').CList; % My favorite colormap preset.
    % CList_Flip = flipud(CList);
    % Class_Names = load("Class_Names.mat").Class_Names; 

    %% Parameter Definition
    delta_t = .1;                                   % Time step for evolution of level set functions.
    lambda_1 = 1;                                   % Weight for data fitting term.
    lambda_2 = 1;                                   % Weight for data fitting term.
    h = 1;                                          % Spatial step size.
    epsilon = 1;                                    % Parameter for Heaviside function approximation
    nu = .001 * 255 * 255;                          % Weight for curvature term, scaled by image intensity range.
    fun_n = 2;                                      % Number of level sets (**Fixed in this script with value 2 for four-phase segmentation**).
    mu = 0.5;                                       % Coefficient for the distance regularizing term.
    numIter = 20;                                   % The number of iterations per evolution step.
    evolveStep = 70;                                % Evolve the level set functions over 'evolveStep' steps.
    r = 64;                                         % Radius of initial circles for level set functions.
    Corner_Threshold_Ratio = 0.3;                   % Threshold cutting ratio of image preprocessing in corner representation.
    Estimation_Resolution = 256;                    % Resize scale of images in this work.
    Max_Range = 4;                                  % Max range bin of RTM, unit: m.
    Max_Doppler = 20;                               % Max Doppler bin of DTM, unit: Hz.
    Max_Time = 4;                                   % Max slow time of both RTM and DTM, unit: s.
    
    %% Readin the Datas for Recognition
    % Load and preprocess the input image.
    Img = imread(Data_Path);  
    U = double(Img(:,:,2));      
    U = imresize(U, [Estimation_Resolution Estimation_Resolution]); 
    % Detect corner points for nearest / farest initialization.
    [far_pixel, near_pixel, II_Normalized, locations, farRow, farCol, nearRow, nearCol] = Corner_Representation(U,Corner_Threshold_Ratio);

    % Plot the results of corner representation.
    figure(1);
    set(gcf, 'Position', [50, 450, 400, 500]);  % Set window size to 4:5 ratio.
    imagesc(1-II_Normalized);
    
    % Set up views and labels of ACM numerical solving visualization.
    colormap(CList);
    xticks([0 Estimation_Resolution/4 Estimation_Resolution/2 Estimation_Resolution*3/4 Estimation_Resolution]);
    xticklabels({'','','','',''});
    yticks([0 Estimation_Resolution/4 Estimation_Resolution/2 Estimation_Resolution*3/4 Estimation_Resolution]);
    yticklabels({'','','','',''});
    set(gca,'FontName','TsangerYuMo W03','FontSize',14,'Fontweight','normal');
    % xlabel('Time (s)','FontSize',15,'FontName',"TsangerYuMo W03","FontWeight","normal");
    % ylabel('Doppler (Hz)','FontSize',15,'FontName',"TsangerYuMo W03","FontWeight","normal");
    % title('Corner Representation',"FontSize",16,"FontName","TsangerYuMo W03","FontWeight","bold");    
    hold on;
    plot(locations(:,1), locations(:,2), 'o', 'Color', '#BEE9FF', 'MarkerSize', 10, 'LineWidth', 2);
    % plot(centroid(1), centroid(2), 'co', 'MarkerSize', 10, 'LineWidth', 2);
    % text(centroid(1) + 10, centroid(2), 'Centroid', 'FontName', "TsangerYuMo W03", 'Color', 'b', 'FontSize', 14);
    plot(farCol, farRow, '*', 'Color', [0.25 0.25 0.25], 'MarkerSize', 18, 'LineWidth', 2);
    % text(farCol + 10, farRow, 'Farthest', 'FontName', "TsangerYuMo W03", 'Color', 'g', 'FontSize', 14);
    plot(nearCol, nearRow, '*', 'Color', [1 1 1], 'MarkerSize', 18, 'LineWidth', 2);
    % text(nearCol + 10, nearRow, 'Nearest', 'FontName', "TsangerYuMo W03", 'Color', 'y', 'FontSize', 14);
    % title('Corner Representation',"FontSize",16,"FontName","TsangerYuMo W03","FontWeight","bold");
    hold off;
    pause(0.1);
    exportgraphics(gcf,strcat(Storing_Path,'\Corner_SimH_DTM_',Current_Image_Name),'Resolution', 800);
    
    % Determine image dimensions.
    [nrow, ncol] = size(U);  
    
    %% ACM-Based Micro-Doppler Signature Extraction
    % Initialize the level set functions.
    ic = nrow / 2;  % Center row index.
    jc = ncol / 2;  % Center column index.
    phi = initial_sdf2circle(nrow, ncol, ic, jc, r, fun_n, far_pixel, near_pixel);  % Initialize two level set functions.
    
    % Prepare image data for computations.
    I = double(U);  
    I_max = max(max(I));
    I_min = min(min(I));
    I_Normalized = (I-I_min)/(I_max-I_min);
    
    % Visualize the initial level set functions.
    figure(2);       
    set(gcf, 'Position', [500, 450, 400, 500]);  % Set window size to 4:5 ratio.
    imagesc(1-I_Normalized);
    
    % Set up views and labels of ACM numerical solving visualization.
    colormap(CList);
    xticks([0 Estimation_Resolution/4 Estimation_Resolution/2 Estimation_Resolution*3/4 Estimation_Resolution]);
    xticklabels({'','','','',''});
    yticks([0 Estimation_Resolution/4 Estimation_Resolution/2 Estimation_Resolution*3/4 Estimation_Resolution]);
    yticklabels({'','','','',''});
    set(gca,'FontName','TsangerYuMo W03','FontSize',14,'Fontweight','normal');
    % xlabel('Time (s)','FontSize',15,'FontName',"TsangerYuMo W03","FontWeight","normal");
    % ylabel('Doppler (Hz)','FontSize',15,'FontName',"TsangerYuMo W03","FontWeight","normal");
    % title('RTM for ACM Solving',"FontSize",16,"FontName","TsangerYuMo W03","FontWeight","bold"); 

    % Contour visualization of the level set phi_1 input to the function.
    [Contour_Edge,Contour_Obj] = contour(phi(:,:,1),[0 0],'Color',[1 1 1],'LineWidth',2);
    % Set up views and labels of ACM numerical solving visualization.
    colormap(CList);
    xticks([0 Estimation_Resolution/4 Estimation_Resolution/2 Estimation_Resolution*3/4 Estimation_Resolution]);
    xticklabels({'','','','',''});
    yticks([0 Estimation_Resolution/4 Estimation_Resolution/2 Estimation_Resolution*3/4 Estimation_Resolution]);
    yticklabels({'','','','',''});
    set(gca,'FontName','TsangerYuMo W03','FontSize',14,'Fontweight','normal');
    % xlabel('Time (s)','FontSize',15,'FontName',"TsangerYuMo W03","FontWeight","normal");
    % ylabel('Doppler (Hz)','FontSize',15,'FontName',"TsangerYuMo W03","FontWeight","normal");
    % title('ACM Feature Extraction',"FontSize",16,"FontName","TsangerYuMo W03","FontWeight","bold");

    % Contour visualization of the level set phi_2 input to the function.
    contour(phi(:,:,2),[0 0],'Color',[0.25 0.25 0.25],'LineWidth',2);
    % Set up views and labels of ACM numerical solving visualization.
    colormap(CList);
    xticks([0 Estimation_Resolution/4 Estimation_Resolution/2 Estimation_Resolution*3/4 Estimation_Resolution]);
    xticklabels({'','','','',''});
    yticks([0 Estimation_Resolution/4 Estimation_Resolution/2 Estimation_Resolution*3/4 Estimation_Resolution]);
    yticklabels({'','','','',''});
    set(gca,'FontName','TsangerYuMo W03','FontSize',14,'Fontweight','normal');
    % xlabel('Time (s)','FontSize',15,'FontName',"TsangerYuMo W03","FontWeight","normal");
    % ylabel('Doppler (Hz)','FontSize',15,'FontName',"TsangerYuMo W03","FontWeight","normal");
    % title('ACM Feature Extraction',"FontSize",16,"FontName","TsangerYuMo W03","FontWeight","bold");
    
    % Main loop for evolving the level set functions.
    for k = 1:evolveStep
        phi = EVOLUTION_4PHASE_DR(I, phi, nu, lambda_1, lambda_2, mu, delta_t, epsilon, numIter);  % Update the level set functions.
        if mod(k, 2) == 0  
            pause(.01);     
            imagesc(1-I_Normalized);  
            colormap(gray);  

            hold on;            
            phi_1 = phi(:,:,1); % Extract the first level set function.
            phi_2 = phi(:,:,2); % Extract the second level set function.

            % Contour visualization of the level set phi_1 input to the function.
            [Contour_Edge,Contour_Obj] = contour(phi_1,[0 0],'Color',[1 1 1],'LineWidth',2);
            % Set up views and labels of ACM numerical solving visualization.
            colormap(CList);
            xticks([0 Estimation_Resolution/4 Estimation_Resolution/2 Estimation_Resolution*3/4 Estimation_Resolution]);
            xticklabels({'','','','',''});
            yticks([0 Estimation_Resolution/4 Estimation_Resolution/2 Estimation_Resolution*3/4 Estimation_Resolution]);
            yticklabels({'','','','',''});
            set(gca,'FontName','TsangerYuMo W03','FontSize',14,'Fontweight','normal');
            % xlabel('Time (s)','FontSize',15,'FontName',"TsangerYuMo W03","FontWeight","normal");
            % ylabel('Doppler (Hz)','FontSize',15,'FontName',"TsangerYuMo W03","FontWeight","normal");
            % title('ACM Feature Extraction',"FontSize",16,"FontName","TsangerYuMo W03","FontWeight","bold");
        
            % Contour visualization of the level set phi_2 input to the function.
            contour(phi_2,[0 0],'Color',[0.25 0.25 0.25],'LineWidth',2);
            % Set up views and labels of ACM numerical solving visualization.
            colormap(CList);
            xticks([0 Estimation_Resolution/4 Estimation_Resolution/2 Estimation_Resolution*3/4 Estimation_Resolution]);
            xticklabels({'','','','',''});
            yticks([0 Estimation_Resolution/4 Estimation_Resolution/2 Estimation_Resolution*3/4 Estimation_Resolution]);
            yticklabels({'','','','',''});
            set(gca,'FontName','TsangerYuMo W03','FontSize',14,'Fontweight','normal');
            % xlabel('Time (s)','FontSize',15,'FontName',"TsangerYuMo W03","FontWeight","normal");
            % ylabel('Doppler (Hz)','FontSize',15,'FontName',"TsangerYuMo W03","FontWeight","normal");
            % title('ACM Feature Extraction',"FontSize",16,"FontName","TsangerYuMo W03","FontWeight","bold");

            hold off;          
        end
    end
    pause(0.1);
    exportgraphics(gcf,strcat(Storing_Path,'\ACM_SimH_DTM_',Current_Image_Name),'Resolution', 800);
    
    %% Visualization of Feature Extraction
    % Normalization of two level sets.
    phi_1_max = max(max(phi_1));
    phi_1_min = min(min(phi_1));
    phi_1_Normalized = (phi_1-phi_1_min)/(phi_1_max-phi_1_min);
    phi_2_max = max(max(phi_2));
    phi_2_min = min(min(phi_2));
    phi_2_Normalized = (phi_2-phi_2_min)/(phi_2_max-phi_2_min);
    
    % Visualize the final level set functions as 3D meshes.
    figure(3);      
    set(gcf, 'Position', [950, 450, 400, 500]);  % Set window size to 4:5 ratio.
    mesh(flip(phi_1_Normalized));             
    
    % Set up views and labels of phi_1 visualization.
    colormap(CList_Flip);
    view(45,25); % Set to 3D view, 45 degree viewing azimuth angle, 45 degree viewing pitch angle.
    zticks([0 1]);
    set(gca,'zticklabel',[]);
    zticklabels({'',''});
    xticks([0 Estimation_Resolution/4 Estimation_Resolution/2 Estimation_Resolution*3/4 Estimation_Resolution]);
    xticklabels({'','','','',''});
    yticks([0 Estimation_Resolution/4 Estimation_Resolution/2 Estimation_Resolution*3/4 Estimation_Resolution]);
    yticklabels({'','','','',''});
    set(gca,'FontName','TsangerYuMo W03','FontSize',14,'Fontweight','normal');
    % xlabel('Time (s)','FontSize',15,'FontName',"TsangerYuMo W03","FontWeight","normal");
    % ylabel('Doppler (Hz)','FontSize',15,'FontName',"TsangerYuMo W03","FontWeight","normal");
    % zlabel('Amp','FontSize',15,'FontName',"TsangerYuMo W03","FontWeight","normal");
    % title('Micro-Doppler Signature',"FontSize",16,"FontName","TsangerYuMo W03","FontWeight","bold");
    pause(0.1);
    exportgraphics(gcf,strcat(Storing_Path,'\mD_SimH_DTM_',Current_Image_Name),'Resolution', 800);
    
    figure(4);   
    set(gcf, 'Position', [1400, 450, 400, 500]);  % Set window size to 4:5 ratio.
    mesh(flip(phi_2_Normalized)); 
    
    % Set up views and labels of phi_2 visualization.
    colormap(CList_Flip);
    view(45,25); % Set to 3D view, 45 degree viewing azimuth angle, 45 degree viewing pitch angle.
    zticks([0 1]);
    set(gca,'zticklabel',[]);
    zticklabels({'',''});
    xticks([0 Estimation_Resolution/4 Estimation_Resolution/2 Estimation_Resolution*3/4 Estimation_Resolution]);
    xticklabels({'','','','',''});
    yticks([0 Estimation_Resolution/4 Estimation_Resolution/2 Estimation_Resolution*3/4 Estimation_Resolution]);
    yticklabels({'','','','',''});
    set(gca,'FontName','TsangerYuMo W03','FontSize',14,'Fontweight','normal');
    % xlabel('Time (s)','FontSize',15,'FontName',"TsangerYuMo W03","FontWeight","normal");
    % ylabel('Doppler (Hz)','FontSize',15,'FontName',"TsangerYuMo W03","FontWeight","normal");
    % zlabel('Amp','FontSize',15,'FontName',"TsangerYuMo W03","FontWeight","normal");
    % title('Noise Background',"FontSize",16,"FontName","TsangerYuMo W03","FontWeight","bold");
    
    %% PointCloud Generation
    % PointCloud = Select_Points_for_Columns(phi_1_Normalized,1);
    PointCloud = Contour_Edge;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generating the Visualizations of Measured RTM Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for Image_Name = 1:length(Class_Names)
    Current_Image_Name = strcat(Class_Names(Image_Name),'.png'); % Define the name of the current image.
    Data_Path = fullfile('Visualizations','RW_RTM_Images',Current_Image_Name); % Define the path of the current image.
    Storing_Path = fullfile('Visualizations','RW_RTM_12_Activities');

    %% Re-Initialization of Matlab Script
    % clear all;
    close all;
    % clc;
    % CList = load('JoeyBG_CList.mat').CList; % My favorite colormap preset.
    % CList_Flip = flipud(CList);
    % Class_Names = load("Class_Names.mat").Class_Names; 

    %% Parameter Definition
    delta_t = .1;                                   % Time step for evolution of level set functions.
    lambda_1 = 1;                                   % Weight for data fitting term.
    lambda_2 = 1;                                   % Weight for data fitting term.
    h = 1;                                          % Spatial step size.
    epsilon = 1;                                    % Parameter for Heaviside function approximation
    nu = .001 * 255 * 255;                          % Weight for curvature term, scaled by image intensity range.
    fun_n = 2;                                      % Number of level sets (**Fixed in this script with value 2 for four-phase segmentation**).
    mu = 0.5;                                       % Coefficient for the distance regularizing term.
    numIter = 30;                                   % The number of iterations per evolution step.
    evolveStep = 70;                                % Evolve the level set functions over 'evolveStep' steps.
    r = 32;                                         % Radius of initial circles for level set functions.
    Corner_Threshold_Ratio = 0.3;                   % Threshold cutting ratio of image preprocessing in corner representation.
    Estimation_Resolution = 256;                    % Resize scale of images in this work.
    Max_Range = 4;                                  % Max range bin of RTM, unit: m.
    Max_Doppler = 20;                               % Max Doppler bin of DTM, unit: Hz.
    Max_Time = 4;                                   % Max slow time of both RTM and DTM, unit: s.
    
    %% Readin the Datas for Recognition
    % Load and preprocess the input image.
    Img = imread(Data_Path);  
    U = double(Img(:,:,2));      
    U = imresize(U, [Estimation_Resolution Estimation_Resolution]); 
    % Detect corner points for nearest / farest initialization.
    [far_pixel, near_pixel, II_Normalized, locations, farRow, farCol, nearRow, nearCol] = Corner_Representation(U,Corner_Threshold_Ratio);

    % Plot the results of corner representation.
    figure(1);
    set(gcf, 'Position', [50, 450, 400, 500]);  % Set window size to 4:5 ratio.
    imagesc(1-II_Normalized);
    
    % Set up views and labels of ACM numerical solving visualization.
    colormap(CList_Flip);
    xticks([0 Estimation_Resolution/4 Estimation_Resolution/2 Estimation_Resolution*3/4 Estimation_Resolution]);
    xticklabels({'','','','',''});
    yticks([0 Estimation_Resolution/4 Estimation_Resolution/2 Estimation_Resolution*3/4 Estimation_Resolution]);
    yticklabels({'','','','',''});
    set(gca,'FontName','TsangerYuMo W03','FontSize',14,'Fontweight','normal');
    % xlabel('Time (s)','FontSize',15,'FontName',"TsangerYuMo W03","FontWeight","normal");
    % ylabel('Range (m)','FontSize',15,'FontName',"TsangerYuMo W03","FontWeight","normal");
    % title('Corner Representation',"FontSize",16,"FontName","TsangerYuMo W03","FontWeight","bold");    
    hold on;
    plot(locations(:,1), locations(:,2), 'o', 'Color', '#BEE9FF', 'MarkerSize', 10, 'LineWidth', 2);
    % plot(centroid(1), centroid(2), 'co', 'MarkerSize', 10, 'LineWidth', 2);
    % text(centroid(1) + 10, centroid(2), 'Centroid', 'FontName', "TsangerYuMo W03", 'Color', 'b', 'FontSize', 14);
    plot(farCol, farRow, '*', 'Color', [0.25 0.25 0.25], 'MarkerSize', 18, 'LineWidth', 2);
    % text(farCol + 10, farRow, 'Farthest', 'FontName', "TsangerYuMo W03", 'Color', 'g', 'FontSize', 14);
    plot(nearCol, nearRow, '*', 'Color', [1 1 1], 'MarkerSize', 18, 'LineWidth', 2);
    % text(nearCol + 10, nearRow, 'Nearest', 'FontName', "TsangerYuMo W03", 'Color', 'y', 'FontSize', 14);
    % title('Corner Representation',"FontSize",16,"FontName","TsangerYuMo W03","FontWeight","bold");
    hold off;
    pause(0.1);
    exportgraphics(gcf,strcat(Storing_Path,'\Corner_RW_RTM_',Current_Image_Name),'Resolution', 800);
    
    % Determine image dimensions.
    [nrow, ncol] = size(U);  
    
    %% ACM-Based Micro-Doppler Signature Extraction
    % Initialize the level set functions.
    ic = nrow / 2;  % Center row index.
    jc = ncol / 2;  % Center column index.
    phi = initial_sdf2circle(nrow, ncol, ic, jc, r, fun_n, far_pixel, near_pixel);  % Initialize two level set functions.
    
    % Prepare image data for computations.
    I = double(U);  
    I_max = max(max(I));
    I_min = min(min(I));
    I_Normalized = (I-I_min)/(I_max-I_min);
    
    % Visualize the initial level set functions.
    figure(2);       
    set(gcf, 'Position', [500, 450, 400, 500]);  % Set window size to 4:5 ratio.
    imagesc(1-I_Normalized);
    
    % Set up views and labels of ACM numerical solving visualization.
    colormap(CList_Flip);
    xticks([0 Estimation_Resolution/4 Estimation_Resolution/2 Estimation_Resolution*3/4 Estimation_Resolution]);
    xticklabels({'','','','',''});
    yticks([0 Estimation_Resolution/4 Estimation_Resolution/2 Estimation_Resolution*3/4 Estimation_Resolution]);
    yticklabels({'','','','',''});
    set(gca,'FontName','TsangerYuMo W03','FontSize',14,'Fontweight','normal');
    % xlabel('Time (s)','FontSize',15,'FontName',"TsangerYuMo W03","FontWeight","normal");
    % ylabel('Range (m)','FontSize',15,'FontName',"TsangerYuMo W03","FontWeight","normal");
    % title('RTM for ACM Solving',"FontSize",16,"FontName","TsangerYuMo W03","FontWeight","bold"); 

    % Contour visualization of the level set phi_1 input to the function.
    [Contour_Edge,Contour_Obj] = contour(phi(:,:,1),[0 0],'Color',[1 1 1],'LineWidth',2);
    % Set up views and labels of ACM numerical solving visualization.
    colormap(CList_Flip);
    xticks([0 Estimation_Resolution/4 Estimation_Resolution/2 Estimation_Resolution*3/4 Estimation_Resolution]);
    xticklabels({'','','','',''});
    yticks([0 Estimation_Resolution/4 Estimation_Resolution/2 Estimation_Resolution*3/4 Estimation_Resolution]);
    yticklabels({'','','','',''});
    set(gca,'FontName','TsangerYuMo W03','FontSize',14,'Fontweight','normal');
    % xlabel('Time (s)','FontSize',15,'FontName',"TsangerYuMo W03","FontWeight","normal");
    % ylabel('Range (m)','FontSize',15,'FontName',"TsangerYuMo W03","FontWeight","normal");
    % title('ACM Feature Extraction',"FontSize",16,"FontName","TsangerYuMo W03","FontWeight","bold");

    % Contour visualization of the level set phi_2 input to the function.
    contour(phi(:,:,2),[0 0],'Color',[0.25 0.25 0.25],'LineWidth',2);
    % Set up views and labels of ACM numerical solving visualization.
    colormap(CList_Flip);
    xticks([0 Estimation_Resolution/4 Estimation_Resolution/2 Estimation_Resolution*3/4 Estimation_Resolution]);
    xticklabels({'','','','',''});
    yticks([0 Estimation_Resolution/4 Estimation_Resolution/2 Estimation_Resolution*3/4 Estimation_Resolution]);
    yticklabels({'','','','',''});
    set(gca,'FontName','TsangerYuMo W03','FontSize',14,'Fontweight','normal');
    % xlabel('Time (s)','FontSize',15,'FontName',"TsangerYuMo W03","FontWeight","normal");
    % ylabel('Range (m)','FontSize',15,'FontName',"TsangerYuMo W03","FontWeight","normal");
    % title('ACM Feature Extraction',"FontSize",16,"FontName","TsangerYuMo W03","FontWeight","bold");
    
    % Main loop for evolving the level set functions.
    for k = 1:evolveStep
        phi = EVOLUTION_4PHASE_DR(I, phi, nu, lambda_1, lambda_2, mu, delta_t, epsilon, numIter);  % Update the level set functions.
        if mod(k, 2) == 0  
            pause(.01);     
            imagesc(1-I_Normalized);  
            colormap(gray);  

            hold on;            
            phi_1 = phi(:,:,1); % Extract the first level set function.
            phi_2 = phi(:,:,2); % Extract the second level set function.

            % Contour visualization of the level set phi_1 input to the function.
            [Contour_Edge,Contour_Obj] = contour(phi_1,[0 0],'Color',[1 1 1],'LineWidth',2);
            % Set up views and labels of ACM numerical solving visualization.
            colormap(CList_Flip);
            xticks([0 Estimation_Resolution/4 Estimation_Resolution/2 Estimation_Resolution*3/4 Estimation_Resolution]);
            xticklabels({'','','','',''});
            yticks([0 Estimation_Resolution/4 Estimation_Resolution/2 Estimation_Resolution*3/4 Estimation_Resolution]);
            yticklabels({'','','','',''});
            set(gca,'FontName','TsangerYuMo W03','FontSize',14,'Fontweight','normal');
            % xlabel('Time (s)','FontSize',15,'FontName',"TsangerYuMo W03","FontWeight","normal");
            % ylabel('Range (m)','FontSize',15,'FontName',"TsangerYuMo W03","FontWeight","normal");
            % title('ACM Feature Extraction',"FontSize",16,"FontName","TsangerYuMo W03","FontWeight","bold");
        
            % Contour visualization of the level set phi_2 input to the function.
            contour(phi_2,[0 0],'Color',[0.25 0.25 0.25],'LineWidth',2);
            % Set up views and labels of ACM numerical solving visualization.
            colormap(CList_Flip);
            xticks([0 Estimation_Resolution/4 Estimation_Resolution/2 Estimation_Resolution*3/4 Estimation_Resolution]);
            xticklabels({'','','','',''});
            yticks([0 Estimation_Resolution/4 Estimation_Resolution/2 Estimation_Resolution*3/4 Estimation_Resolution]);
            yticklabels({'','','','',''});
            set(gca,'FontName','TsangerYuMo W03','FontSize',14,'Fontweight','normal');
            % xlabel('Time (s)','FontSize',15,'FontName',"TsangerYuMo W03","FontWeight","normal");
            % ylabel('Range (m)','FontSize',15,'FontName',"TsangerYuMo W03","FontWeight","normal");
            % title('ACM Feature Extraction',"FontSize",16,"FontName","TsangerYuMo W03","FontWeight","bold");

            hold off;          
        end
    end
    pause(0.1);
    exportgraphics(gcf,strcat(Storing_Path,'\ACM_RW_RTM_',Current_Image_Name),'Resolution', 800);
    
    %% Visualization of Feature Extraction
    % Normalization of two level sets.
    phi_1_max = max(max(phi_1));
    phi_1_min = min(min(phi_1));
    phi_1_Normalized = (phi_1-phi_1_min)/(phi_1_max-phi_1_min);
    phi_2_max = max(max(phi_2));
    phi_2_min = min(min(phi_2));
    phi_2_Normalized = (phi_2-phi_2_min)/(phi_2_max-phi_2_min);
    
    % Visualize the final level set functions as 3D meshes.
    figure(3);      
    set(gcf, 'Position', [950, 450, 400, 500]);  % Set window size to 4:5 ratio.
    mesh(flip(phi_1_Normalized));             
    
    % Set up views and labels of phi_1 visualization.
    colormap(CList_Flip);
    view(45,25); % Set to 3D view, 45 degree viewing azimuth angle, 45 degree viewing pitch angle.
    zticks([0 1]);
    set(gca,'zticklabel',[]);
    zticklabels({'',''});
    xticks([0 Estimation_Resolution/4 Estimation_Resolution/2 Estimation_Resolution*3/4 Estimation_Resolution]);
    xticklabels({'','','','',''});
    yticks([0 Estimation_Resolution/4 Estimation_Resolution/2 Estimation_Resolution*3/4 Estimation_Resolution]);
    yticklabels({'','','','',''});
    set(gca,'FontName','TsangerYuMo W03','FontSize',14,'Fontweight','normal');
    % xlabel('Time (s)','FontSize',15,'FontName',"TsangerYuMo W03","FontWeight","normal");
    % ylabel('Range (m)','FontSize',15,'FontName',"TsangerYuMo W03","FontWeight","normal");
    % zlabel('Amp','FontSize',15,'FontName',"TsangerYuMo W03","FontWeight","normal");
    % title('Micro-Doppler Signature',"FontSize",16,"FontName","TsangerYuMo W03","FontWeight","bold");
    pause(0.1);
    exportgraphics(gcf,strcat(Storing_Path,'\mD_RW_RTM_',Current_Image_Name),'Resolution', 800);
    
    figure(4);   
    set(gcf, 'Position', [1400, 450, 400, 500]);  % Set window size to 4:5 ratio.
    mesh(flip(phi_2_Normalized)); 
    
    % Set up views and labels of phi_2 visualization.
    colormap(CList_Flip);
    view(45,25); % Set to 3D view, 45 degree viewing azimuth angle, 45 degree viewing pitch angle.
    zticks([0 1]);
    set(gca,'zticklabel',[]);
    zticklabels({'',''});
    xticks([0 Estimation_Resolution/4 Estimation_Resolution/2 Estimation_Resolution*3/4 Estimation_Resolution]);
    xticklabels({'','','','',''});
    yticks([0 Estimation_Resolution/4 Estimation_Resolution/2 Estimation_Resolution*3/4 Estimation_Resolution]);
    yticklabels({'','','','',''});
    set(gca,'FontName','TsangerYuMo W03','FontSize',14,'Fontweight','normal');
    % xlabel('Time (s)','FontSize',15,'FontName',"TsangerYuMo W03","FontWeight","normal");
    % ylabel('Range (m)','FontSize',15,'FontName',"TsangerYuMo W03","FontWeight","normal");
    % zlabel('Amp','FontSize',15,'FontName',"TsangerYuMo W03","FontWeight","normal");
    % title('Noise Background',"FontSize",16,"FontName","TsangerYuMo W03","FontWeight","bold");
    
    %% PointCloud Generation
    % PointCloud = Select_Points_for_Columns(phi_1_Normalized,1);
    PointCloud = Contour_Edge;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generating the Visualizations of Measured DTM Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for Image_Name = 1:length(Class_Names)
    Current_Image_Name = strcat(Class_Names(Image_Name),'.png'); % Define the name of the current image.
    Data_Path = fullfile('Visualizations','RW_DTM_Images',Current_Image_Name); % Define the path of the current image.
    Storing_Path = fullfile('Visualizations','RW_DTM_12_Activities');

    %% Re-Initialization of Matlab Script
    % clear all;
    close all;
    % clc;
    % CList = load('JoeyBG_CList.mat').CList; % My favorite colormap preset.
    % CList_Flip = flipud(CList);
    % Class_Names = load("Class_Names.mat").Class_Names; 

    %% Parameter Definition
    delta_t = .1;                                   % Time step for evolution of level set functions.
    lambda_1 = 1;                                   % Weight for data fitting term.
    lambda_2 = 1;                                   % Weight for data fitting term.
    h = 1;                                          % Spatial step size.
    epsilon = 1;                                    % Parameter for Heaviside function approximation
    nu = .001 * 255 * 255;                          % Weight for curvature term, scaled by image intensity range.
    fun_n = 2;                                      % Number of level sets (**Fixed in this script with value 2 for four-phase segmentation**).
    mu = 0.5;                                       % Coefficient for the distance regularizing term.
    numIter = 50;                                   % The number of iterations per evolution step.
    evolveStep = 70;                                % Evolve the level set functions over 'evolveStep' steps.
    r = 32;                                         % Radius of initial circles for level set functions.
    Corner_Threshold_Ratio = 0.1;                   % Threshold cutting ratio of image preprocessing in corner representation.
    Estimation_Resolution = 256;                    % Resize scale of images in this work.
    Max_Range = 4;                                  % Max range bin of RTM, unit: m.
    Max_Doppler = 20;                               % Max Doppler bin of DTM, unit: Hz.
    Max_Time = 4;                                   % Max slow time of both RTM and DTM, unit: s.
    
    %% Readin the Datas for Recognition
    % Load and preprocess the input image.
    Img = imread(Data_Path);  
    U = double(Img(:,:,2));      
    U = imresize(U, [Estimation_Resolution Estimation_Resolution]); 
    % Detect corner points for nearest / farest initialization.
    [far_pixel, near_pixel, II_Normalized, locations, farRow, farCol, nearRow, nearCol] = Corner_Representation(U,Corner_Threshold_Ratio);

    % Plot the results of corner representation.
    figure(1);
    set(gcf, 'Position', [50, 450, 400, 500]);  % Set window size to 4:5 ratio.
    imagesc(1-II_Normalized);
    
    % Set up views and labels of ACM numerical solving visualization.
    colormap(CList_Flip);
    xticks([0 Estimation_Resolution/4 Estimation_Resolution/2 Estimation_Resolution*3/4 Estimation_Resolution]);
    xticklabels({'','','','',''});
    yticks([0 Estimation_Resolution/4 Estimation_Resolution/2 Estimation_Resolution*3/4 Estimation_Resolution]);
    yticklabels({'','','','',''});
    set(gca,'FontName','TsangerYuMo W03','FontSize',14,'Fontweight','normal');
    % xlabel('Time (s)','FontSize',15,'FontName',"TsangerYuMo W03","FontWeight","normal");
    % ylabel('Doppler (Hz)','FontSize',15,'FontName',"TsangerYuMo W03","FontWeight","normal");
    % title('Corner Representation',"FontSize",16,"FontName","TsangerYuMo W03","FontWeight","bold");    
    hold on;
    plot(locations(:,1), locations(:,2), 'o', 'Color', '#BEE9FF', 'MarkerSize', 10, 'LineWidth', 2);
    % plot(centroid(1), centroid(2), 'co', 'MarkerSize', 10, 'LineWidth', 2);
    % text(centroid(1) + 10, centroid(2), 'Centroid', 'FontName', "TsangerYuMo W03", 'Color', 'b', 'FontSize', 14);
    plot(farCol, farRow, '*', 'Color', [0.25 0.25 0.25], 'MarkerSize', 18, 'LineWidth', 2);
    % text(farCol + 10, farRow, 'Farthest', 'FontName', "TsangerYuMo W03", 'Color', 'g', 'FontSize', 14);
    plot(nearCol, nearRow, '*', 'Color', [1 1 1], 'MarkerSize', 18, 'LineWidth', 2);
    % text(nearCol + 10, nearRow, 'Nearest', 'FontName', "TsangerYuMo W03", 'Color', 'y', 'FontSize', 14);
    % title('Corner Representation',"FontSize",16,"FontName","TsangerYuMo W03","FontWeight","bold");
    hold off;
    pause(0.1);
    exportgraphics(gcf,strcat(Storing_Path,'\Corner_RW_DTM_',Current_Image_Name),'Resolution', 800);
    
    % Determine image dimensions.
    [nrow, ncol] = size(U);  
    
    %% ACM-Based Micro-Doppler Signature Extraction
    % Initialize the level set functions.
    ic = nrow / 2;  % Center row index.
    jc = ncol / 2;  % Center column index.
    phi = initial_sdf2circle(nrow, ncol, ic, jc, r, fun_n, far_pixel, near_pixel);  % Initialize two level set functions.
    
    % Prepare image data for computations.
    I = double(U);  
    I_max = max(max(I));
    I_min = min(min(I));
    I_Normalized = (I-I_min)/(I_max-I_min);
    
    % Visualize the initial level set functions.
    figure(2);       
    set(gcf, 'Position', [500, 450, 400, 500]);  % Set window size to 4:5 ratio.
    imagesc(1-I_Normalized);
    
    % Set up views and labels of ACM numerical solving visualization.
    colormap(CList_Flip);
    xticks([0 Estimation_Resolution/4 Estimation_Resolution/2 Estimation_Resolution*3/4 Estimation_Resolution]);
    xticklabels({'','','','',''});
    yticks([0 Estimation_Resolution/4 Estimation_Resolution/2 Estimation_Resolution*3/4 Estimation_Resolution]);
    yticklabels({'','','','',''});
    set(gca,'FontName','TsangerYuMo W03','FontSize',14,'Fontweight','normal');
    % xlabel('Time (s)','FontSize',15,'FontName',"TsangerYuMo W03","FontWeight","normal");
    % ylabel('Doppler (Hz)','FontSize',15,'FontName',"TsangerYuMo W03","FontWeight","normal");
    % title('RTM for ACM Solving',"FontSize",16,"FontName","TsangerYuMo W03","FontWeight","bold"); 

    % Contour visualization of the level set phi_1 input to the function.
    [Contour_Edge,Contour_Obj] = contour(phi(:,:,1),[0 0],'Color',[1 1 1],'LineWidth',2);
    % Set up views and labels of ACM numerical solving visualization.
    colormap(CList_Flip);
    xticks([0 Estimation_Resolution/4 Estimation_Resolution/2 Estimation_Resolution*3/4 Estimation_Resolution]);
    xticklabels({'','','','',''});
    yticks([0 Estimation_Resolution/4 Estimation_Resolution/2 Estimation_Resolution*3/4 Estimation_Resolution]);
    yticklabels({'','','','',''});
    set(gca,'FontName','TsangerYuMo W03','FontSize',14,'Fontweight','normal');
    % xlabel('Time (s)','FontSize',15,'FontName',"TsangerYuMo W03","FontWeight","normal");
    % ylabel('Doppler (Hz)','FontSize',15,'FontName',"TsangerYuMo W03","FontWeight","normal");
    % title('ACM Feature Extraction',"FontSize",16,"FontName","TsangerYuMo W03","FontWeight","bold");

    % Contour visualization of the level set phi_2 input to the function.
    contour(phi(:,:,2),[0 0],'Color',[0.25 0.25 0.25],'LineWidth',2);
    % Set up views and labels of ACM numerical solving visualization.
    colormap(CList_Flip);
    xticks([0 Estimation_Resolution/4 Estimation_Resolution/2 Estimation_Resolution*3/4 Estimation_Resolution]);
    xticklabels({'','','','',''});
    yticks([0 Estimation_Resolution/4 Estimation_Resolution/2 Estimation_Resolution*3/4 Estimation_Resolution]);
    yticklabels({'','','','',''});
    set(gca,'FontName','TsangerYuMo W03','FontSize',14,'Fontweight','normal');
    % xlabel('Time (s)','FontSize',15,'FontName',"TsangerYuMo W03","FontWeight","normal");
    % ylabel('Doppler (Hz)','FontSize',15,'FontName',"TsangerYuMo W03","FontWeight","normal");
    % title('ACM Feature Extraction',"FontSize",16,"FontName","TsangerYuMo W03","FontWeight","bold");
    
    % Main loop for evolving the level set functions.
    for k = 1:evolveStep
        phi = EVOLUTION_4PHASE_DR(I, phi, nu, lambda_1, lambda_2, mu, delta_t, epsilon, numIter);  % Update the level set functions.
        if mod(k, 2) == 0  
            pause(.01);     
            imagesc(1-I_Normalized);  
            colormap(gray);  

            hold on;            
            phi_1 = phi(:,:,1); % Extract the first level set function.
            phi_2 = phi(:,:,2); % Extract the second level set function.

            % Contour visualization of the level set phi_1 input to the function.
            [Contour_Edge,Contour_Obj] = contour(phi_1,[0 0],'Color',[1 1 1],'LineWidth',2);
            % Set up views and labels of ACM numerical solving visualization.
            colormap(CList_Flip);
            xticks([0 Estimation_Resolution/4 Estimation_Resolution/2 Estimation_Resolution*3/4 Estimation_Resolution]);
            xticklabels({'','','','',''});
            yticks([0 Estimation_Resolution/4 Estimation_Resolution/2 Estimation_Resolution*3/4 Estimation_Resolution]);
            yticklabels({'','','','',''});
            set(gca,'FontName','TsangerYuMo W03','FontSize',14,'Fontweight','normal');
            % xlabel('Time (s)','FontSize',15,'FontName',"TsangerYuMo W03","FontWeight","normal");
            % ylabel('Doppler (Hz)','FontSize',15,'FontName',"TsangerYuMo W03","FontWeight","normal");
            % title('ACM Feature Extraction',"FontSize",16,"FontName","TsangerYuMo W03","FontWeight","bold");
        
            % Contour visualization of the level set phi_2 input to the function.
            contour(phi_2,[0 0],'Color',[0.25 0.25 0.25],'LineWidth',2);
            % Set up views and labels of ACM numerical solving visualization.
            colormap(CList_Flip);
            xticks([0 Estimation_Resolution/4 Estimation_Resolution/2 Estimation_Resolution*3/4 Estimation_Resolution]);
            xticklabels({'','','','',''});
            yticks([0 Estimation_Resolution/4 Estimation_Resolution/2 Estimation_Resolution*3/4 Estimation_Resolution]);
            yticklabels({'','','','',''});
            set(gca,'FontName','TsangerYuMo W03','FontSize',14,'Fontweight','normal');
            % xlabel('Time (s)','FontSize',15,'FontName',"TsangerYuMo W03","FontWeight","normal");
            % ylabel('Doppler (Hz)','FontSize',15,'FontName',"TsangerYuMo W03","FontWeight","normal");
            % title('ACM Feature Extraction',"FontSize",16,"FontName","TsangerYuMo W03","FontWeight","bold");

            hold off;          
        end
    end
    pause(0.1);
    exportgraphics(gcf,strcat(Storing_Path,'\ACM_RW_DTM_',Current_Image_Name),'Resolution', 800);
    
    %% Visualization of Feature Extraction
    % Normalization of two level sets.
    phi_1_max = max(max(phi_1));
    phi_1_min = min(min(phi_1));
    phi_1_Normalized = (phi_1-phi_1_min)/(phi_1_max-phi_1_min);
    phi_2_max = max(max(phi_2));
    phi_2_min = min(min(phi_2));
    phi_2_Normalized = (phi_2-phi_2_min)/(phi_2_max-phi_2_min);
    
    % Visualize the final level set functions as 3D meshes.
    figure(3);      
    set(gcf, 'Position', [950, 450, 400, 500]);  % Set window size to 4:5 ratio.
    mesh(flip(phi_1_Normalized));             
    
    % Set up views and labels of phi_1 visualization.
    colormap(CList_Flip);
    view(45,25); % Set to 3D view, 45 degree viewing azimuth angle, 45 degree viewing pitch angle.
    zticks([0 1]);
    set(gca,'zticklabel',[]);
    zticklabels({'',''});
    xticks([0 Estimation_Resolution/4 Estimation_Resolution/2 Estimation_Resolution*3/4 Estimation_Resolution]);
    xticklabels({'','','','',''});
    yticks([0 Estimation_Resolution/4 Estimation_Resolution/2 Estimation_Resolution*3/4 Estimation_Resolution]);
    yticklabels({'','','','',''});
    set(gca,'FontName','TsangerYuMo W03','FontSize',14,'Fontweight','normal');
    % xlabel('Time (s)','FontSize',15,'FontName',"TsangerYuMo W03","FontWeight","normal");
    % ylabel('Doppler (Hz)','FontSize',15,'FontName',"TsangerYuMo W03","FontWeight","normal");
    % zlabel('Amp','FontSize',15,'FontName',"TsangerYuMo W03","FontWeight","normal");
    % title('Micro-Doppler Signature',"FontSize",16,"FontName","TsangerYuMo W03","FontWeight","bold");
    pause(0.1);
    exportgraphics(gcf,strcat(Storing_Path,'\mD_RW_DTM_',Current_Image_Name),'Resolution', 800);
    
    figure(4);   
    set(gcf, 'Position', [1400, 450, 400, 500]);  % Set window size to 4:5 ratio.
    mesh(flip(phi_2_Normalized)); 
    
    % Set up views and labels of phi_2 visualization.
    colormap(CList_Flip);
    view(45,25); % Set to 3D view, 45 degree viewing azimuth angle, 45 degree viewing pitch angle.
    zticks([0 1]);
    set(gca,'zticklabel',[]);
    zticklabels({'',''});
    xticks([0 Estimation_Resolution/4 Estimation_Resolution/2 Estimation_Resolution*3/4 Estimation_Resolution]);
    xticklabels({'','','','',''});
    yticks([0 Estimation_Resolution/4 Estimation_Resolution/2 Estimation_Resolution*3/4 Estimation_Resolution]);
    yticklabels({'','','','',''});
    set(gca,'FontName','TsangerYuMo W03','FontSize',14,'Fontweight','normal');
    % xlabel('Time (s)','FontSize',15,'FontName',"TsangerYuMo W03","FontWeight","normal");
    % ylabel('Doppler (Hz)','FontSize',15,'FontName',"TsangerYuMo W03","FontWeight","normal");
    % zlabel('Amp','FontSize',15,'FontName',"TsangerYuMo W03","FontWeight","normal");
    % title('Noise Background',"FontSize",16,"FontName","TsangerYuMo W03","FontWeight","bold");
    
    %% PointCloud Generation
    % PointCloud = Select_Points_for_Columns(phi_1_Normalized,1);
    PointCloud = Contour_Edge;
end