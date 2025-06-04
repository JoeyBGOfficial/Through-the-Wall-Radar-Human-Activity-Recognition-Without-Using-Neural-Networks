%% Function for Feature Extraction and Pointcloud Generation (Measured, RTM)
% Author: JoeyBG.
% Improved By: JoeyBG.
% Time: 2025/5/29.
% Affiliation: Beijing Institute of Technology.
%
% Information:
% This script first shows through-the-wall radar human motion echo model and data processing methods to 
% generate both range-time map (RTM) and Doppler-time map (DTM).
% Then, this script also implements a four-phase level set method to segment an image into four regions
% using two level set functions. It incorporates distance regularization to maintain the
% smoothness of the level set functions during evolution. The script loads an image, initializes
% the level set functions based on detected features, evolves them over iterations, and visualizes
% both the intermediate and final segmentation results.
% Finally, this script generates model-based sparse point cloud of the contour of ACM feature extraction.
%
% Some Key Parameters:
%   - delta_t: Time step for the evolution of level set functions.
%   - lambda_1, lambda_2: Weights for the data fitting terms (inside and outside regions).
%   - h: Spatial step size (not utilized in this implementation).
%   - epsilon: Parameter for the Heaviside function approximation.
%   - nu: Weight for the curvature regularization term.
%   - fun_n: Number of level set functions (2, enabling four-phase segmentation).
%   - numIter: Number of iterations per evolution step.
%   - mu: Coefficient for the distance regularizing term.
%   - numIter: Maximum iterations of ACM estimation.
%   - evolveStep: The steps of evolution of levet set functions.
%   - r: The radius of level set initialization based on two recorded coordinates.
%   - Corner_Threshold_Ratio: Threshold cutting ratio of image preprocessing in corner representation.
%   - Estimation_Resolution: The resize scale of input image.
%   - Max_Range: Maximum range bin in RTM.
%   - Max_Doppler: Maximum positive Doppler bin in DTM.
%   - Max_Time: Maximum slow time index in both RTM and DTM.
%
% The script processes the images in 'Example_Datas/' as examples, and uses corner detection to guide
% the initialization of the level set functions. It visualizes the evolution every two steps and
% displays the final results as 3D meshes.
%
% Dependencies:
%   - initial_sdf2circle.m: Initializes the level set functions as signed distance fields.
%   - Corner_Representation.m: Detects corner points and computes pixel positions.
%   - EVOLUTION_4PHASE_DR.m: Evolves the level set functions with distance regularization.
%
% Note: Ensure the pre-defined image path is in the 'Example_Datas/' directory before running the script as an example.

%% Function Body
function PointCloud = Feature_Extraction_RWRTM(Data_Path)
    %% Initialization of Matlab Script
    % clear all;
    close all;
    % clc;
    CList = load('JoeyBG_CList.mat').CList; % My favorite colormap preset.
    CList_Flip = flipud(CList);
    Class_Names = load("Class_Names.mat").Class_Names; 

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
    xticklabels({'0',num2str(Max_Time/4),num2str(Max_Time/2),num2str(Max_Time*3/4),num2str(Max_Time)});
    yticks([0 Estimation_Resolution/4 Estimation_Resolution/2 Estimation_Resolution*3/4 Estimation_Resolution]);
    yticklabels({num2str(Max_Range),num2str(Max_Range*3/4),num2str(Max_Range/2),num2str(Max_Range/4),'0'});
    set(gca,'FontName','TsangerYuMo W03','FontSize',14,'Fontweight','normal');
    xlabel('Time (s)','FontSize',15,'FontName',"TsangerYuMo W03","FontWeight","normal");
    ylabel('Range (m)','FontSize',15,'FontName',"TsangerYuMo W03","FontWeight","normal");
    title('Corner Representation',"FontSize",16,"FontName","TsangerYuMo W03","FontWeight","bold");    
    hold on;
    plot(locations(:,1), locations(:,2), 'o', 'Color', '#BEE9FF', 'MarkerSize', 10, 'LineWidth', 2);
    % plot(centroid(1), centroid(2), 'co', 'MarkerSize', 10, 'LineWidth', 2);
    % text(centroid(1) + 10, centroid(2), 'Centroid', 'FontName', "TsangerYuMo W03", 'Color', 'b', 'FontSize', 14);
    plot(farCol, farRow, '*', 'Color', [0.25 0.25 0.25], 'MarkerSize', 18, 'LineWidth', 2);
    % text(farCol + 10, farRow, 'Farthest', 'FontName', "TsangerYuMo W03", 'Color', 'g', 'FontSize', 14);
    plot(nearCol, nearRow, '*', 'Color', [1 1 1], 'MarkerSize', 18, 'LineWidth', 2);
    % text(nearCol + 10, nearRow, 'Nearest', 'FontName', "TsangerYuMo W03", 'Color', 'y', 'FontSize', 14);
    title('Corner Representation',"FontSize",16,"FontName","TsangerYuMo W03","FontWeight","bold");
    hold off;
    
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
    xticklabels({'0',num2str(Max_Time/4),num2str(Max_Time/2),num2str(Max_Time*3/4),num2str(Max_Time)});
    yticks([0 Estimation_Resolution/4 Estimation_Resolution/2 Estimation_Resolution*3/4 Estimation_Resolution]);
    yticklabels({num2str(Max_Range),num2str(Max_Range*3/4),num2str(Max_Range/2),num2str(Max_Range/4),'0'});
    set(gca,'FontName','TsangerYuMo W03','FontSize',14,'Fontweight','normal');
    xlabel('Time (s)','FontSize',15,'FontName',"TsangerYuMo W03","FontWeight","normal");
    ylabel('Range (m)','FontSize',15,'FontName',"TsangerYuMo W03","FontWeight","normal");
    title('RTM for ACM Solving',"FontSize",16,"FontName","TsangerYuMo W03","FontWeight","bold"); 

    % Contour visualization of the level set phi_1 input to the function.
    [Contour_Edge,Contour_Obj] = contour(phi(:,:,1),[0 0],'Color',[1 1 1],'LineWidth',2);
    % Set up views and labels of ACM numerical solving visualization.
    colormap(CList_Flip);
    xticks([0 Estimation_Resolution/4 Estimation_Resolution/2 Estimation_Resolution*3/4 Estimation_Resolution]);
    xticklabels({'0',num2str(Max_Time/4),num2str(Max_Time/2),num2str(Max_Time*3/4),num2str(Max_Time)});
    yticks([0 Estimation_Resolution/4 Estimation_Resolution/2 Estimation_Resolution*3/4 Estimation_Resolution]);
    yticklabels({num2str(Max_Range),num2str(Max_Range*3/4),num2str(Max_Range/2),num2str(Max_Range/4),'0'});
    set(gca,'FontName','TsangerYuMo W03','FontSize',14,'Fontweight','normal');
    xlabel('Time (s)','FontSize',15,'FontName',"TsangerYuMo W03","FontWeight","normal");
    ylabel('Range (m)','FontSize',15,'FontName',"TsangerYuMo W03","FontWeight","normal");
    title('ACM Feature Extraction',"FontSize",16,"FontName","TsangerYuMo W03","FontWeight","bold");

    % Contour visualization of the level set phi_2 input to the function.
    contour(phi(:,:,2),[0 0],'Color',[0.25 0.25 0.25],'LineWidth',2);
    % Set up views and labels of ACM numerical solving visualization.
    colormap(CList_Flip);
    xticks([0 Estimation_Resolution/4 Estimation_Resolution/2 Estimation_Resolution*3/4 Estimation_Resolution]);
    xticklabels({'0',num2str(Max_Time/4),num2str(Max_Time/2),num2str(Max_Time*3/4),num2str(Max_Time)});
    yticks([0 Estimation_Resolution/4 Estimation_Resolution/2 Estimation_Resolution*3/4 Estimation_Resolution]);
    yticklabels({num2str(Max_Range),num2str(Max_Range*3/4),num2str(Max_Range/2),num2str(Max_Range/4),'0'});
    set(gca,'FontName','TsangerYuMo W03','FontSize',14,'Fontweight','normal');
    xlabel('Time (s)','FontSize',15,'FontName',"TsangerYuMo W03","FontWeight","normal");
    ylabel('Range (m)','FontSize',15,'FontName',"TsangerYuMo W03","FontWeight","normal");
    title('ACM Feature Extraction',"FontSize",16,"FontName","TsangerYuMo W03","FontWeight","bold");
    
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
            xticklabels({'0',num2str(Max_Time/4),num2str(Max_Time/2),num2str(Max_Time*3/4),num2str(Max_Time)});
            yticks([0 Estimation_Resolution/4 Estimation_Resolution/2 Estimation_Resolution*3/4 Estimation_Resolution]);
            yticklabels({num2str(Max_Range),num2str(Max_Range*3/4),num2str(Max_Range/2),num2str(Max_Range/4),'0'});
            set(gca,'FontName','TsangerYuMo W03','FontSize',14,'Fontweight','normal');
            xlabel('Time (s)','FontSize',15,'FontName',"TsangerYuMo W03","FontWeight","normal");
            ylabel('Range (m)','FontSize',15,'FontName',"TsangerYuMo W03","FontWeight","normal");
            title('ACM Feature Extraction',"FontSize",16,"FontName","TsangerYuMo W03","FontWeight","bold");
        
            % Contour visualization of the level set phi_2 input to the function.
            contour(phi_2,[0 0],'Color',[0.25 0.25 0.25],'LineWidth',2);
            % Set up views and labels of ACM numerical solving visualization.
            colormap(CList_Flip);
            xticks([0 Estimation_Resolution/4 Estimation_Resolution/2 Estimation_Resolution*3/4 Estimation_Resolution]);
            xticklabels({'0',num2str(Max_Time/4),num2str(Max_Time/2),num2str(Max_Time*3/4),num2str(Max_Time)});
            yticks([0 Estimation_Resolution/4 Estimation_Resolution/2 Estimation_Resolution*3/4 Estimation_Resolution]);
            yticklabels({num2str(Max_Range),num2str(Max_Range*3/4),num2str(Max_Range/2),num2str(Max_Range/4),'0'});
            set(gca,'FontName','TsangerYuMo W03','FontSize',14,'Fontweight','normal');
            xlabel('Time (s)','FontSize',15,'FontName',"TsangerYuMo W03","FontWeight","normal");
            ylabel('Range (m)','FontSize',15,'FontName',"TsangerYuMo W03","FontWeight","normal");
            title('ACM Feature Extraction',"FontSize",16,"FontName","TsangerYuMo W03","FontWeight","bold");

            hold off;          
        end
    end
    
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
    zticklabels({'0','1'});
    xticks([0 Estimation_Resolution/4 Estimation_Resolution/2 Estimation_Resolution*3/4 Estimation_Resolution]);
    xticklabels({'0',num2str(Max_Time/4),num2str(Max_Time/2),num2str(Max_Time*3/4),num2str(Max_Time)});
    yticks([0 Estimation_Resolution/4 Estimation_Resolution/2 Estimation_Resolution*3/4 Estimation_Resolution]);
    yticklabels({'0',num2str(Max_Range/4),num2str(Max_Range/2),num2str(Max_Range*3/4),num2str(Max_Range)});
    set(gca,'FontName','TsangerYuMo W03','FontSize',14,'Fontweight','normal');
    xlabel('Time (s)','FontSize',15,'FontName',"TsangerYuMo W03","FontWeight","normal");
    ylabel('Range (m)','FontSize',15,'FontName',"TsangerYuMo W03","FontWeight","normal");
    zlabel('Amp','FontSize',15,'FontName',"TsangerYuMo W03","FontWeight","normal");
    title('Micro-Doppler Signature',"FontSize",16,"FontName","TsangerYuMo W03","FontWeight","bold");
    
    figure(4);   
    set(gcf, 'Position', [1400, 450, 400, 500]);  % Set window size to 4:5 ratio.
    mesh(flip(phi_2_Normalized)); 
    
    % Set up views and labels of phi_2 visualization.
    colormap(CList_Flip);
    view(45,25); % Set to 3D view, 45 degree viewing azimuth angle, 45 degree viewing pitch angle.
    zticks([0 1]);
    set(gca,'zticklabel',[]);
    zticklabels({'0','1'});
    xticks([0 Estimation_Resolution/4 Estimation_Resolution/2 Estimation_Resolution*3/4 Estimation_Resolution]);
    xticklabels({'0',num2str(Max_Time/4),num2str(Max_Time/2),num2str(Max_Time*3/4),num2str(Max_Time)});
    yticks([0 Estimation_Resolution/4 Estimation_Resolution/2 Estimation_Resolution*3/4 Estimation_Resolution]);
    yticklabels({'0',num2str(Max_Range/4),num2str(Max_Range/2),num2str(Max_Range*3/4),num2str(Max_Range)});
    set(gca,'FontName','TsangerYuMo W03','FontSize',14,'Fontweight','normal');
    xlabel('Time (s)','FontSize',15,'FontName',"TsangerYuMo W03","FontWeight","normal");
    ylabel('Range (m)','FontSize',15,'FontName',"TsangerYuMo W03","FontWeight","normal");
    zlabel('Amp','FontSize',15,'FontName',"TsangerYuMo W03","FontWeight","normal");
    title('Noise Background',"FontSize",16,"FontName","TsangerYuMo W03","FontWeight","bold");
    
    %% PointCloud Generation
    % PointCloud = Select_Points_for_Columns(phi_1_Normalized,1);
    PointCloud = Contour_Edge;
end