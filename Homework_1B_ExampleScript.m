% Homework_1B_ExampleScript.m
% SD2709 Underwater Technology
% Homework: Module 1 - Part B
% Ocean Properties Analysis Using MATLAB
% Author: Axel Magnusson
% Date: 1/9/2025

% Description:
% This script demonstrates how to analyze ocean properties using MATLAB.
% It allows the user to choose between plotting real Argo float data or synthetic data.

% Clear the workspace and command window
clear;
clc;

% User choice for data type
choice = menu('Choose the type of data to plot:', 'Real Argo Data', 'Synthetic Data');

switch choice
    case 1
        % Plot real Argo float data
        argoFilename = 'PR_PF_4903884.csv';  % Replace with your Argo CSV file name
        plotArgoFloatData(argoFilename);
        
    case 2
        % Plot synthetic data
        plotSyntheticData();
        
    otherwise
        disp('No valid choice made. Exiting...');
        return;
end
