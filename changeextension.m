clear
close all
clc

addpath(genpath('J:\Data\Matlab'));             
%%

% Select the directory to search
directory = uigetdir;
% List all .swc files
fileList = dir([directory, '\*.swc']); 
% Loop through each .swc file, copy it and give new extension: .txt
for i = 1:numel(fileList)
    file = fullfile(directory, fileList(i).name);
    [tempDir, tempFile] = fileparts(file); 
    status = copyfile(file, fullfile(tempDir, [tempFile, '.txt']));
    % Delete the .out file; but don't do this until you've backed up your data!!
    % delete(file)  
end