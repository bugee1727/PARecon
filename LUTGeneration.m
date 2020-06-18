%% %%%%%%%%%%%%%%%%%%%%%%%%%%
% Lookup table generation
% Input file: 
%   - "param.mat"
%   - Correction of acquision/system parameters 
% Output file: 
%   - "LUTMat.npy"
% Function required:
%   - "LUTGenerator.mat"
%   - Converting a 2D (raw data) array to a 3D array
% Reference: 
%   - Kim et al. "Deep-learning Image Reconstruction for 
%           Real-time Photoacoustic System." IEEE TMI, 2020
% Minwoo Kim, mkim180@uw.edu, 05/18/2020

clear all
load('table\param.mat'); % load parameters

outMat = LUTFunc(param);
[inx1,inx2,val] = find(outMat);

LUTMat = [inx1,inx2,val];

writeNPY(LUTMat,['table/LUTMat.npy']);
