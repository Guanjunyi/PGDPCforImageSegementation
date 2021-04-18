% Please kindly cite the paper Junyi Guan, Sheng Li, Xiongxiong He, and Jiajia Chen,
%"Peak-graph-based fast density peak clustering for image segmentation,"
% IEEE SIGNAL PROCESSING LETTERS, 2021,Doi:10.1109/LSP.2021.3072794
% The code was written by Junyi Guan in 2020.
clear all; close all; clc;
%% reading image data
imageName = '3063';
% imageName = '61034';
% imageName = '100007';
% imageName = '198004';
% imageName = '279005';
imagename=['.../PGDPCforImageSegementation-master/test/images/',imageName,'.jpg'];
fig_ori = imread(imagename);
%% get groundtruth
gtname=['.../PGDPCforImageSegementation-master/test/groundTruth/',imageName,'.mat'];
gtfile = load(gtname);
gt = gtfile.groundTruth;
%% setting number of superpixels
SupN = 500;
%% PGDPC for image Segmentation
[Seg_M,PRI,VOI,GCE,BDE,NC] = PGDPCForSEG(fig_ori,SupN,gt);
