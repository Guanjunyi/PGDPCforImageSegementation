% Please kindly cite the paper Junyi Guan, Sheng Li, Xiongxiong He, and Jiajia Chen, 
%"Peak-graph-based fast density peak clustering for image segmentation," 
% IEEE SIGNAL PROCESSING LETTERS, 2021,Doi:?
% The code was written by Junyi Guan in 2020.
clear all; close all; clc;
%% loading data
load dataset/jain
%% parameter setting
data = jain(:,1:2);
k = 10;
%% PGDPC clustering
PGDPC(data,k);