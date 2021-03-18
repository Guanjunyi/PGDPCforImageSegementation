% Please kindly cite the paper Junyi Guan, Sheng Li, Xiongxiong He, and Jiajia Chen,
%"Peak-graph-based fast density peak clustering for image segmentation,"
% IEEE SIGNAL PROCESSING LETTERS, 2021,Doi:?
% The code was written by Junyi Guan in 2020.
function [Seg_M,PRI,VOI,GCE,BDE,NC] = PGDPCForSEG(fig_ori,supN,gt)
%% parameter setting
k = round(supN/50);
%% gaussian filtering
sigma = 5.0;fig_gf = gausFT(fig_ori,sigma);
%% superpixels generated
[Lp,NSup] = genesups(fig_ori,supN);
%% preproces image data
[PTs Lp] = preproces(fig_gf,Lp,NSup);
%% PGDPC clustering
[CL_sp,~] = PGDPC(PTs,k);
%% NC final number of clusters
NC = max(CL_sp);
%% Assign lable to each pixel
CL = assignlable(CL_sp,Lp);
%% build segmentation label matix
Seg_M = bulidSegLab(fig_ori,CL);
%% evaluation:PRI,VOI
[PRI,VOI,GCE,BDE] = match_segmentations(Seg_M,gt);
%% Label_image
[fs,~,] = Label_image(fig_ori,Seg_M);
%% Label_image of gt
Gt_M = gt{1,5}.Segmentation;
[fs2,~,] = Label_image(fig_ori,Gt_M);
%% show result
figure(3)
subplot(2,2,1); imshow(fig_ori);
title('input');
subplot(2,2,2); imshow(fs);
title('segmentation result');
subplot(2,2,3);BW = boundarymask(Seg_M);imshow(imoverlay(fs,BW,'w'),'InitialMagnification',100);
title('boundary line');
subplot(2,2,4);BW_gt = boundarymask(Gt_M);imshow(imoverlay(fs2,BW_gt,'g'),'InitialMagnification',100);
title('groundtruth');




