function [Seg_M] = bulidSegLab(fig,CL)
[n,m,~] = size(fig);
 Seg_M = reshape(CL,[m,n]);
 Seg_M = Seg_M';
end