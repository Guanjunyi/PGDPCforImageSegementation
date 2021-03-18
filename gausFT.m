function [fig_gf] = gausFT(fig,sigma)
[n,m,~]=size(fig);
gauF_w=max(n,m);
gausFilter=fspecial('gaussian',[gauF_w gauF_w],sigma);fig_gf=imfilter(fig,gausFilter,'replicate');
end