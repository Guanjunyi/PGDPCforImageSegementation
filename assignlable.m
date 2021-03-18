%% Assign lable to each pixel
function [CL] = assignlable(CL_sup,LSup)
pixelnum = length(LSup);
CL = zeros(pixelnum,1);
NCL = max(CL_sup);
for i = 1:NCL
    Same_sups = find(CL_sup==i);
    for j = 1:length(Same_sups)
        CL(LSup==Same_sups(j))=i;
    end
end