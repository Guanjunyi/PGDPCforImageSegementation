function [PTs Lp] = preproces(fig_gf,Lp,NSup)

Lp = Lp';
Lp = Lp(:);
[n,m,~] = size(fig_gf);
X_Y = (1:(n*m))';
x = ceil(X_Y/m);
y = mod((X_Y),m);y(y==0) = m;
PTs = zeros(NSup,5);

%% RGB
RGB = 0;
if RGB
    r1 = fig_gf(:,:,1)';r1 = r1(:);
    g1 = fig_gf(:,:,2)';g1 = g1(:);
    b1 = fig_gf(:,:,3)';b1 = b1(:);
    r_rate = [r1./g1 r1./b1];
    r_rate = min(r_rate');
    r_rate(isnan(r_rate)) = 0;
    r_rate(r_rate==inf) = 0;
    points_ori = [r1 g1 b1 x y];
else
%% 正常情况 lab
    fig_lab = rgb2lab(fig_gf);
    l1 = fig_lab(:,:,1)';l1 = l1(:);
    a1 = fig_lab(:,:,2)'; a1 = a1(:);
    b1 = fig_lab(:,:,3)';b1 = b1(:);
    points_ori = [l1 a1 b1 x y];
end

for i = 1:NSup
    P_inSub = points_ori(Lp==i,:);
    %         PTs(i,:) = median(PPP);
    PTs(i,:) = mean(P_inSub);
    sup_points{1,i}=(Lp==i);end
end