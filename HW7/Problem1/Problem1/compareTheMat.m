coMat = textread('coVarMat.txt');

sample = textread('mulSampleMat.txt');

h = figure(1);
subplot(2,1,1);
imagesc(coMat);
title('covariance matrix');
colorbar;
subplot(2,1,2);
imagesc(sample);
title('covariance matrix of 10000 draws');
set(h, 'Position', [100, 100, 500, 1000]);
colorbar;