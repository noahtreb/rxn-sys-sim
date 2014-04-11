fileName = 'Schlogl2-new.nc';
speciesId = 1;

t = ncread(fileName, 'time');
d1 = ncread(fileName, 'initFwdData');
d2 = ncread(fileName, 'initRevData');
d3 = ncread(fileName, 'fwdData');
d4 = ncread(fileName, 'revData');
numDataSavePts = ncread(fileName, 'numDataSavePts');

yMax = 1000; %= max([max(max(d1(speciesId, :, :)))...
    %max(max(d2(speciesId, :, :)))...
    %max(max(max(d3(speciesId, :, :, :))))...
    %max(max(max(d4(speciesId, :, :, :))))]);

figure;
subplot(1, 2, 1);
plot(t, squeeze(d1(speciesId, :, :)));
ylim([0 yMax]);

subplot(1, 2, 2);
plot(t, squeeze(d2(speciesId, :, :)));
ylim([0 yMax]);

%%

for i = 1:numDataSavePts
    figure;
    subplot(1, 2, 1);
    plot(t, squeeze(d3(speciesId, :, :, i)));
    ylim([0 yMax]);
    
    subplot(1, 2, 2);
    plot(t, squeeze(d4(speciesId, :, :, i)));
    ylim([0 yMax]);
end