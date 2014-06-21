%% Read the data.

fileName = 'Schlogl1-new6.nc';
speciesId = 2;

t = ncread(fileName, 'time');
d1 = ncread(fileName, 'initFwdData');
d2 = ncread(fileName, 'initRevData');
d3 = ncread(fileName, 'fwdData');
d4 = ncread(fileName, 'revData');
a1 = ncread(fileName, 'initFwdAbsCurr');
a2 = ncread(fileName, 'initRevAbsCurr');
a3 = ncread(fileName, 'fwdAbsCurr');
a4 = ncread(fileName, 'revAbsCurr');
numDataSavePts = ncread(fileName, 'numDataSavePts');

yMax = max([max(max(d1(speciesId, :, :)))...
    max(max(d2(speciesId, :, :)))...
    max(max(max(d3(speciesId, :, :, :))))...
    max(max(max(d4(speciesId, :, :, :))))]);

%% Produce intial data plots.

figure;
subplot(1, 2, 1);
plot(t, squeeze(d1(speciesId, :, :)));
ylim([0 yMax]);

subplot(1, 2, 2);
plot(t, squeeze(d2(speciesId, :, :)));
ylim([0 yMax]);

%% Produce refined data plots.

for i = 1:numDataSavePts
    figure;
    subplot(1, 2, 1);
    plot(t, squeeze(d3(speciesId, :, :, i)));
    ylim([0 yMax]);
    
    subplot(1, 2, 2);
    plot(t, squeeze(d4(speciesId, :, :, i)));
    ylim([0 yMax]);
end

%% Produce initial plots of the absorbing current.

figure;
subplot(1, 2, 1);
plot(t, a1(speciesId, :), 'or');

subplot(1, 2, 2);
plot(t, a2(speciesId, :), 'or');