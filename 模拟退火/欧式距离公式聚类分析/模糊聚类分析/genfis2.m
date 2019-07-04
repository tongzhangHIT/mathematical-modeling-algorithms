clc
clear all
close all
tripdata
subplot(211),plot(datin);
subplot(212),plot(datout);
fismat=genfis2(datin,datout,0.5);
fuzout=evalfis(datin,fismat);
trnRMSE=norm(fuzout-datout)/sqrt(length(fuzout));
trnRMSE=0.5276;
figure,
plot(datout,'o');
hold on
plot(fuzout)