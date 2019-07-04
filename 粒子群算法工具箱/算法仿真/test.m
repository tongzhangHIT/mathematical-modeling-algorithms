clc;
clear all;
x_range=[-2,2];
y_range=[-2,2];
range=[x_range;y_range];
MaX_V=0.2*(range(:,2)-range(:,1));%粒子群搜寻速度
n=2;
PSOparams=[100 150 20 2 2 0.9 0.4 1500 1e-25 250 NaN 0 0];

pso_Trelea_vectorized('test_func',n,MaX_V,range,1,PSOparams); 