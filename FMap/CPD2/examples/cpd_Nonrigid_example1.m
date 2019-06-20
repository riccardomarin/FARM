%  Nonrigid Example1. Coherent Point Drift (CPD).
%  Registration of 2D fish point sets without noise and outliers.

clear all; close all; clc;
addpath('../data/')
addpath('../core/mex')
addpath('../core/utils')
addpath('../core/Rigid')
addpath('../core/Nonrigid')
addpath('../core/')

load cpd_data2D_fish.mat

opt.method='nonrigid';
[Transform, C]=cpd_register(X,Y, opt);

figure,cpd_plot_iter(X, Y); title('Before');
figure,cpd_plot_iter(X, Transform.Y);  title('After registering Y to X');
