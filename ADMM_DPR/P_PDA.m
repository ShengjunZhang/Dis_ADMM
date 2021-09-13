clear; clc; close all;

load('opt_data.mat');

T =3000;
beta=1;
[sq_grad_Prox_PDA, xminuxbar_Prox_PDA, ~, time_Prox_PDA] =  Prox_PDA_beta(x0, edge_index,A,B,d,n_agents,y_all, a_Re_all, a_Im_all, y, a_Re, a_Im,T, beta);
