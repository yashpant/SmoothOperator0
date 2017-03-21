clc;clear all;close all;
%% Run planar example in robustness maximization (SOP-R) mode:
example2_toy_1shot;

%% Run building example in SOP-R mode:
BuildingMain_u_variable; %will make some garbage plots

%% Run building example in online ShrinkingHorizon mode
Bldg_Receding; 

%% Run dual quadrotor example
CaseMain_u_var;