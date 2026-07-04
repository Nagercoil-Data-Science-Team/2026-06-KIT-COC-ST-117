clc;
clear;
close all;

%% Number of Experiments

N = 500;

%% Input Parameters

Temperature = 20 + (60-20).*rand(N,1);

Flow_Rate = 1 + (20-1).*rand(N,1);

Flow_Ratio = 1 + (5-1).*rand(N,1);

Ionizable_Lipid = 30 + (60-30).*rand(N,1);

DSPC = 5 + (20-5).*rand(N,1);

Cholesterol = 20 + (50-20).*rand(N,1);

PEG_Lipid = 0.5 + (5-0.5).*rand(N,1);

Drug_Loading = 1 + (20-1).*rand(N,1);

Surfactant = 0.1 + (5-0.1).*rand(N,1);

Solvent_Ratio = 10 + (90-10).*rand(N,1);

%% Store Table

LNP_Data = table(Temperature,...
                 Flow_Rate,...
                 Flow_Ratio,...
                 Ionizable_Lipid,...
                 DSPC,...
                 Cholesterol,...
                 PEG_Lipid,...
                 Drug_Loading,...
                 Surfactant,...
                 Solvent_Ratio);

writetable(LNP_Data,'LNP_Input_Data.xlsx');

disp('Dataset Created Successfully');