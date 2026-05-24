%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                  %
%                Nonlinear Model DH                %
% Author: Huascar Mirko Montecinos Cortez          %
% Technological Institute of Aeronautics - ITA     %
%                                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Technological Institute of Aeronautics
% Electronic Devices and Systems (EEC-D)
% Copyright 2026 Regents of the Technological Institute of Aeronautics.
% All rights reserved.

close all
clear all
bdclose all
clc

%% Load information
load(["quad_model_v3.mat"])
alt = [1.75 2 2.25 2.5 2.75 3 3.25 3.5 3.75 4];
indice = 1;

%%
open('quad_model_v4')

for ii = 1:1:10
    alt(ii);
    sim('quad_model_v4');
    indice = indice + 1;
end


    