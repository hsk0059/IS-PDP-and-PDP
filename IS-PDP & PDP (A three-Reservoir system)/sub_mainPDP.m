%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Importance Sampling - Parallel Dynamic Programming (IS-PDP) Toolbox
%
% This toolbox purpose is to solve the curse of dimensionality of multi-reservoir 
% impoundment operation, in an example of  three-reservoir system (Li-Yuan,A-hai,
% Jin-An-Qiao) to close the gap between research studies and real-world applications. 
% The toolbox, called Importance Sampling-Parallel Dynamic Programming(IS-PDP), allows 
% users to seek for an suboptimal trajectory of deterministic problem through the 
% IS-PDP method, which could balance computational efficiency and solution quality.

% Developers: Shaokun He, Shenglian Guo & Kebing Chen
% contact email: he_shaokun@whu.edu.cn
% Contact address: State Key Laboratory of Water Resources and Hydropower Engineering 
% Science, Wuhan University, Wuhan, 430072, China
% Year first available: 2019

% A synthetic Yangtze dataset of similar complexity was generated because the authors 
% cannot get permission to share the Yangtze dataset.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% define some basis of the state combinations
% num: number of reservoirs
% ndis: number of discrete state variables of each reservoir
% ncom:  number of combinations of discrete state variables of three reservoirs
num=3;ndis=30; ncom=ndis^num;

for jj=1:ncom
    path{jj,1}=zeros(lenthfile,92);
end

%% the maximum benifit at current stage
for jj=1:ncom
    Etemp{jj,1}=zeros(lenthfile,1);% the maximum benifit at current stage i
    Etemp1{jj,1}=zeros(lenthfile,1);%%the maximum benifit at last stage i-1
end
PDP;