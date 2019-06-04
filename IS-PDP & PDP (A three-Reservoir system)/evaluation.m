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

%% calculate objective value
function E=evaluation(i,j,k,x,y,runoffly,runoffla,runoffaj,runoffjl,...
           a3zv_ly,a2zv_ly,a1zv_ly,a0zv_ly,a3vz_ly,a2vz_ly,a1vz_ly,a0vz_ly,a4Qzd_ly,a3Qzd_ly,a2Qzd_ly,a1Qzd_ly,a0Qzd_ly,...
           a3zv_ah,a2zv_ah,a1zv_ah,a0zv_ah,a3vz_ah,a2vz_ah,a1vz_ah,a0vz_ah,a6Qzd_ah,a5Qzd_ah,a4Qzd_ah,a3Qzd_ah,a2Qzd_ah,a1Qzd_ah,a0Qzd_ah,a3zQmax_ah,a2zQmax_ah,a1zQmax_ah,a0zQmax_ah,...
           a3zv_jaq,a2zv_jaq,a1zv_jaq,a0zv_jaq,a3vz_jaq,a2vz_jaq,a1vz_jaq,a0vz_jaq,a6Qzd_jaq,a5Qzd_jaq,a4Qzd_jaq,a3Qzd_jaq,a2Qzd_jaq,a1Qzd_jaq,a0Qzd_jaq,a3zQmax_jaq,a2zQmax_jaq,a1zQmax_jaq,a0zQmax_jaq)

%% 1 LY Reservoir
%% Augmented Lagrangian function: E_ly
[E_ly,qoutly]=LYtraopt(i,j,k,x,y,runoffly,runoffla,...
                       a3zv_ly,a2zv_ly,a1zv_ly,a0zv_ly,a3vz_ly,a2vz_ly,a1vz_ly,a0vz_ly,a4Qzd_ly,a3Qzd_ly,a2Qzd_ly,a1Qzd_ly,a0Qzd_ly);

%% 2 AH Reservoir
runoffah=qoutly+runoffla;
[E_ah,qoutah]=ahtraopt(i,j,k,x,y,runoffah,runoffaj,...
                       a3zv_ah,a2zv_ah,a1zv_ah,a0zv_ah,a3vz_ah,a2vz_ah,a1vz_ah,a0vz_ah,a6Qzd_ah,a5Qzd_ah,a4Qzd_ah,a3Qzd_ah,a2Qzd_ah,a1Qzd_ah,a0Qzd_ah,a3zQmax_ah,a2zQmax_ah,a1zQmax_ah,a0zQmax_ah);
 
%% 3 JAQ Reservoir
runoffjaq=qoutah+runoffaj;
[E_jaq,qoutjaq]=jaqtraopt(i,j,k,x,y,runoffjaq,runoffjl,...
                          a3zv_jaq,a2zv_jaq,a1zv_jaq,a0zv_jaq,a3vz_jaq,a2vz_jaq,a1vz_jaq,a0vz_jaq,a6Qzd_jaq,a5Qzd_jaq,a4Qzd_jaq,a3Qzd_jaq,a2Qzd_jaq,a1Qzd_jaq,a0Qzd_jaq,a3zQmax_jaq,a2zQmax_jaq,a1zQmax_jaq,a0zQmax_jaq);

 E=E_ly+E_ah+E_jaq;

end