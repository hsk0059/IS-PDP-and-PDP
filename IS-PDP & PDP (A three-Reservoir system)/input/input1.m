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

%% input1: input daily runoff, reservoir characteristics of each reservoir,initial solution;

%%%%%%%%%%%%%%%%%%%%%%%%%%%% 1 Li-Yuan (LY) Reservoir %%%%%%%%%%%%%%%%%%%%%
%% dialy runoff: Aug.1-Oct.31
global runoff_ly
runoff_ly=load('1LY8-1_10-31.txt');

%% Reservoir Characteristics
global zv_ly zq_ly QZDown_ly XZ_ly lose_ly 
zv_ly=load('lyZV.txt'); %zv
zq_ly=load('lyZQ.txt'); %zq
QZDown_ly=load('lyQZDown.txt');%QZDown
XZ_ly=load('lyXZ.txt');%hN
lose_ly=load('lylose.txt');%lose

%% Parameters: cerve fitting
global  a3zv_ly  a2zv_ly a1zv_ly  a0zv_ly 
global  a3vz_ly a2vz_ly a1vz_ly a0vz_ly
global  a4Qzd_ly a3Qzd_ly  a2Qzd_ly a1Qzd_ly a0Qzd_ly
global  a4zdQ_ly a3zdQ_ly a2zdQ_ly a1zdQ_ly a0zdQ_ly
global  a3zQmax_ly a2zQmax_ly a1zQmax_ly a0zQmax_ly
global  a6Qmaxz_ly a5Qmaxz_ly a4Qmaxz_ly a3Qmaxz_ly a2Qmaxz_ly a1Qmaxz_ly a0Qmaxz_ly 
 
p=polyfit(zv_ly(23:27,1),zv_ly(23:27,2),3);
a3zv_ly=p(1);a2zv_ly=p(2);a1zv_ly=p(3);a0zv_ly=p(4);
p=polyfit(zv_ly(23:27,2),zv_ly(23:27,1),3);
a3vz_ly=p(1); a2vz_ly=p(2); a1vz_ly=p(3); a0vz_ly=p(4); 
 
p=polyfit(QZDown_ly(5:11,1),QZDown_ly(5:11,2),4);
a4Qzd_ly=p(1); a3Qzd_ly=p(2); a2Qzd_ly=p(3);a1Qzd_ly=p(4);a0Qzd_ly=p(5);
p=polyfit(QZDown_ly(5:11,2),QZDown_ly(5:11,1),4);
a4zdQ_ly=p(1);  a3zdQ_ly=p(2);  a2zdQ_ly=p(3); a1zdQ_ly=p(4);a0zdQ_ly=p(5);
 
p=polyfit(zq_ly(2:12,1),zq_ly(2:12,2),3);
a3zQmax_ly=p(1);a2zQmax_ly=p(2);a1zQmax_ly=p(3);a0zQmax_ly=p(4);
p=polyfit(zq_ly(2:12,2),zq_ly(2:12,1),6);
a6Qmaxz_ly=p(1); a5Qmaxz_ly=p(2); a4Qmaxz_ly=p(3); a3Qmaxz_ly =p(4); a2Qmaxz_ly=p(5); a1Qmaxz_ly=p(6); a0Qmaxz_ly=p(7); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%% 2 A-Hai (AH) Reservoir %%%%%%%%%%%%%%%%%%%%%%%
global runoff_la
runoff_la=load('2la8-1_10-31.txt');%runoff of interval basin between LY-AH

%% Reservoir Characteristics
global zv_ah zq_ah QZDown_ah XZ_ah lose_ah 
zv_ah=load('ahZV.txt'); 
zq_ah=load('ahZQ.txt');
QZDown_ah=load('ahQZDown.txt');
XZ_ah=load('ahXZ.txt');
lose_ah=load('ahlose.txt');

%% Parameters
global  a3zv_ah  a2zv_ah a1zv_ah  a0zv_ah 
global  a3vz_ah a2vz_ah a1vz_ah a0vz_ah
global  a6Qzd_ah a5Qzd_ah a4Qzd_ah a3Qzd_ah  a2Qzd_ah a1Qzd_ah a0Qzd_ah
global  a3zdQ_ah a2zdQ_ah a1zdQ_ah a0zdQ_ah
global  a3zQmax_ah a2zQmax_ah a1zQmax_ah a0zQmax_ah
global  a3Qmaxz_ah a2Qmaxz_ah a1Qmaxz_ah a0Qmaxz_ah 

p=polyfit(zv_ah(17:21,1),zv_ah(17:21,2),3);
a3zv_ah=p(1);a2zv_ah=p(2);a1zv_ah=p(3);a0zv_ah=p(4);
p=polyfit(zv_ah(17:21,2),zv_ah(17:21,1),3);
a3vz_ah=p(1); a2vz_ah=p(2); a1vz_ah=p(3); a0vz_ah=p(4); 

p=polyfit(QZDown_ah(1:11,1),QZDown_ah(1:11,2),6);
a6Qzd_ah=p(1); a5Qzd_ah=p(2); a4Qzd_ah=p(3);a3Qzd_ah=p(4);a2Qzd_ah=p(5);a1Qzd_ah=p(6); a0Qzd_ah=p(7);
p=polyfit(QZDown_ah(1:11,2),QZDown_ah(1:11,1),3);
a3zdQ_ah=p(1);  a2zdQ_ah=p(2);  a1zdQ_ah=p(3); a0zdQ_ah=p(4);

p=polyfit(zq_ah(9:26,1),zq_ah(9:26,2),3);
a3zQmax_ah=p(1);a2zQmax_ah=p(2);a1zQmax_ah=p(3);a0zQmax_ah=p(4);
p=polyfit(zq_ah(9:26,2),zq_ah(9:26,1),3);
a3Qmaxz_ah=p(1); a2Qmaxz_ah=p(2); a1Qmaxz_ah=p(3); a0Qmaxz_ah =p(4); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%% 3 Jin-An-Qiao (JAQ) Reservoir %%%%%%%%%%%%%%%%
global runoff_aj
runoff_aj=load('3aj8-1_10-31.txt');%runoff of interval basin between AH-JAQ

%% Reservoir Characteristics
global zv_jaq zq_jaq QZDown_jaq XZ_jaq lose_jaq 
zv_jaq=load('jaqZV.txt');
zq_jaq=load('jaqZQ.txt');
QZDown_jaq=load('jaqQZDown.txt');
XZ_jaq=load('jaqXZ.txt');
lose_jaq=load('jaqlose.txt');
 
%% Parameters
global  a3zv_jaq  a2zv_jaq a1zv_jaq  a0zv_jaq 
global  a3vz_jaq a2vz_jaq a1vz_jaq a0vz_jaq
global  a6Qzd_jaq a5Qzd_jaq a4Qzd_jaq a3Qzd_jaq  a2Qzd_jaq a1Qzd_jaq a0Qzd_jaq
global   a5zdQ_jaq a4zdQ_jaq a3zdQ_jaq  a2zdQ_jaq a1zdQ_jaq a0zdQ_jaq
global  a3zQmax_jaq a2zQmax_jaq a1zQmax_jaq a0zQmax_jaq
global  a6Qmaxz_jaq a5Qmaxz_jaq a4Qmaxz_jaq a3Qmaxz_jaq a2Qmaxz_jaq a1Qmaxz_jaq a0Qmaxz_jaq 
 
p=polyfit(zv_jaq(:,1),zv_jaq(:,2),3);
a3zv_jaq=p(1);a2zv_jaq=p(2);a1zv_jaq=p(3);a0zv_jaq=p(4);
p=polyfit(zv_jaq(:,2),zv_jaq(:,1),3);
a3vz_jaq=p(1); a2vz_jaq=p(2); a1vz_jaq=p(3); a0vz_jaq=p(4); 
 
p=polyfit(QZDown_jaq(1:21,1),QZDown_jaq(1:21,2),6);
a6Qzd_jaq=p(1); a5Qzd_jaq=p(2); a4Qzd_jaq=p(3);a3Qzd_jaq=p(4);a2Qzd_jaq=p(5);a1Qzd_jaq=p(6); a0Qzd_jaq=p(7);
p=polyfit(QZDown_jaq(1:21,2),QZDown_jaq(1:21,1),5);
a5zdQ_jaq=p(1); a4zdQ_jaq=p(2);a3zdQ_jaq=p(3);  a2zdQ_jaq=p(4);  a1zdQ_jaq=p(5); a0zdQ_jaq=p(6);
 
p=polyfit(zq_jaq(:,1),zq_jaq(:,2),3);
a3zQmax_jaq=p(1);a2zQmax_jaq=p(2);a1zQmax_jaq=p(3);a0zQmax_jaq=p(4);
p=polyfit(zq_jaq(:,2),zq_jaq(:,1),6);
a6Qmaxz_jaq=p(1); a5Qmaxz_jaq=p(2); a4Qmaxz_jaq=p(3);  a3Qmaxz_jaq=p(4); a2Qmaxz_jaq=p(5); a1Qmaxz_jaq=p(6); a0Qmaxz_jaq =p(7);

%%%%%%%%%%%%%%%%%%%%%%%%%%%% 4 Long-Kai-Kou (LKK) Reservoir %%%%%%%%%%%%%%%
global runoff_jl
runoff_jl=load('4jl8-1_10-31.txt');%runoff of interval basin between JAQ-LKK

clc;

%% define initial solution
z_ly=load ('z1_ly.mat');z_ly=z_ly.z1_ly';z_ah=load ('z1_ah.mat');z_ah=z_ah.z1_ah';
z_jaq=load ('z1_jaq.mat');z_jaq=z_jaq.z1_jaq';

%%  EE: the objective value of initial solution
     for i=1:92
         
         runoffly=runoff_ly(:,i);runoffla=runoff_la(:,i);runoffaj=runoff_aj(:,i);runoffjl=runoff_jl(:,i);
%%
          if(i==1)
              x{1,1}=ones(lenthfile,1)*[zbegin_ly,zbegin_ah,zbegin_jaq];
              y0{1,1}=[z_ly(:,i),z_ah(:,i),z_jaq(:,i)];
              j=1;k=1;
                EE=evaluation(i,j,k,x,y0,runoffly,runoffla,runoffaj,runoffjl,...
                    a3zv_ly,a2zv_ly,a1zv_ly,a0zv_ly,a3vz_ly,a2vz_ly,a1vz_ly,a0vz_ly,a4Qzd_ly,a3Qzd_ly,a2Qzd_ly,a1Qzd_ly,a0Qzd_ly,...
                    a3zv_ah,a2zv_ah,a1zv_ah,a0zv_ah,a3vz_ah,a2vz_ah,a1vz_ah,a0vz_ah,a6Qzd_ah,a5Qzd_ah,a4Qzd_ah,a3Qzd_ah,a2Qzd_ah,a1Qzd_ah,a0Qzd_ah,a3zQmax_ah,a2zQmax_ah,a1zQmax_ah,a0zQmax_ah,...
                    a3zv_jaq,a2zv_jaq,a1zv_jaq,a0zv_jaq,a3vz_jaq,a2vz_jaq,a1vz_jaq,a0vz_jaq,a6Qzd_jaq,a5Qzd_jaq,a4Qzd_jaq,a3Qzd_jaq,a2Qzd_jaq,a1Qzd_jaq,a0Qzd_jaq,a3zQmax_jaq,a2zQmax_jaq,a1zQmax_jaq,a0zQmax_jaq);
          else
              x{1,1}=[z_ly(:,i-1),z_ah(:,i-1),z_jaq(:,i-1)];
              y0{1,1}=[z_ly(:,i),z_ah(:,i),z_jaq(:,i)];
              j=1;k=1;
                EE=EE+evaluation(i,j,k,x,y0,runoffly,runoffla,runoffaj,runoffjl,...
                        a3zv_ly,a2zv_ly,a1zv_ly,a0zv_ly,a3vz_ly,a2vz_ly,a1vz_ly,a0vz_ly,a4Qzd_ly,a3Qzd_ly,a2Qzd_ly,a1Qzd_ly,a0Qzd_ly,...
                        a3zv_ah,a2zv_ah,a1zv_ah,a0zv_ah,a3vz_ah,a2vz_ah,a1vz_ah,a0vz_ah,a6Qzd_ah,a5Qzd_ah,a4Qzd_ah,a3Qzd_ah,a2Qzd_ah,a1Qzd_ah,a0Qzd_ah,a3zQmax_ah,a2zQmax_ah,a1zQmax_ah,a0zQmax_ah,...
                        a3zv_jaq,a2zv_jaq,a1zv_jaq,a0zv_jaq,a3vz_jaq,a2vz_jaq,a1vz_jaq,a0vz_jaq,a6Qzd_jaq,a5Qzd_jaq,a4Qzd_jaq,a3Qzd_jaq,a2Qzd_jaq,a1Qzd_jaq,a0Qzd_jaq,a3zQmax_jaq,a2zQmax_jaq,a1zQmax_jaq,a0zQmax_jaq);
          end
     end

%% in this case, the authors only choose some points whose Manhattan Distance to initial solution is less than 4
abc3=load ('abc3.mat');abc3=abc3.abc;