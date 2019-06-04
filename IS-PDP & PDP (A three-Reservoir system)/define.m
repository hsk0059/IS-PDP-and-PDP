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

%% 64 years:1951-2014£¬92 datas per year:8.1-10.31.
global lenthfile data dt cl
lenthfile = 64; data=92; dt=24*3600; cl=1.0E8; 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%% 1 Li-Yuan (LY) Reservoir %%%%%%%%%%%%%%%%%%%%%
%% basic statistics of Li-Yuan (LY) reservoir
% vmin: reservoir storage corresponding to annual top of buffer pool;
% vdesign:reservoir storage corresponding to design flood water level;
% vip:reservoir storage corresponding to top of inactive pool;
global vmin_ly vdesign_ly  vip_ly
vmin_ly=5.5394;  vdesign_ly=7.2710; vip_ly=5.5394;
 
%% 
% zmin: annual top of buffer pool;  zdesign:design flood water level;
% zip: top of inactive pool
global zmin_ly zdesign_ly zip_ly
zmin_ly=1605; zdesign_ly=1618; zip_ly=1605;
 
%% qmin :minimum outflow; dtqmax: maximum rate of water level fluctuation
global dtqmax_ly qmin_ly
qmin_ly=2000; dtqmax_ly=2000;
 
%% 
% zbegin: initial water level of the beginning impoundment period:
% pp: firm output;
% zend: end water level of the final impoundment period:
global zbegin_ly pp_ly zend_ly
zbegin_ly=1605; pp_ly=1102; zend_ly=1618;

%% upper and lower bound of optimization method
lowh_ly=[[1605:13/60:1618]';ones(32,1)*zend_ly];
highh_ly=[ones(15,1)*1610;ones(16,1)*1613;ones(15,1)*1617;ones(15,1)*1618;ones(31,1)*zend_ly];

%%%%%%%%%%%%%%%%%%%%%%%%%%%% 2 A-Hai (AH) Reservoir %%%%%%%%%%%%%%%%%%%%%
%% basic statistics of A-Hai (AH) reservoir
global vmin_ah vdesign_ah  vip_ah
vmin_ah=5.8983;  vdesign_ah=8.0576; vip_ah=5.6682;
 
%% 
global zmin_ah zdesign_ah zip_ah
zmin_ah=1493.3; zdesign_ah=1504; zip_ah=1492;
 
%%
global dtqmax_ah qmin_ah
qmin_ah=2700; dtqmax_ah=2000;
 
%% 
global zbegin_ah pp_ah zend_ah
zbegin_ah=1493.3; pp_ah=914; zend_ah=1504;

%%
lowh_ah=[[1493.3:10.7/60:1504]';ones(31,1)*zend_ah];
highh_ah=[ones(15,1)*1497;ones(16,1)*1500;ones(15,1)*1503;ones(15,1)*1504;ones(31,1)*zend_ah];

%%%%%%%%%%%%%%%%%%%%%%%%%%%% 3 Jin-An-Qiao (JAQ) Reservoir %%%%%%%%%%%%%%%%
%% basic statistics of Jin-An-Qiao (JAQ) reservoir
global vmin_jaq vdesign_jaq  vip_jaq
vmin_jaq=6.9072;  vdesign_jaq=8.4701; vip_jaq=5.0045;
 
%% 
global zmin_jaq zdesign_jaq zip_jaq
zmin_jaq=1410; zdesign_jaq=1418; zip_jaq=1398;
 
%%
global dtqmax_jaq qmin_jaq
qmin_jaq=3000; dtqmax_jaq=2000;
 
%% 
global zbegin_jaq pp_jaq zend_jaq
zbegin_jaq=1410; pp_jaq=1351.3; zend_jaq=1418;

%%
lowh_jaq=[[1410:8/60:1418]';ones(31,1)*zend_jaq];
highh_jaq=[ones(15,1)*1415;ones(16,1)*1417;ones(15,1)*1418;ones(15,1)*1418;ones(31,1)*zend_jaq];