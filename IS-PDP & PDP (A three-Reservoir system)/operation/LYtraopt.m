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

function [E_ly,qout_ly]=LYtraopt(i,j,k,x,y,runoffly,runoffla,...
                        a3zv_ly,a2zv_ly,a1zv_ly,a0zv_ly,a3vz_ly,a2vz_ly,a1vz_ly,a0vz_ly,a4Qzd_ly,a3Qzd_ly,a2Qzd_ly,a1Qzd_ly,a0Qzd_ly)
%% reiterate some variables in define.m file
lenthfile=size(runoffly,1);dt=24*3600;
npre_ly=1102; % firm output;
%% reiterate some variables in input1.m file
XZ_ly=[85.9,1489.8;116,2280;300,2280];


%% define minimum, maximum outflow, respectively.

%% control the minimum outflow of each 10-day(or 11-day)
        if ((i>=1)&&(i<=11)) %8.1-8.10 
            qoutmin_ly=1500;
        elseif((i>=12)&&(i<=32))%8.11-8.20
            qoutmin_ly=2000;
        elseif((i>=33)&&(i<=63))%10.11-10.31
            qoutmin_ly=2500;
        else
            qoutmin_ly=1300;
        end

%% maximum downstream outflow 
 qmax_ly=8500;
%% maximum discharge capacity
% q0_ly=a3zQmax_ly*x^3+a2zQmax_ly*x^2+a1zQmax_ly*x^1+a0zQmax_ly;
q0_ly=10000;

qoutmax_ly=min(qmax_ly,q0_ly);

%% water balance equation
beginv_ly=a3zv_ly*x{j,1}(:,1).^3+a2zv_ly*x{j,1}(:,1).^2+a1zv_ly*x{j,1}(:,1).^1+a0zv_ly;
endv_ly=a3zv_ly*y{k,1}(:,1).^3+a2zv_ly*y{k,1}(:,1).^2+a1zv_ly*y{k,1}(:,1).^1+a0zv_ly;
qout_ly=runoffly+(beginv_ly-endv_ly)*1.0E8/dt;

%% calculate output
vup_ly=(beginv_ly+endv_ly)/2;
zup_ly=a3vz_ly*vup_ly.^3+a2vz_ly*vup_ly.^2+a1vz_ly*vup_ly.^1+a0vz_ly;
zdown_ly=a4Qzd_ly*qout_ly.^4+a3Qzd_ly*qout_ly.^3+a2Qzd_ly*qout_ly.^2+a1Qzd_ly*qout_ly.^1+a0Qzd_ly;
hlost_ly=0;
    
hpower_ly=zup_ly-zdown_ly-hlost_ly;
    
%% KN,KF:efficency coefficient£»Qelse:discharge for other use£»   
Qelse=120;KN=8.6;KF=0.933;

N1=KN*(qout_ly-Qelse).*hpower_ly/1.0E3; 
N2=interp1(XZ_ly(:,1),XZ_ly(:,2),hpower_ly);

power_ly=min(N1,KF*N2);
power1_ly=power_ly;

%% penalty factors
subject=zeros(lenthfile,1);
subject1=zeros(lenthfile,1);
subject2=zeros(lenthfile,1);
subject3=zeros(lenthfile,1);
subject4=zeros(lenthfile,1);

if(i<62)
    
    for jj=1:lenthfile

            %% penalty for maximum discharge limit
                if (qout_ly(jj)>qoutmax_ly)
                    subject1(jj)=(qout_ly(jj)-qoutmax_ly)*1.0E8;     
                end

            %% penalty for minimum discharge limit
            if(qout_ly(jj)<qoutmin_ly)
                subject2(jj)=(qoutmin_ly-qout_ly(jj))*1.0E8;
              %% for a few dry years which can not guarantee minimum discharge, we choose smaller penalty  
                if(runoffly(jj)<qoutmin_ly+40)
                    subject2(jj)=(qoutmin_ly-qout_ly(jj))*0.01*24;
                end    

            end
            
            if(qout_ly(jj)+runoffla(jj)<0)
               subject2(jj)=-1*(qout_ly(jj)+runoffla(jj))*1.0E28;
            end

            if(qout_ly(jj)<0)
               subject2(jj)=-(qout_ly(jj))*1.0E28;
            end

           %% penalty for hydropower generation limits
            if(power1_ly(jj)<npre_ly)
                subject3(jj)=(npre_ly-power1_ly(jj))*1.0E8;
              %% for a few dry years which can not guarantee firm output, we choose smaller penalty 
                if(runoffly(jj)<qoutmin_ly)
                   subject3(jj)=(npre_ly-power1_ly(jj))*0.01*24;
                end      
            end

            %% penalty for water level fluctuation rate            
            deltaz(jj)=y{k,1}(jj,1)-x{j,1}(jj,1);% water level fluctuation at adjacent periods

            if(i<32)
            %% control earlier impoundment process slower 
                if(deltaz(jj)<0)
                    subject4(jj)=-(deltaz(jj))*1.0E28*24;           
                elseif(deltaz(jj)>0.5)
                    subject4(jj)=(deltaz(jj)-0.5)*1E28*24;
                end
            %% control post impoundment process faster 
            else
                if(deltaz(jj)<0)
                    subject4(jj)=-(deltaz(jj))*1.0E28*24; 
                elseif(deltaz(jj)>0.6)
                    subject4(jj)=(deltaz(jj)-0.6)*1E28*24;
                end   
            end
    end
end

subject=subject+subject1+subject2+subject3+subject4;                
E_ly=(power_ly.*24-subject)/1E5;
end
