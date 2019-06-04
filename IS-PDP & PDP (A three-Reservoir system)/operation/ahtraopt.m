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

function [E_ah,qout_ah]=ahtraopt(i,j,k,x,y,runoffah,runoffaj,...
                        a3zv_ah,a2zv_ah,a1zv_ah,a0zv_ah,a3vz_ah,a2vz_ah,a1vz_ah,a0vz_ah,a6Qzd_ah,a5Qzd_ah,a4Qzd_ah,a3Qzd_ah,a2Qzd_ah,a1Qzd_ah,a0Qzd_ah,a3zQmax_ah,a2zQmax_ah,a1zQmax_ah,a0zQmax_ah)
%% reiterate some variables in define.m file
lenthfile=size(runoffah,1); dt=24*3600;
npre_ah=914; % firm output;
%% reiterate some variables in input1.m file
XZ_ah=[63,1309.1;77,2000;95,2000];


%% define minimum, maximum outflow, respectively.

%% control the minimum outflow of each 10-day(or 11-day)
        if ((i>=1)&&(i<=11)) %8.1-8.10 
            qoutmin_ah=2700;
        elseif((i>=12)&&(i<=32))%8.11-8.20
            qoutmin_ah=2000;
        elseif((i>=33)&&(i<=63))%10.11-10.31
            qoutmin_ah=2500;
        else
            qoutmin_ah=1000;
        end
 
%% maximum downstream outflow
 qmax_ah=9500;
%% maximum discharge capacity
q0_ah=a3zQmax_ah*y{k,1}(:,2).^3+a2zQmax_ah*y{k,1}(:,2).^2+a1zQmax_ah*y{k,1}(:,2).^1+a0zQmax_ah;


qoutmax_ah=min(qmax_ah,q0_ah);
 
%% water balance equation
beginv_ah=a3zv_ah*x{j,1}(:,2).^3+a2zv_ah*x{j,1}(:,2).^2+a1zv_ah*x{j,1}(:,2).^1+a0zv_ah;
endv_ah=a3zv_ah*y{k,1}(:,2).^3+a2zv_ah*y{k,1}(:,2).^2+a1zv_ah*y{k,1}(:,2).^1+a0zv_ah;
qout_ah=runoffah+(beginv_ah-endv_ah)*1.0E8/dt; 
 
%% calculate output
vup_ah=(beginv_ah+endv_ah)/2;
zup_ah=a3vz_ah*vup_ah.^3+a2vz_ah*vup_ah.^2+a1vz_ah*vup_ah.^1+a0vz_ah;
zdown_ah=a6Qzd_ah*qout_ah.^6+a5Qzd_ah*qout_ah.^5+a4Qzd_ah*qout_ah.^4+a3Qzd_ah*qout_ah.^3+a2Qzd_ah*qout_ah.^2+a1Qzd_ah*qout_ah.^1+a0Qzd_ah;
hlost_ah=0;
    
hpower_ah=zup_ah-zdown_ah-hlost_ah;
    
%% KN,KF:efficency coefficient£»Qelse:discharge for other use£»   
Qelse=120;KN=8.6;KF=0.933;
 
N1=KN*(qout_ah-Qelse).*hpower_ah/1.0E3; 
N2=interp1(XZ_ah(:,1),XZ_ah(:,2),hpower_ah);
 
power_ah=min(N1,KF*N2);
power1_ah=power_ah;
 
%% penalty factors
subject=zeros(lenthfile,1);
subject1=zeros(lenthfile,1);
subject2=zeros(lenthfile,1);
subject3=zeros(lenthfile,1);
subject4=zeros(lenthfile,1);

if(i<62)
    
    for jj=1:lenthfile
        
       %% penalty for maximum discharge limit
        if (qout_ah(jj)>qoutmax_ah(jj))
            subject1(jj)=(qout_ah(jj)-qoutmax_ah(jj))*1.0E8;    
          %% for a few wet years which make maximum discharge constraint
            %  violate water level fluctuation constraint,we choose smaller penalty  
            if(runoffah(jj)>q0_ah(jj))
                    subject1(jj)=(qout_ah(jj)-q0_ah(jj))*0.01*24;
            end
        end

       %% penalty for minimum discharge limit
        if(qout_ah(jj)<qoutmin_ah)
            subject2(jj)=(qoutmin_ah-qout_ah(jj))*1.0E8;
              %% for a few dry years which can not guarantee minimum discharge, we choose smaller penalty      
            if(runoffah(jj)<qoutmin_ah+50)
                subject2(jj)=(qoutmin_ah-qout_ah(jj))*0.01*24;
            end
        end    

        if(qout_ah(jj)+runoffaj(jj)<0)
            subject2(jj)=-(qout_ah(jj)+runoffaj(jj))*1.0E28;
        end  

        if(qout_ah(jj)<0)
           subject2(jj)=-(qout_ah(jj))*1.0E28;
        end
     
       %% penalty for hydropower generation limits
        if(power1_ah(jj)<npre_ah)
            subject3(jj)=(npre_ah-power1_ah(jj))*1.0E8;
           %% for a few dry years which can not guarantee firm output, we choose smaller penalty   
            if(runoffah(jj)<qoutmin_ah)
               subject3(jj)=(npre_ah-power1_ah(jj))*0.01*24;
            end
        end


       %% penalty for water level fluctuation rate 
        deltaz(jj)=y{k,1}(jj,2)-x{j,1}(jj,2);% water level fluctuation at adjacent periods

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
 E_ah=(power_ah.*24-subject)/1E5;
end