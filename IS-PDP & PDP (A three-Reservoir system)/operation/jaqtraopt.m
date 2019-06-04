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

function [E_jaq,qout_jaq]=jaqtraopt(i,j,k,x,y,runoffjaq,runoffjl,...
                          a3zv_jaq,a2zv_jaq,a1zv_jaq,a0zv_jaq,a3vz_jaq,a2vz_jaq,a1vz_jaq,a0vz_jaq,a6Qzd_jaq,a5Qzd_jaq,a4Qzd_jaq,a3Qzd_jaq,a2Qzd_jaq,a1Qzd_jaq,a0Qzd_jaq,a3zQmax_jaq,a2zQmax_jaq,a1zQmax_jaq,a0zQmax_jaq)
%% reiterate some variables in define.m file
lenthfile=size(runoffjaq,1); dt=24*3600;
npre_jaq=1351.3;
%% reiterate some variables in input1.m file
XZ_jaq=[94.7,1638.1;111,2400;125.9,2400];


%% define minimum, maximum outflow, respectively.
 
%% control the minimum outflow of each 10-day(or 11-day)
        if ((i>=1)&&(i<=11)) %8.1-8.10 
            qoutmin_jaq=2700;
        elseif((i>=12)&&(i<=32))%8.11-8.20
            qoutmin_jaq=2000;
        elseif((i>=33)&&(i<=63))%10.11-10.31
            qoutmin_jaq=2500;
        else
            qoutmin_jaq=1000;
        end
 
%% maximum downstream outflow
qmax_jaq=10000;
%% maximum discharge capacity
q0_jaq=a3zQmax_jaq*y{k,1}(:,3).^3+a2zQmax_jaq*y{k,1}(:,3).^2+a1zQmax_jaq*y{k,1}(:,3).^1+a0zQmax_jaq;


qoutmax_jaq=min(qmax_jaq,q0_jaq);
 
%% water balance equation
beginv_jaq=a3zv_jaq*x{j,1}(:,3).^3+a2zv_jaq*x{j,1}(:,3).^2+a1zv_jaq*x{j,1}(:,3).^1+a0zv_jaq;
endv_jaq=a3zv_jaq*y{k,1}(:,3).^3+a2zv_jaq*y{k,1}(:,3).^2+a1zv_jaq*y{k,1}(:,3).^1+a0zv_jaq;
qout_jaq=runoffjaq+(beginv_jaq-endv_jaq)*1.0E8/dt;
 
%% calculate output
vup_jaq=(beginv_jaq+endv_jaq)/2;
zup_jaq=a3vz_jaq*vup_jaq.^3+a2vz_jaq*vup_jaq.^2+a1vz_jaq*vup_jaq.^1+a0vz_jaq;
zdown_jaq=a6Qzd_jaq*qout_jaq.^6+a5Qzd_jaq*qout_jaq.^5+a4Qzd_jaq*qout_jaq.^4+a3Qzd_jaq*qout_jaq.^3+a2Qzd_jaq*qout_jaq.^2+a1Qzd_jaq*qout_jaq.^1+a0Qzd_jaq;
hlost_jaq=0;
    
hpower_jaq=zup_jaq-zdown_jaq-hlost_jaq;
    
%% KN,KF:efficency coefficient£»Qelse:discharge for other use£»    
Qelse=120;KN=8.4;KF=0.933;
 
N1=KN*(qout_jaq-Qelse).*hpower_jaq/1.0E3; 
N2=interp1(XZ_jaq(:,1),XZ_jaq(:,2),hpower_jaq);
 
power_jaq=min(N1,KF*N2);
power1_jaq=power_jaq;
 
%% penalty factors
subject=zeros(lenthfile,1);
subject1=zeros(lenthfile,1);
subject2=zeros(lenthfile,1);
subject3=zeros(lenthfile,1);
subject4=zeros(lenthfile,1); 

if(i<62)
    
    for jj=1:lenthfile
        
       %% penalty for maximum discharge limit
            if (qout_jaq(jj)>qoutmax_jaq(jj))
                subject1(jj)=(qout_jaq(jj)-qoutmax_jaq(jj))*1.0E8;    
          %% for a few wet years which make maximum discharge constraint
            %  violate water level fluctuation constraint,we choose smaller penalty
                if(runoffjaq(jj)>q0_jaq(jj))
                    subject1(jj)=(qout_jaq(jj)-q0_jaq(jj))*0.01*24;
                end
            end

       %% penalty for minimum discharge limit
            if(qout_jaq(jj)<qoutmin_jaq)
                subject2(jj)=(qoutmin_jaq-qout_jaq(jj))*1.0E8;
              %% for a few dry years which can not guarantee minimum discharge, we choose smaller penalty    
                if(runoffjaq(jj)<qoutmin_jaq+25)
                    subject2(jj)=(qoutmin_jaq-qout_jaq(jj))*0.01*24; 
                end 
            end
            
            if(qout_jaq(jj)+runoffjl(jj)<0)
                subject2(jj)=-(qout_jaq(jj)+runoffjl(jj))*1.0E28;
            end  

            if(qout_jaq(jj)<0)
                subject2(jj)=-(qout_jaq(jj))*1.0E28;
            end
             
       %% penalty for hydropower generation limits
            if(power1_jaq(jj)<npre_jaq)
                subject3(jj)=(npre_jaq-power1_jaq(jj))*1.0E8;
           %% for a few dry years which can not guarantee firm output, we choose smaller penalty     
                if(runoffjaq(jj)<qoutmin_jaq)
                   subject3(jj)=(npre_jaq-power1_jaq(jj))*0.01*24;
                end      
            end

       %% penalty for water level fluctuation rate
            deltaz(jj)=y{k,1}(jj,3)-x{j,1}(jj,3);% water level fluctuation at adjacent periods

            if(i<32)
       %% control earlier impoundment process faster
                if(deltaz(jj)<0)
                    subject4(jj)=-(deltaz(jj))*1.0E28*24;          
                elseif(deltaz(jj)>0.4)
                    subject4(jj)=(deltaz(jj)-0.4)*1E28*24;
                end
       %% control post impoundment process slower
            else
                if(deltaz(jj)<0)
                    subject4(jj)=-(deltaz(jj))*1.0E28*24;      
                elseif(deltaz(jj)>0.3)
                    subject4(jj)=(deltaz(jj)-0.3)*1E28*24;
                end   
            end
    end
end

subject=subject+subject1+subject2+subject3+subject4;      
E_jaq=(power_jaq.*24-subject)/1E5;
end
 
