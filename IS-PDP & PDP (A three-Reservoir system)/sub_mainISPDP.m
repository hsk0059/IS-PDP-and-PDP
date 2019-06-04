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
% ndis: number of discrete state variables of each reservoir
% ncom:  number of combinations of discrete state variables of three reservoirs
% the number of combinations whose Manhattan Distance to initial solution is less than 4 is 27
ndis=30; ncom=27;  

for jj=1:ncom
    path{jj,1}=zeros(lenthfile,data);
end

%% the maximum benifit at current stage
for jj=1:27
    Etemp{jj,1}=zeros(lenthfile,1);% the maximum benifit at current stage i
    Etemp1{jj,1}=zeros(lenthfile,1);%the maximum benifit at last stage i-1
end

%% iter: 1,2 or 3; 3: the year has converged; 1 or 2: the year has not converged;
iter=ones(lenthfile,1); iter1=iter;

%% lenthfile1: total years; fini: the indicator to stop 'While' cycle
lenthfile1=lenthfile;fini=1;
m1=[1:1:lenthfile1]';tt=0; % tt number of iterations

%%  save initial objective value
E1{1,1}=EE;
%% z1:the matrix of lenthfile1*data, recording the optimal trajectory of all years at current iteration 
%% z: recording the optimal trajectory at current iteration of some years  which do not converge 
%% z2: temporal substitute of z
z1_ly=z_ly;z1_ah=z_ah;z1_jaq=z_jaq;

while(fini)
% x: combinations at i-1 stage;    y: combinations at i stage;  
% hsk1 hsk2 hsk3 hsk4 Err Err1 Err2 :intermediate variables;
    clear x y hsk1 hsk2 hsk3 hsk4 Err Err1 Err2 
    
    lenthfile=size(runoff_ly,1);% lenthfile: now record some years which do not converge at current iteration
    tt=tt+1; EE(:,tt)=E1{1,1};
    
    %% run IS-PDP algorithm
    ISPDP;
    
    E1{1,1}(m1,:)=E{1,1};
%%  change the form of the path
   for i=1:92
        for mm=1:27
            for jj=1:lenthfile
                path1{jj,1}(mm,i)=path{mm,1}(jj,i);
            end
        end
    end
%%  determine the new 'initial solution' according to the new path form
for m=1:lenthfile

        path2=path1{m,1};

        aa(1)=1;

          for i=92:-1:2
              aa(94-i)=path2(aa(93-i),i);
          end

          bb=fliplr(aa);

           for i=1:92
              z_ly(m,i)=y{bb(i),i}(m,1);z_ah(m,i)=y{bb(i),i}(m,2);z_jaq(m,i)=y{bb(i),i}(m,3);
           end

          %%  use scaling factor to control search step if the initial solution is optimal
          if(sum(bb==1)==92)
              iter(m)=iter(m)+1;
          end       


end             
%%
m3=-1;%index of year
iter2=iter;
for m=1:lenthfile
        %% judge whether the exact year should stop iteration;
            if(iter(m)==2)
                m3=m3+1;
                
                m2=m1(m-m3);%the actual year
                
                iter1(m2,1)=2;
                
                runoff_ly(m-m3,:)=[];runoff_la(m-m3,:)=[];runoff_aj(m-m3,:)=[];runoff_jl(m-m3,:)=[];
                
                z1_ly(m2,:)=z_ly(m-m3,:);z1_ah(m2,:)=z_ah(m-m3,:);z1_jaq(m2,:)=z_jaq(m-m3,:);
                               
                z_ly(m-m3,:)=[];z_ah(m-m3,:)=[];z_jaq(m-m3,:)=[];
                
                m1(m-m3)=[]; iter2(m-m3)=[];
                
            end
end
iter=iter2;
%% stop criteria?
    if(sum(iter1==2)==lenthfile1)
        fini=0;
    end
end
