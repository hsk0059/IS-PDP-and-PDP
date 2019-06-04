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

for i=1:data    
    %% determine the state vector at stage 1
        if(i==1)
            kmax=1;
            x{1,1}=ones(lenthfile,1)*[zbegin_ly,zbegin_ah,zbegin_jaq];
            y{1,1}=[z_ly(:,i),z_ah(:,i),z_jaq(:,i)];
       
        elseif(i<61)%
            kmax=ncom; % from ith stage (i>2), the number of combinations are ncom(ncom=27);
%%
            table=abc3;% for three resevoirs
            
            for jj=1:kmax
                y{jj,i}(:,1:3)=[z_ly(:,i),z_ah(:,i),z_jaq(:,i)];
            end

          %% sample all combinations whose Manhattan Distance to initial solution is less than 4
            hsk4(1,:)=table(end,:); 
            hsk4(2:7,:)=table(21:26,:);
            hsk4(8:19,:)=table(9:20,:);
            hsk4(20:27,:)=table(1:8,:);
%%
           for jj=1:kmax
               hsk3{jj,1}=ones(lenthfile,1)*hsk4(jj,:);

               for jjj=1:lenthfile
                   hdis_ly(jjj)=(highh_ly(i)-lowh_ly(i))/(ndis-1)*(0.5)^(iter(jjj)-1);%discrete value of one-unit of LY
                   hdis_ah(jjj)=(highh_ah(i)-lowh_ah(i))/(ndis-1)*(0.5)^(iter(jjj)-1);%discrete value of one-unit of AH
                   hdis_jaq(jjj)=(highh_jaq(i)-lowh_jaq(i))/(ndis-1)*(0.5)^(iter(jjj)-1);%discrete value of one-unit of JAQ

                   hsk1{jjj,1}=zeros(2);hsk1{jjj,1}(1,1)=hdis_ly(jjj);hsk1{jjj,1}(2,2)=hdis_ah(jjj);hsk1{jjj,1}(3,3)=hdis_jaq(jjj);
                  %% Err:disturbed value
                   Err{jj,1}(jjj,:)=hsk3{jj,1}(jjj,:)*hsk1{jjj,1};
               end

               y{jj,i}=y{jj,i}+Err{jj,1};

               y{jj,i}(:,1)=max(y{jj,i}(:,1),lowh_ly(i));y{jj,i}(:,1)=min(y{jj,i}(:,1),highh_ly(i));

               y{jj,i}(:,2)=max(y{jj,i}(:,2),lowh_ah(i)); y{jj,i}(:,2)=min(y{jj,i}(:,2),highh_ah(i));
               
               y{jj,i}(:,3)=max(y{jj,i}(:,3),lowh_jaq(i)); y{jj,i}(:,3)=min(y{jj,i}(:,3),highh_jaq(i));

           end                 
      
        else
            kmax=1;

            y{1,i}=[z_ly(:,i),z_ah(:,i),z_jaq(:,i)];
        end

       Etemp1= Etemp;
       
       for j=1:kmax
            E{j,1}=-ones(lenthfile,1)*1E30;
       end
       
%%
      runoffly=runoff_ly(:,i);runoffla=runoff_la(:,i);runoffaj=runoff_aj(:,i);runoffjl=runoff_jl(:,i);

%%           
        if(i==1)
            t=y;k=1;j=1;
            x{1,1}=ones(lenthfile,1)*[zbegin_ly,zbegin_ah,zbegin_jaq];
            Etemp{k,1}=evaluation(i,j,k,x,t,runoffly,runoffla,runoffaj,runoffjl,...
                        a3zv_ly,a2zv_ly,a1zv_ly,a0zv_ly,a3vz_ly,a2vz_ly,a1vz_ly,a0vz_ly,a4Qzd_ly,a3Qzd_ly,a2Qzd_ly,a1Qzd_ly,a0Qzd_ly,...
                        a3zv_ah,a2zv_ah,a1zv_ah,a0zv_ah,a3vz_ah,a2vz_ah,a1vz_ah,a0vz_ah,a6Qzd_ah,a5Qzd_ah,a4Qzd_ah,a3Qzd_ah,a2Qzd_ah,a1Qzd_ah,a0Qzd_ah,a3zQmax_ah,a2zQmax_ah,a1zQmax_ah,a0zQmax_ah,...
                        a3zv_jaq,a2zv_jaq,a1zv_jaq,a0zv_jaq,a3vz_jaq,a2vz_jaq,a1vz_jaq,a0vz_jaq,a6Qzd_jaq,a5Qzd_jaq,a4Qzd_jaq,a3Qzd_jaq,a2Qzd_jaq,a1Qzd_jaq,a0Qzd_jaq,a3zQmax_jaq,a2zQmax_jaq,a1zQmax_jaq,a0zQmax_jaq);

                       for m=1:lenthfile
                            if(Etemp{k,1}(m,1)>E{k,1}(m,1))
                                E{k,1}(m,1)=Etemp{k,1}(m,1);%when i=1,t is the same;
                            end
                       end

        elseif(i==2)
           x{1,1}=y{1,i-1}; 

           for jj=1:kmax
                t{jj,1}=y{jj,2};
           end

           j=1;

            parfor k=1:kmax  
                    Etemp{k,1}=Etemp1{1,1}+evaluation(i,j,k,x,t,runoffly,runoffla,runoffaj,runoffjl,...
                                a3zv_ly,a2zv_ly,a1zv_ly,a0zv_ly,a3vz_ly,a2vz_ly,a1vz_ly,a0vz_ly,a4Qzd_ly,a3Qzd_ly,a2Qzd_ly,a1Qzd_ly,a0Qzd_ly,...
                                a3zv_ah,a2zv_ah,a1zv_ah,a0zv_ah,a3vz_ah,a2vz_ah,a1vz_ah,a0vz_ah,a6Qzd_ah,a5Qzd_ah,a4Qzd_ah,a3Qzd_ah,a2Qzd_ah,a1Qzd_ah,a0Qzd_ah,a3zQmax_ah,a2zQmax_ah,a1zQmax_ah,a0zQmax_ah,...
                                a3zv_jaq,a2zv_jaq,a1zv_jaq,a0zv_jaq,a3vz_jaq,a2vz_jaq,a1vz_jaq,a0vz_jaq,a6Qzd_jaq,a5Qzd_jaq,a4Qzd_jaq,a3Qzd_jaq,a2Qzd_jaq,a1Qzd_jaq,a0Qzd_jaq,a3zQmax_jaq,a2zQmax_jaq,a1zQmax_jaq,a0zQmax_jaq);

                     for m=1:lenthfile
                        if(Etemp{k,1}(m,1)>E{k,1}(m,1))
                            E{k,1}(m,1)=Etemp{k,1}(m,1);                   
                            path{k,1}(m,i)=j;%% record the optimal path at last stage 
                        end
                    end
            end
        else

            if(i<61)
                jmax=27;
                kmax=27;
            elseif(i<62)
                jmax=27;
                kmax=1;
            else
                jmax=1;
                kmax=1;
            end

            for jj=1:jmax
                x{jj,1}=y{jj,i-1};
            end

            for jj=1:kmax
                t{jj,1}=y{jj,i};
            end
           
            parfor k=1:kmax                            
                for j=1:jmax
                        Etemp{k,1}=Etemp1{j,1}+evaluation(i,j,k,x,t,runoffly,runoffla,runoffaj,runoffjl,...
                                    a3zv_ly,a2zv_ly,a1zv_ly,a0zv_ly,a3vz_ly,a2vz_ly,a1vz_ly,a0vz_ly,a4Qzd_ly,a3Qzd_ly,a2Qzd_ly,a1Qzd_ly,a0Qzd_ly,...
                                    a3zv_ah,a2zv_ah,a1zv_ah,a0zv_ah,a3vz_ah,a2vz_ah,a1vz_ah,a0vz_ah,a6Qzd_ah,a5Qzd_ah,a4Qzd_ah,a3Qzd_ah,a2Qzd_ah,a1Qzd_ah,a0Qzd_ah,a3zQmax_ah,a2zQmax_ah,a1zQmax_ah,a0zQmax_ah,...
                                    a3zv_jaq,a2zv_jaq,a1zv_jaq,a0zv_jaq,a3vz_jaq,a2vz_jaq,a1vz_jaq,a0vz_jaq,a6Qzd_jaq,a5Qzd_jaq,a4Qzd_jaq,a3Qzd_jaq,a2Qzd_jaq,a1Qzd_jaq,a0Qzd_jaq,a3zQmax_jaq,a2zQmax_jaq,a1zQmax_jaq,a0zQmax_jaq);

                     for m=1:lenthfile                     
                        if(Etemp{k,1}(m,1)>E{k,1}(m,1))
                            E{k,1}(m,1)=Etemp{k,1}(m,1);                   
                            path{k,1}(m,i)=j;%% record the optimal path at last stage 
                        end
                    end
               end
            end
      end
      Etemp=E;
    end