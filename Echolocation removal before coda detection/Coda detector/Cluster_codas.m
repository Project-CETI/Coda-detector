function [L_hat_inds,L_hat,U_max]=Cluster_codas(sim,IPI,Locs,Pks)

%   AUTHOR:         Guy Gubnitsky
%   DATE:           February 2024
%   DESCRIPTION:

%   This function gets a vector of M transients and outputs a vector
%   contains the indices of suspected sperm whale coda clicks 
%   and a  vector containing the detection scores. 
%   The function is based on the clustering method described in the paper:
%   "Automatic Detection and Annotation of Sperm Whale Codas: a test case for
%    the Dominica Island"

%   INPUT:
%   > sim                   - Matrix of MXM elements with similarity measures between every pair of identified transients
%   > IPI                   - Vector of 1XM with IPI estimations of each identified transient.
%   > Locs                  - Vector of 1XM with time of arrival in seconds for M identified transients
%   > Pks                   - Vector of 1XM containing the peaks of M identified transients

%   OUTPUT:
%   > L_hat_inds            - Cell array containing vectors of 1XW with indices of suspected sperm whale coda clicks
%   > L_hat                 - Cell array containing vectors of KXWk with K detected codas with Wk time of arrivals for the clicks' of the k'th coda. 
%   > U_max                 - Vector of 1XK representing the utility score associated with each detected coda.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Import parameters
load('EigV.mat');
load('EigV_3D.mat');

% load('Params_mix.mat');
load('Params.mat');
load('Params3D_rev.mat');


%% Main parameters:

Reduce_complexity=1; %Optional: use ICI constraints to reduce program run time.
if Reduce_complexity
    ICI_Min=0.05; % Minimum allowed ICI 
    ICI_Max=0.6;   % Maximum allowed ICI
    Consi_max=[0 0 1 0.83 0.65 0.5 0.63 0.5 0.71 0.4];  % Maximum allowed of median of Consistency
    Consi_3_max=[0 0 1 0.8 0.7 0.93 0.99 0.83 1.1 0.74];  % Maximum allowed of median of Consistency
end
    alpha2=15;%3/5;    % Normalizing factor to weight the penalty over the cluster’s rank
    alpha1=1e3;    % Normalizing factor to weight the penalty over the cluster’s temporal likelihood
    alpha3=0.1;
    C_lim=1e5;  % Capacity limit- limit over the maximum matrix rank matlab can handle. Larger matrices will divided by loops.
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%%%%%%%%% Clustering and Detection%%%%%%%%%%%%%%%

% Variables initialization
k=0;   % initialize clustering iteration index
Ls=Locs;
L=length(Ls);
S=[sim ; zeros(1,length(Locs))]; % intialize feature similarity matrix
L_hat={};  U_max=[]; L_hat_inds={};

while(1)

    k=k+1;    

    %% Clustering initialization:

        Score=-100*ones(1,8);
        Corresponing_clusters={};
        All_clusters=[]; All_U=[];
    
    %% Clustering

    
    if length(Locs)<8
        Max_rank=length(Locs);
    else
        Max_rank=8; %set the maximum cluster size to 8 clicks (since we dont have enough codas of higher number of clicks in our current legacy database, such that rhythmic features can be reliably estimated)
    end

    for i=3:Max_rank  % evaluate likelihood scores on clusters of size ranges between 3 clicks and Max_rank clicks
         
         v = 1:length(Locs);
         AllCombs = nchoosek(v,i);
         AllCombs_bank(i).comb=AllCombs; % store all index combination of clicks in clusters of size i
      %% Optional for reducing runtime: ICI extreame values calculation for removing irrelevant clusters 
         ICI_g_max=max(diff(Locs(AllCombs),1,2),[],2); % calculate the max ICI value of each candidate cluster
         ICI_g_min=min(diff(Locs(AllCombs),1,2),[],2); % calculate the min ICI value of each candidate cluster
         ICI_mat=diff(Locs(AllCombs),1,2); 
         ICI_n1=ICI_mat(:,2:end); ICI_n0=ICI_mat(:,1:end-1);
      %% Optional for reducing runtime: ICI consistency calculations for removing irrelevant clusters 
         Consi=median(abs(log(ICI_n1./ICI_n0)),2);  
         Consi_3=abs(log(ICI_n1(:,end)./ICI_n0(:,end)));
      %% Pick candidate clusters
         if Reduce_complexity
             G_ind=find(ICI_g_min>ICI_Min & ICI_g_max<ICI_Max & Consi<Consi_max(i) & Consi_3<Consi_3_max(i)); 
         else
             G_ind=[1:size(AllCombs,1)];
         end
         AllCombs_mat=mat2inds(AllCombs,length(Locs)); % Convert clicks indices to a binary association hypothesis
         
         %% Calculation of the structural likelihood L_s
         
         Rank=i;                                   % clusters' rank (i.e., the number of clicks in a cluster)
         N_nodes=0.5*(Rank-1)*(Rank);              % number of nodes in the clusters         
         Capacity=size(AllCombs_mat,1)*size(AllCombs_mat,2)*size(S,1)*size(S,2);
     %% Use matrix calculation if the number of candidate cluster exceeds the limit specified in the parameter: C_lim
         if Capacity>C_lim
            C_sum=[];
            C_slices=floor(C_lim/(size(S,1)^3));
            C_iterations=floor(size(AllCombs_mat,1)/C_slices);
            for Ci=1:C_iterations
                C_sum=[C_sum ; diag(AllCombs_mat((Ci-1)*C_slices+1:Ci*C_slices,:)*S*AllCombs_mat((Ci-1)*C_slices+1:Ci*C_slices,:)')];
            end
                C_sum=[C_sum ; diag(AllCombs_mat(Ci*C_slices+1:end,:)*S*AllCombs_mat(Ci*C_slices+1:end,:)')];
        else
         C_sum=diag(AllCombs_mat*S*AllCombs_mat'); % sum of weights in each cluster
        end
 
        %% Calculation of the structural likelihood L_s
        L_s=C_sum/N_nodes;                        % structural likelihood
         
        %% Calculation of the temporal likelihood L_t
                 
        L_t=zeros(size(L_s));
        ICI_seq=diff(Locs(AllCombs),[],2);

        %% Apply PCA projection to extract the rhthmic features from the observed ICI measurements of clicks in each candidate cluster
       
        if Rank==3
            p=2;
            observation=ICI_seq*EigV{Rank};
            NOC=size(Params(Rank).m,2);
        else
            p=3;
            observation=ICI_seq*EigV_3D{Rank};
            NOC=size(Params3D_rev(Rank).m,2);
        end
        
        %% Evaluate the temporal likelihood using MGGD (Multivariate-Generalized-Gaussian-Distribution) model
            
        for iter=1:length(G_ind)
            for q=1:NOC
                if Rank==3
                    L_t(G_ind(iter))=L_t(G_ind(iter))+MGGD_fun(observation(G_ind(iter),:),Params(Rank).Sigma{q},Params(Rank).Beta(q),Params(Rank).mu(q,:),Params(Rank).m(q),p)/Params(Rank).Z_max(q);
                else
                    L_t(G_ind(iter))=L_t(G_ind(iter))+MGGD_fun(observation(G_ind(iter),:),Params3D_rev(Rank).Sigma{q},Params3D_rev(Rank).Beta(q),Params3D_rev(Rank).mu(q,:),Params3D_rev(Rank).m(q),p)/Params3D_rev(Rank).Z_max(q);
                end
            end
        end
      
%%  calculation of the utility function

        U=L_s(G_ind)-exp(-alpha1*L_t(G_ind))-alpha2*exp(-Rank);
        
        Chosen_clusters=AllCombs(G_ind,:);
        Chosen_clusters_indexed=mat2inds(Chosen_clusters,length(Locs));       
        All_clusters=[All_clusters ; Chosen_clusters_indexed];
        All_U=[All_U ; U];
    end

 %% Pick the cluster with the highest utility score

    L_hat(k)={[Locs(All_clusters(All_U==max(All_U),:)==1) ; Pks(All_clusters(All_U==max(All_U),:)==1)]};
    L_hat_inds(k)={find(All_clusters(All_U==max(All_U),:)==1)};

    if ~isempty(All_U)
        %% Discard the previous cluster to meet the ortogonality constraint
        U_max(k)=max(All_U);
        Locs(All_clusters(All_U==max(All_U),:)==1)=[];
        Pks(All_clusters(All_U==max(All_U),:)==1)=[];        
        IPI(All_clusters(All_U==max(All_U),:)==1)=[];
        S(All_clusters(All_U==max(All_U),:)==1,:)=[];
        S(:,All_clusters(All_U==max(All_U),:)==1)=[];

        L(k+1)=length(Locs);

        if L(k+1)==L(k) || L(k+1)<4
            break;
        end   
    else
        break;
    end
end


  
        
end
