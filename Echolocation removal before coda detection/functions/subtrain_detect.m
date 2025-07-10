function Detected_subtrains=subtrain_detect(roi_max,ref_ToAs,F_weights,locs,buffer,F_ds,All_objs,mode,Tag_flag)

% locs=locs_tmp
% ref_ToAs=ref_ToAs_tmp;
%             locs=locs_left
%             ref_ToAs=ref_ToAs_left

            %% Feature extraction
            if size(buffer,2)==1
                buffer=buffer';
            end

            W=3e-3; roi=[5e-3 roi_max]; %roi=[5e-3 40e-3]; %roi=[5e-3 18e-3];
            Click=[]; Click_multi=[]; EnV=[]; Noise=[];
            Click_multi_L1=length(buffer(int32((locs(1)-roi(1))*F_ds):int32((locs(1)+roi(2))*F_ds)));
                for iter=1:length(locs) 
                    if iter>1
                        L1=length(buffer(int32((locs(iter)-W)*F_ds):int32((locs(iter)+W)*F_ds)));
                        L2=length(Click(iter-1,:));
                        Click_TMP=buffer(int32((locs(iter)-W)*F_ds):int32((locs(iter)+W)*F_ds));
                        if L1<L2
                            Click(iter,:)=zeros(1,L2);
                            Click(iter,1:L1)=Click_TMP;
                        elseif L1>L2
                            Click(iter,:)=zeros(1,L2);
                            Click(iter,:)=Click_TMP(1:L2);                            
                        else
                            Click(iter,:)=Click_TMP;
                        end
                    else
                        Click(iter,:)=buffer(int32((locs(iter)-W)*F_ds):int32((locs(iter)+W)*F_ds));
                    end
                     Click_multi_tmp=buffer(int32((locs(iter)-roi(1))*F_ds):int32((locs(iter)+roi(2))*F_ds));
                     if Click_multi_L1==length(Click_multi_tmp)
                         Click_multi(iter,:)=Click_multi_tmp;
                     elseif Click_multi_L1<length(Click_multi_tmp)
                         Click_multi(iter,:)=Click_multi_tmp(1:Click_multi_L1);
                     elseif Click_multi_L1>length(Click_multi_tmp)
                         Click_multi_tmp2=zeros(1,Click_multi_L1);
                         Click_multi_tmp2(1:length(Click_multi_tmp))=Click_multi_tmp;
                         Click_multi(iter,:)=Click_multi_tmp2;
                     end

                     % [E_up,~] = envelope(Click(iter,:)); EnV(iter,:)=E_up;
                     % Noise(iter,:)=[buffer(int32((locs(iter)-2*W)*F_ds):int32((locs(iter)-W)*F_ds))' buffer(int32((locs(iter)+roi(2))*F_ds):int32((locs(iter)+roi(2)+W)*F_ds))']; 
                     Slant_Delay(iter)=estimate_slant_delays(Click_multi(iter,:),F_ds);
                end
                % Noise(:,1)=[];
        
            %%
            % Pd_clicks=wavelet_click_classifier(buffer,locs,F_ds,trainedModel)+eps;
            ICI_min=0.4; %0.34;
            ICI_max=2.5; %2.6;
            Nd=length(locs);
            L_tot = 1e9*ones(Nd,Nd);
            P_transition=zeros(Nd,Nd);
            ICI=zeros(Nd,Nd);
            for i=1:Nd-1
               for j=i+1:Nd               
                        if locs(j)-locs(i)>ICI_min && locs(j)-locs(i)<ICI_max
                           ICI(i,j)= locs(j)-locs(i);
                           % [~,~,Spatial,Orientation,~,Temporal]=extract_click_pair_features(Click(i,:),Click(j,:),Click_multi(i,:),Click_multi(j,:),Noise(i,:),Noise(j,:),EnV(i,:),EnV(j,:),F_ds);
                           % [Spatial,Orientation,Temporal]=extract_click_pair_features_V2(Click(i,:),Click(j,:),Click_multi(i,:),Click_multi(j,:),F_ds);
                            [~,~,~,Membership]=extract_click_pair_features_V3(F_weights,ICI(i,j),Click(i,:),Click(j,:),Click_multi(i,:),Click_multi(j,:),F_ds,All_objs,Tag_flag); 
                            P_transition(i,j)=Membership;
                            L_tot(i,j)=-log(Membership);%-log(1+lone_p); %-log(Pd_clicks(i,1))-log(Pd_clicks(j,1));
                            % L_tot(i,j)=-Estimate_log_likelihood_V2(Temporal,Orientation,Spatial,Params_Orientation,Params_Spatial);
                        
                        end          
               end
            end

            if mode
            Pfd=eps*zeros(1,Nd);
            for i=1:Nd
                Pi=max(P_transition(i,:));
                Pj=max(P_transition(:,i));
                Pfd(i)=max([Pi Pj]);
            end

            for i=1:Nd-1
               for j=i+1:Nd
                    if locs(j)-locs(i)>ICI_min && locs(j)-locs(i)<ICI_max
                       L_tot(i,j)=L_tot(i,j)-log(Pfd(i));%+log(Pfd(j))); 
                    end
               end
            end
            % 
             [C,cost] = do_assignment (L_tot,-log(1-Pfd), 1e9);
            else
                lone_p=-0.07;
                [C,cost] = do_assignment (L_tot,-log(-lone_p), 1e9);
            end
            %  % [C,cost] = do_assignment (L_tot,Pd_clicks, 1e9);
            % % % lone_p=-1;
            % % % solve LA problem
            % % % cd('C:\Users\User\Desktop\source separation\code\Benchmarks')
            % [C,cost] = do_assignment (L_tot,-log(-lone_p), 1e9);
            % % % cd('C:\Users\User\Desktop\source separation\Annotaion_GUI\Click_level_features')
            % % % unravel the chains. When finished, nchain will be the number of chains and chain_nr(i) will
            % % % be the assignment of click i to a chain
            % % 
            nchain = 0;
            chain_nr = zeros(Nd,1);
            for i = 1:Nd
                j = find(C(i,1:Nd));
                if (~isempty(j))
                    if (chain_nr(i) && ~any(chain_nr(j)))
                        chain_nr(j) = chain_nr(i);
                    else
                        if(chain_nr(i)==0 && ~any(chain_nr(j)))
                            nchain = nchain + 1;
                            chain_nr(j) = nchain;
                            chain_nr(i) = nchain;
                        end
                    end
                end
            end
        
             % figure;
             % chain_nr=chain_nr+1;

             El=[]; c=0;
                for i=1:nchain+1
                     ind=find(chain_nr==i);          
                     L_S=0;
                     for x=1:length(ind)-1
                         L_S=L_S+L_tot(ind(x),ind(x+1));
                         tmp(x)=L_tot(ind(x),ind(x+1));
                     end
                     if min(diff(locs(ind)))<ICI_min
                         El=[El ind'];
                     end
                     if median(ref_ToAs(ind))<3 && ~Tag_flag
                        El=[El ind'];
                     end
                     if std(rmoutliers(ref_ToAs(ind)))>1.5 && ~Tag_flag
                        El=[El ind'];
                     end                     
                     ICI_ind=diff(locs(ind));
                     if length(ind)==3 & abs(log(ICI_ind(2)/ICI_ind(1)))>0.1
                        El=[El ind'];
                     end
                     % if length(ind)==3 & mean(ICI_ind)<0.7
                     %    El=[El ind'];
                     % end
                     if length(ind)<2 %& min(diff(locs(ind)))<1.5
                         El=[El ind'];
                     elseif isscalar(ind)
                         El=[El ind'];
                     else
                         c=c+1;
                         Detection_inds(c)={ind};
                     end        
                end
                El=unique(El);
                All_traces=chain_nr; All_traces(El)=[];
                All_traces(find(All_traces==0))=[];
                Detected_traces=unique(All_traces)';
                 Gr={};
                 if ~isempty(Detected_traces)
                    for i=Detected_traces
                            Gr(i)={find(chain_nr==i)};
                    end
                    
                    Detected_subtrains= Gr(~cellfun(@isempty, Gr));
                 else
                     Detected_subtrains={};
                 end


                    % subplot(2,1,1);
                    % for i=1:length(Detected_subtrains)
                    %     ind=Detected_subtrains{i};              
                    %     hold on; plot(t1+locs(ind),pks(ind),'*','LineWidth',1.5);
                    %     legendInfo{i}=num2str(i);
                    % end
                    % 
                    % subplot(2,1,2); 
                    % for i=1:length(Detected_subtrains)
                    %     ind=Detected_subtrains{i};              
                    %     plot(t1+locs(ind),ref_ToAs(ind),'*','LineWidth',2); hold on; grid on;
                    %     legendInfo{i}=num2str(i); 
                    % end

end