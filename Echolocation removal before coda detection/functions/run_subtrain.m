   
function [Detections,subtarin_evaluation_results,tp_class,fp_class,fn_class,N_GT,TP_detector,TP_detector_Baggenstos,FP_analysis]=run_subtrain(trainedModel,roi,F_weights,Buffer_length,ToAs,Refs,test,lone_p,F_ds,All_objs,snr_lim,Probs)

        subtarin_evaluation_results=[];
        FP_analysis=[];
        tp_class=[];
        fp_class=[];
        fn_class=[];
        N_GT=[];
        TP_detector=[];
        TP_detector_Baggenstos=[];

Eval_flag=1;
    Noise_flag=0;
    Pre_detector_flag=1;
    FP_clicks=0;
    FP_clicks_Baggenstos=0;
    TP_detector=zeros(1,6);
    TP_detector_Baggenstos=zeros(1,6);
    subtarin_evaluation_results=zeros(3,6*2);
    TP=zeros(1,6); TP_Baggenstos=zeros(1,6); 
    FP=zeros(1,6); FP_Baggenstos=zeros(1,6); 
    IDSw=zeros(1,6); IDSw_Baggenstos=zeros(1,6); 
    MD=zeros(1,6); MD_Baggenstos=zeros(1,6); 
    N_GT=zeros(1,6);
    for NOW_buffer=1:6
       tp_class(NOW_buffer).vec=[]; tp_class(NOW_buffer).vec_Baggenstos=[]; tp_class(NOW_buffer).vec_Spectral=[]; tp_class(NOW_buffer).vec_ML=[];
       fp_class(NOW_buffer).vec=[]; fp_class(NOW_buffer).vec_Baggenstos=[]; fp_class(NOW_buffer).vec_Spectral=[]; fp_class(NOW_buffer).vec_ML=[];
       fn_class(NOW_buffer).vec=[]; fn_class(NOW_buffer).vec_Baggenstos=[]; fn_class(NOW_buffer).vec_Spectral=[]; fn_class(NOW_buffer).vec_ML=[];
    end
    
    
    if length(ToAs)==6
        ToA_1=ToAs{1}; ToA_2=ToAs{2}; ToA_3=ToAs{3}; ToA_4=ToAs{4}; ToA_5=ToAs{5}; ToA_6=ToAs{6};
        ref_1=Refs{1}; ref_2=Refs{2}; ref_3=Refs{3}; ref_4=Refs{4}; ref_5=Refs{5}; ref_6=Refs{6};
    elseif length(ToAs)==5
        ToA_1=ToAs{1}; ToA_2=ToAs{2}; ToA_3=ToAs{3}; ToA_4=ToAs{4}; ToA_5=ToAs{5};
        ref_1=Refs{1}; ref_2=Refs{2}; ref_3=Refs{3}; ref_4=Refs{4}; ref_5=Refs{5};
    elseif length(ToAs)==4
        ToA_1=ToAs{1}; ToA_2=ToAs{2}; ToA_3=ToAs{3}; ToA_4=ToAs{4};
        ref_1=Refs{1}; ref_2=Refs{2}; ref_3=Refs{3}; ref_4=Refs{4};
    elseif length(ToAs)==3
        ToA_1=ToAs{1}; ToA_2=ToAs{2}; ToA_3=ToAs{3};
        ref_1=Refs{1}; ref_2=Refs{2}; ref_3=Refs{3}; 
    elseif length(ToAs)==2
        ToA_1=ToAs{1}; ToA_2=ToAs{2}; 
        ref_1=Refs{1}; ref_2=Refs{2};
    elseif isscalar(ToAs)
        ToA_1=ToAs{1}; 
        ref_1=Refs{1}; 
    end


    % roi=15e-3;
    t_index=0;
    t1=-Buffer_length;
    t2=0;    
    TP_bench=0;
    Total_clicks=0;    
    plot_flag=0;
    Thresh=4;

    while(1)       
        
        t1=t1+Buffer_length;
        t2=t2+Buffer_length;
        t_index=t_index+1;
        t1
        % if t2>333
        %     break;
        % end
        if length(ToAs)==6
            GT_1=ToA_1(ToA_1>t1 & ToA_1<t2); GT_1_ref=ref_1(ToA_1>t1 & ToA_1<t2);
            GT_2=ToA_2(ToA_2>t1 & ToA_2<t2); GT_2_ref=ref_2(ToA_2>t1 & ToA_2<t2);
            GT_3=ToA_3(ToA_3>t1 & ToA_3<t2); GT_3_ref=ref_3(ToA_3>t1 & ToA_3<t2);
            GT_4=ToA_4(ToA_4>t1 & ToA_4<t2); GT_4_ref=ref_4(ToA_4>t1 & ToA_4<t2); 
            GT_5=ToA_5(ToA_5>t1 & ToA_5<t2); GT_5_ref=ref_5(ToA_5>t1 & ToA_5<t2);
            GT_6=ToA_6(ToA_6>t1 & ToA_6<t2); GT_6_ref=ref_6(ToA_6>t1 & ToA_6<t2);
        elseif length(ToAs)==5
            GT_1=ToA_1(ToA_1>t1 & ToA_1<t2); GT_1_ref=ref_1(ToA_1>t1 & ToA_1<t2);
            GT_2=ToA_2(ToA_2>t1 & ToA_2<t2); GT_2_ref=ref_2(ToA_2>t1 & ToA_2<t2);
            GT_3=ToA_3(ToA_3>t1 & ToA_3<t2); GT_3_ref=ref_3(ToA_3>t1 & ToA_3<t2);
            GT_4=ToA_4(ToA_4>t1 & ToA_4<t2); GT_4_ref=ref_4(ToA_4>t1 & ToA_4<t2); 
            GT_5=ToA_5(ToA_5>t1 & ToA_5<t2); GT_5_ref=ref_5(ToA_5>t1 & ToA_5<t2); 
            GT_6=[]; GT_6_ref=[];
        elseif length(ToAs)==4
            GT_1=ToA_1(ToA_1>t1 & ToA_1<t2); GT_1_ref=ref_1(ToA_1>t1 & ToA_1<t2);
            GT_2=ToA_2(ToA_2>t1 & ToA_2<t2); GT_2_ref=ref_2(ToA_2>t1 & ToA_2<t2);
            GT_3=ToA_3(ToA_3>t1 & ToA_3<t2); GT_3_ref=ref_3(ToA_3>t1 & ToA_3<t2);
            GT_4=ToA_4(ToA_4>t1 & ToA_4<t2); GT_4_ref=ref_4(ToA_4>t1 & ToA_4<t2);  
            GT_5=[]; GT_5_ref=[];
            GT_6=[]; GT_6_ref=[];
        elseif length(ToAs)==3
            GT_1=ToA_1(ToA_1>t1 & ToA_1<t2); GT_1_ref=ref_1(ToA_1>t1 & ToA_1<t2);
            GT_2=ToA_2(ToA_2>t1 & ToA_2<t2); GT_2_ref=ref_2(ToA_2>t1 & ToA_2<t2);
            GT_3=ToA_3(ToA_3>t1 & ToA_3<t2); GT_3_ref=ref_3(ToA_3>t1 & ToA_3<t2);
            GT_4=[]; GT_4_ref=[];
            GT_5=[]; GT_5_ref=[];
            GT_6=[]; GT_6_ref=[];
        elseif length(ToAs)==2
            GT_1=ToA_1(ToA_1>t1 & ToA_1<t2); GT_1_ref=ref_1(ToA_1>t1 & ToA_1<t2);
            GT_2=ToA_2(ToA_2>t1 & ToA_2<t2); GT_2_ref=ref_2(ToA_2>t1 & ToA_2<t2);
            GT_3=[]; GT_3_ref=[];
            GT_4=[]; GT_4_ref=[];
            GT_5=[]; GT_5_ref=[];
            GT_6=[]; GT_6_ref=[];
        elseif isscalar(ToAs)
            GT_1=ToA_1(ToA_1>t1 & ToA_1<t2); GT_1_ref=ref_1(ToA_1>t1 & ToA_1<t2);
            GT_2=[]; GT_2_ref=[];
            GT_3=[]; GT_3_ref=[];
            GT_4=[]; GT_4_ref=[];
            GT_5=[]; GT_5_ref=[];
            GT_6=[]; GT_6_ref=[];
        end

        GT_cell={GT_1,GT_2,GT_3, GT_4, GT_5, GT_6};
        GT_cell=GT_cell(~cellfun(@isempty, GT_cell));
        GT=[GT_1 GT_2 GT_3 GT_4 GT_5 GT_6]; 

       

    % GT={GT_1,GT_2,GT_3,GT_4,GT_5};  
    % GT = GT(~cellfun(@isempty, GT));
    % NOW=0;
    % for X=1:length(GT)
    %     Consi_x=[];
    %     ICI_X=diff(GT{X});
    %     for X2=1:length(ICI_X)-1
    %         Consi_x(X2)=abs(log(ICI_X(X2+1)/ICI_X(X2)));
    %     end
    %     if length(GT{X})>2 & min(Consi_x)<0.08
    %        NOW=NOW+1;
    %     end
    % end
    
        if t2*F_ds>length(test)
            break;
        end
        buffer=test(t1*F_ds+1:t2*F_ds);
        t_buffer=t1+[0:1/F_ds:(1/F_ds)*(length(buffer)-1)];
        tb=[0:1/F_ds:(1/F_ds)*(length(buffer)-1)];
        % SNR_min=10;
        if plot_flag
            figure; subplot(2,1,1); plot(t_buffer,buffer); hold on; grid on;    
            plot(GT_1,zeros(1,length(GT_1)),'*','LineWidth',2);
            plot(GT_2,zeros(1,length(GT_2)),'*','LineWidth',2);
            plot(GT_3,zeros(1,length(GT_3)),'*','LineWidth',2);
            plot(GT_4,zeros(1,length(GT_4)),'*','LineWidth',2);
            plot(GT_5,zeros(1,length(GT_5)),'*','LineWidth',2);
            plot(GT_6,zeros(1,length(GT_6)),'*','LineWidth',2);
        end

        if isempty(GT) 
             Detected_subtrains={};
             Detected_subtrains_bench={};
        else
            if Pre_detector_flag 
                locs=[GT_1 GT_2 GT_3 GT_4 GT_5 GT_6]-t1;
                locs_ML=locs;
                % El_inds=1.1*roi;
                El_inds=25e-3;
                locs(locs<El_inds | locs>(Buffer_length-El_inds))=[];
                if ~isempty(locs)
                    pks=Peaks_extract(buffer,locs,F_ds);
                else
                    pks=[];
                end
            else
                [pks,locs] =findpeaks(buffer,F_ds,'MinPeakDistance',4.5e-3,'MinPeakHeight',snr_lim);
                
                % [pks,locs] =findpeaks(buffer,F_ds,'MinPeakDistance',4e-3,'MinPeakHeight',2e-3);
                % [pks,locs] =findpeaks(buffer,F_ds,'MinPeakDistance',4e-3,'MinPeakHeight',0.015);
        
               El_inds=1.1*roi;
               pks(locs<El_inds | locs>(Buffer_length-El_inds))=[];
               locs(locs<El_inds | locs>(Buffer_length-El_inds))=[];
           
                pulse_len=1e-3;
                Eliminate=[]; SNR2=[];
                % SNR_min=10; SNR=[];
                SNR_min=2.3;
                for i=1:length(locs)
                    % l_cand=locs(i+1:end)'-locs(i);
                    click=buffer((locs(i)-pulse_len)*F_ds:(locs(i)+pulse_len)*F_ds);
                    noise1=[buffer((locs(i)-2*pulse_len)*F_ds:(locs(i)-pulse_len)*F_ds) ; buffer((locs(i)+pulse_len)*F_ds:(locs(i)+2*pulse_len)*F_ds)];
                    noise2=buffer((locs(i)-3*pulse_len)*F_ds:(locs(i)-pulse_len)*F_ds);
                    noise3=buffer((locs(i)+pulse_len)*F_ds:(locs(i)+3*pulse_len)*F_ds);
                    [~,min_noise]=min([median(abs(noise1)) median(abs(noise2)) median(abs(noise3))]);
                    Noise_ops={noise1,noise2,noise3};
                    noise=Noise_ops{min_noise};
                    [minlen,idx_min]=min([length(click) length(noise)]);
                    if idx_min==1
                        noise=noise(1:minlen);
                    elseif idx_min==2
                        click=click(1:minlen);
                    end                
                    SNR(i)=snr(click,noise);
                    SNR2(i)=log(max(abs(click))/median(abs(noise)));
                    % MP(i)=length(find(l_cand<15e-3));
                    % if MP(i)==1
                    %     Eliminate=[Eliminate i+1];
                    % end                       
                end
                pks(SNR2<SNR_min)=[];
                locs(SNR2<SNR_min)=[];
                % pks(Eliminate)=[];
                % locs(Eliminate)=[];            
                % Pk_max=80;
                % if length(locs)>Pk_max
                %     [pks,I]=sort(pks,'descend');
                %     locs=locs(I);
                %     pks=pks(1:Pk_max);
                %     locs=locs(1:Pk_max);
                %     [locs,I]=sort(locs);
                %     pks=pks(I);
                % end
            end
                if length(locs)>3 & ~Pre_detector_flag
                    [ref_ToAs,locs,pks]=remove_multipass(Thresh,F_ds,t_buffer,buffer,locs,pks,t1,roi);
                else
                    [ref_ToAs,~]=determine_reflection_ToAs(buffer,locs,0,roi,F_ds,0);
                end
                if plot_flag
                    plot(t1+locs,pks,'kx','LineWidth',2);
                    subplot(2,1,2); hold off;
                end
                % roi=15e-3; 
                cd('C:\Users\User\Desktop\source separation\Annotaion_GUI\functions')
                % [ref_ToAs,~]=determine_reflection_ToAs(buffer,locs,0,roi,F_ds,0);
                if plot_flag
                    plot(t1+locs,ref_ToAs,'kx','LineWidth',2);
                    grid on; hold on;
                    plot(GT_1',GT_1_ref,'o','LineWidth',2);
                    plot(GT_2,GT_2_ref,'o','LineWidth',2);
                    plot(GT_3,GT_3_ref,'o','LineWidth',2);
                    plot(GT_4,GT_4_ref,'o','LineWidth',2);
                    plot(GT_5,GT_5_ref,'o','LineWidth',2);
                    plot(GT_6,GT_6_ref,'o','LineWidth',2);
                end
                mode=0;
                if length(locs)>2
                    Detected_subtrains=subtrain_detect(trainedModel,ref_ToAs,F_weights,locs,buffer,lone_p,F_ds,All_objs,mode);
                    Detected_subtrains=Merge_chains(Detected_subtrains,locs);
                    if Eval_flag
                        Detected_subtrains_bench=Baggenstoss_likelihood(locs,buffer,2.5,F_ds,Probs);
                        Detected_subtrains_Spectral=Spectral_cluster(locs,buffer,F_ds,GT_cell);
                        ICI_Mk_max=0;
                        ICI_Mk_min=0;
                        if length(GT_cell)>1
                            Mk=[]; ICI_mu=[];
                            ICI_sigma=[]; Mk_flag=0;
                            for Mk_ind=1:length(GT_cell)
                                Mk=[Mk length(GT_cell{Mk_ind})];
                                if min(Mk)>1
                                    ICI_Mk_ind=diff(GT_cell{Mk_ind});
                                    ICI_Mk_max(Mk_ind)=max(ICI_Mk_ind);
                                    ICI_Mk_min(Mk_ind)=min(ICI_Mk_ind);
                                    ICI_mu=[ICI_mu mean(ICI_Mk_ind)];
                                    std_insert=std(ICI_Mk_ind);
                                    if std_insert==0
                                        std_insert=0.01;
                                    end
                                    ICI_sigma=[ICI_sigma std_insert];
                                else
                                    Mk_flag=1;
                                end
                            end
                            ICI_max_el=max(ICI_Mk_max);
                            ICI_min_el=min(ICI_Mk_min);
                            if  ~Mk_flag
                                Detected_subtrains_ML=Timing_only_ML(locs_ML,ICI_mu,ICI_sigma,ICI_min_el,ICI_max_el,length(GT_cell),Mk);
                            else
                                Detected_subtrains_ML={1:length(cell2mat(GT_cell))};
                            end
                        else
                           Detected_subtrains_ML={1:length(GT_cell{1})};
                        end
                    end
                else
                    Detected_subtrains={};
                    Detected_subtrains_bench={};
                    Detected_subtrains_Spectral={};
                    Detected_subtrains_ML={};
                end
        end


            
% if length(locs)>3            
%     % tic
%     % Detected_subtrains=Global_converge_buffer(NOW,t1,buffer,locs,F_ds,ref_ToAs,F_weights,All_objs);
%     % Detected_subtrains=Merge_chains(Detected_subtrains,locs);
%     % toc
% else
%     Detected_subtrains={};
% end
            % locs_left=locs(~ismember(1:length(locs),cell2mat(Detected_subtrains)));
            % ref_ToAs_left=ref_ToAs(~ismember(1:length(locs),cell2mat(Detected_subtrains)));
            % Detected_subtrains_left=subtrain_detect(ref_ToAs_left,F_weights,locs_left,buffer,lone_p,F_ds,All_objs)
            % 


            % if length(Detected_subtrains)>1
            %     El_inds=Eliminate_ref(Detected_subtrains,locs);
            %     if ~isempty(El_inds)
            %         ref_ToAs(El_inds)=[];
            %         pks(El_inds)=[];
            %         locs(El_inds)=[];
            %        Detected_subtrains=subtrain_detect(ref_ToAs,F_weights,locs,buffer,lone_p,F_ds,All_objs);
            %     end
            % end



            % for K_sub=1:length(Detected_subtrains)
            %     tmp=sort(Detected_subtrains{K_sub}');
            %     chain=locs(tmp)';
            % 
            %     tmp(find(diff(chain)<30e-3)+1)=[];
            %     Detected_subtrains(K_sub)={tmp'};
            % 
            % end

                % 
                % if ~isempty(Detected_subtrains)
                %     % Detected_subtrains=Merge_subtrains(Detected_traces,chain_nr,L_tot,locs);
                %     cd('C:\Users\User\Desktop\source separation\Annotaion_GUI\Click_level_features')
                % 
                %     if length(GT_1)>2
                %         Chosen_score_1=Eval_score(GT_1,Detected_subtrains,locs,t1);
                %         Chosen_score_1_bench=Eval_score(GT_1,Detected_subtrains_bench,locs,t1);                        
                %     else
                %         Chosen_score_1=Eval_score([],Detected_subtrains,locs,t1);
                %         Chosen_score_1_bench=Eval_score([],Detected_subtrains_bench,locs,t1);
                %     end
                %     if length(GT_2)>2
                %         Chosen_score_2=Eval_score(GT_2,Detected_subtrains,locs,t1);
                %         Chosen_score_2_bench=Eval_score(GT_2,Detected_subtrains_bench,locs,t1);
                %     else
                %         Chosen_score_2=Eval_score([],Detected_subtrains,locs,t1);
                %         Chosen_score_2_bench=Eval_score([],Detected_subtrains_bench,locs,t1);
                %     end
                %     if length(GT_3)>2
                %         Chosen_score_3=Eval_score(GT_3,Detected_subtrains,locs,t1);
                %         Chosen_score_3_bench=Eval_score(GT_3,Detected_subtrains_bench,locs,t1);
                %     else
                %         Chosen_score_3=Eval_score([],Detected_subtrains,locs,t1);
                %         Chosen_score_3_bench=Eval_score([],Detected_subtrains_bench,locs,t1);
                %     end
                %     if length(GT_4)>2
                %         Chosen_score_4=Eval_score(GT_4,Detected_subtrains,locs,t1);
                %         Chosen_score_4_bench=Eval_score(GT_4,Detected_subtrains_bench,locs,t1);
                %     else
                %         Chosen_score_4=Eval_score([],Detected_subtrains,locs,t1);
                %         Chosen_score_4_bench=Eval_score([],Detected_subtrains_bench,locs,t1);
                %     end
                %     if length(GT_5)>2
                %         Chosen_score_5=Eval_score(GT_5,Detected_subtrains,locs,t1);
                %         Chosen_score_5_bench=Eval_score(GT_5,Detected_subtrains_bench,locs,t1);
                %     else
                %         Chosen_score_5=Eval_score([],Detected_subtrains,locs,t1);
                %         Chosen_score_5_bench=Eval_score([],Detected_subtrains_bench,locs,t1);
                %     end
                %     if length(GT_6)>2
                %         Chosen_score_6=Eval_score(GT_6,Detected_subtrains,locs,t1);
                %         Chosen_score_6_bench=Eval_score(GT_6,Detected_subtrains_bench,locs,t1);
                %     else
                %         Chosen_score_6=Eval_score([],Detected_subtrains,locs,t1);
                %         Chosen_score_6_bench=Eval_score([],Detected_subtrains_bench,locs,t1);
                %     end
                % 
                %     Total_score=Chosen_score_1+Chosen_score_2+Chosen_score_3+Chosen_score_4+Chosen_score_5+Chosen_score_6;
                %     Total_score_bench=Chosen_score_1_bench+Chosen_score_2_bench+Chosen_score_3_bench+Chosen_score_4_bench+Chosen_score_5_bench+Chosen_score_6_bench;                    
                %     % Pd=Total_score(1)/(length([GT_1 GT_2 GT_3 GT_4]));
                %     TP=TP+Total_score(1);
                %     TP_bench=TP_bench+Total_score_bench(1);
                %     Total_clicks=Total_clicks+length([GT_1 GT_2 GT_3 GT_4 GT_5 GT_6]);
                % end
             if Eval_flag
                if ~isempty(Detected_subtrains)
                    if Noise_flag
                         %% Detection performance- FAR analysis
                         if ~isempty(Detected_subtrains)
                             FP_clicks=FP_clicks+length(cell2mat(Detected_subtrains(:)));
                         end
                         if ~isempty(Detected_subtrains_bench)
                              FP_clicks_Baggenstos=FP_clicks_Baggenstos+length(cell2mat(Detected_subtrains_bench(:)));
                         end
                    else
                        %% Detection performance- Pd analysis
                        GT_clicks=cell2mat(GT_cell);
                        % Detected_clicks=t1+locs(cell2mat(Detected_subtrains));
                        % Detected_clicks_Baggenstos=t1+locs(cell2mat(Detected_subtrains_bench));
                        Detected_clicks=t1+locs(cell2mat(Detected_subtrains(:)));
                        Detected_clicks_Baggenstos=t1+locs(cell2mat(Detected_subtrains_bench(:)));
    
                        TP_detector(length(GT_cell))=TP_detector(length(GT_cell))+sum(ismember(GT_clicks,Detected_clicks));
                        TP_detector_Baggenstos(length(GT_cell))=TP_detector_Baggenstos(length(GT_cell))+sum(ismember(GT_clicks,Detected_clicks_Baggenstos));
    
                        %% Total performance- MOTA analysis
                        Eval_results=Eval_buffer(GT_cell,Detected_subtrains,locs,t1);
                        Eval_results_Baggenstos=Eval_buffer(GT_cell,Detected_subtrains_bench,locs,t1);
                      
                        NOW_buffer=Eval_results(1);
                        TP(NOW_buffer)=TP(NOW_buffer)+Eval_results(2); TP_Baggenstos(NOW_buffer)=TP_Baggenstos(NOW_buffer)+Eval_results_Baggenstos(2);
                        IDSw(NOW_buffer)=IDSw(NOW_buffer)+Eval_results(3); IDSw_Baggenstos(NOW_buffer)=IDSw_Baggenstos(NOW_buffer)+Eval_results_Baggenstos(3);                  
                        MD(NOW_buffer)=MD(NOW_buffer)+Eval_results(4); MD_Baggenstos(NOW_buffer)=MD_Baggenstos(NOW_buffer)+Eval_results_Baggenstos(4);
                        FP(NOW_buffer)=FP(NOW_buffer)+Eval_results(6); FP_Baggenstos(NOW_buffer)=FP_Baggenstos(NOW_buffer)+Eval_results_Baggenstos(6);
                        N_GT(NOW_buffer)=N_GT(NOW_buffer)+Eval_results(5);
    
                        %% Separation performance- F-score analysis
                        [tp,fp,fn]=Eval_buffer_calssification(GT_cell,Detected_subtrains,locs,t1);
                        [tp_Baggenstos,fp_Baggenstos,fn_Baggenstos]=Eval_buffer_calssification(GT_cell,Detected_subtrains_bench,locs,t1);
                        [tp_Spectral,fp_Spectral,fn_Spectral]=Eval_buffer_calssification(GT_cell,Detected_subtrains_Spectral,locs,t1);                       
                        [tp_ML,fp_ML,fn_ML]=Eval_buffer_calssification(GT_cell,Detected_subtrains_ML,locs_ML,t1);
    
                        tp_class(NOW_buffer).vec=[tp_class(NOW_buffer).vec tp]; tp_class(NOW_buffer).vec_Baggenstos=[tp_class(NOW_buffer).vec_Baggenstos tp_Baggenstos]; tp_class(NOW_buffer).vec_Spectral=[tp_class(NOW_buffer).vec_Spectral tp_Spectral]; tp_class(NOW_buffer).vec_ML=[tp_class(NOW_buffer).vec_ML tp_ML];
                        fp_class(NOW_buffer).vec=[fp_class(NOW_buffer).vec fp]; fp_class(NOW_buffer).vec_Baggenstos=[fp_class(NOW_buffer).vec_Baggenstos fp_Baggenstos]; fp_class(NOW_buffer).vec_Spectral=[fp_class(NOW_buffer).vec_Spectral fp_Spectral]; fp_class(NOW_buffer).vec_ML=[fp_class(NOW_buffer).vec_ML fp_ML];
                        fn_class(NOW_buffer).vec=[fn_class(NOW_buffer).vec fn]; fn_class(NOW_buffer).vec_Baggenstos=[fn_class(NOW_buffer).vec_Baggenstos fn_Baggenstos]; fn_class(NOW_buffer).vec_Spectral=[fn_class(NOW_buffer).vec_Spectral fn_Spectral]; fn_class(NOW_buffer).vec_ML=[fn_class(NOW_buffer).vec_ML fn_ML];   
                    end
                end
             end
                subplot(2,2,2);
                % dotH=plot(t1+locs(ind),ref_ToAs(ind),'*','LineWidth',2); hold on;
                for i=1:length(Detected_subtrains)
                    ind=Detected_subtrains{i};              
                    dotH=plot(t1+locs(ind),ref_ToAs(ind),'*','LineWidth',2); hold on; 
                       pause(0.05);  % calls DRAWNOW implicitly
                       set(dotH, 'XData', t1+locs(ind), 'YData', ref_ToAs(ind)); hold on;
                end
                xlabel('Time [sec]'); ylabel('Slant delay [ms]'); 
                grid on; title('Click subtrain detector');

                colors={'rx','gx','cx','mx','yx','bx','r*','g*','c*','m*','y*','b*'};
                if plot_flag
                    % figure;
                    subplot(2,1,1);
                    for i=1:length(Detected_subtrains)
                        ind=Detected_subtrains{i};              
                        hold on; plot(t1+locs(ind),pks(ind),colors{i},'LineWidth',1.5);
                        legendInfo{i}=num2str(i);
                    end

                    subplot(2,1,2); %plot([],[],'x');
                    figure;
                    for i=1:length(Detected_subtrains)
                        ind=Detected_subtrains{i};              
                        plot(t1+locs(ind),ref_ToAs(ind),colors{i},'LineWidth',2); hold on; grid on;
                        legendInfo{i}=num2str(i); 
                    end
                    legend(legendInfo);
                    title(num2str(round(Pd,2)))
                    % pause;
                end

                Confidence=[];
                for i=1:length(Detected_subtrains)
                    ind=Detected_subtrains{i};              
                    ICI_ind=diff(sort(locs(ind)'));
                    Consi=[];
                    for j=1:length(ICI_ind)-1
                        Consi(j)=log(ICI_ind(j+1)/ICI_ind(j));
                    end
                    Confidence(i)=sum(Consi<0.156)/length(ind);
                    % if Confidence(i)<0.5
                    %     Confidence(i)=eps;
                    % elseif Confidence(i)>0.5
                    %     Confidence(i)=1;
                    % end
                end
            
                cd('C:\Users\User\Desktop\source separation\Annotaion_GUI\functions')
                % Detected_subtrains={GT_1,GT_2,GT_3,GT_4};
                % Detected_refs={GT_1_ref,GT_2_ref,GT_3_ref,GT_4_ref};
                % Detected_subtrains = Detected_subtrains(~cellfun(@isempty, Detected_subtrains));
                % Detected_refs = Detected_refs(~cellfun(@isempty, Detected_refs));
                % 
                % for i=1:length(Detected_subtrains)
                % 
                %     ind=Detected_subtrains{i}; 
                %     ind_ref=Detected_refs{i};
                %     % Detections(t_index).Confidence(i)={Confidence};
                %     Detections(t_index).ToAs(i)={ind'};
                %     Detections(t_index).Pkk(i)={Power_estimates(test,ind,F_ds)};
                %     Detections(t_index).ref(i)={ind_ref};
                %     Detections(t_index).ICI(i)={diff(sort(ind))};
                %     Detections(t_index).wav_avg(i)={Waveform_average(ind,F_ds,test,0)};
                % end

                for i=1:length(Detected_subtrains)

                    ind=Detected_subtrains{i}; 
                    % Detections(t_index).Confidence(i)={Confidence};
                    Detections(t_index).ToAs(i)={t1+locs(ind)'};
                    Detections(t_index).Pkk(i)={Power_estimates(test,t1+locs(ind),F_ds)};
                    Detections(t_index).ref(i)={ref_ToAs(ind)};
                    Detections(t_index).ICI(i)={diff(sort(t1+locs(ind)))};
                    Detections(t_index).wav_avg(i)={Waveform_average(t1+locs(ind),F_ds,test,0)};
                end
         
    end
    
    if Eval_flag
        subtarin_evaluation_results=[TP TP_Baggenstos ; FP FP_Baggenstos ; MD MD_Baggenstos ; IDSw IDSw_Baggenstos];
        FP_analysis=[FP_clicks FP_clicks_Baggenstos];
    end
       
    Detections(t_index).Confidence={};
    Detections(t_index).ToAs={};
    Detections(t_index).Pkk={};
    Detections(t_index).ref={};
    Detections(t_index).ICI={};
    Detections(t_index).wav_avg={};

end



