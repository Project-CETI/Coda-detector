function [trajectories,id_j_ToAs]=run_train(Detections,Buffer_Params,test,F_ds,roi,Tag_flag)

%   AUTHOR:         Guy Gubnitsky
%   DATE:           February 2025
%   DESCRIPTION:

% This function combines assigned clicks from consecutive time buffers into complete click trains. 
% This is performed by tracking the click subsets via a MAP approach
% The function outputs structs containing the formed click trains and their
% clicks' arrival times.

%   INPUT:
%   > Detections         - Sturct containing the clicks' attributes of each source in each buffer
%   > roi                - Scalar: region of interest (roi)- defines the time window [in sec] around clicks for analyzing their surface echo
%   > Buffer_Params      - struct containing the GMM parameters of the clicks' similarity attributes in the buffer level
%   > test               - Vector reprsenting a band-passed signal
%   > F_ds               - Scalar: smpale rate of the signal (test)

%   OUTPUT:
%   > trajectories    - Sturct containing the indices of the formed click trains
%   > id_j_ToAs       - Sturct containing the clicks' arrival times within the formed trains

% These two outputs are structured such that
% sol=cell2mat(id_j_ToAs(trajectories{j})'); is a vector containing the
% arrival times of all clicks in train with index j
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    trajectories={};
    id_j_ToAs={};
    l_prev=0;
    transition_arcs=[];
    detection_arcs=[];
    co_1=0;
    co_2=0;
    for i=1:length(Detections)
    
        if i==1
            C_i_en=0;
        else
            Pen=length(Detections(i).ICI)/(length(Detections(i).ICI)+length(Detections(i-1).ICI));
            % Pen=sum(Detections(i).Confidence)/(sum([Detections(i-1).Confidence Detections(i).Confidence]));
            if Pen==0
                Pen=eps;
            end
            C_i_en=-log(Pen);
        end
        if i<length(Detections)
            Pex=length(Detections(i).ICI)/(length(Detections(i).ICI)+length(Detections(i+1).ICI));
            % Pex=sum(Detections(i).Confidence)/(sum([Detections(i).Confidence Detections(i+1).Confidence]));
            if Pex==0
                Pex=eps;
            end        
           C_i_ex=-log(Pex);
        else
           C_i_ex=0;
        end
    
    
        for j=1:length(Detections(i).ICI)
            id_j=l_prev+j;       
            co_1=co_1+1;
            id_j_ToAs(co_1)=Detections(i).ToAs(j);
            if length(id_j_ToAs{co_1})>2
                 C_i=log(0.01/0.99);
                 % C_i=-log(sum(Detections(i).Confidence));
            else
                C_i=log(0.99/0.01);
            end
            detection_arcs(co_1,:)=[id_j C_i_en C_i_ex C_i];
            Pkk_current=median(Detections(i).Pkk{j});
            ref_current=median(Detections(i).ref{j});
            ICI_current=median(Detections(i).ICI{j});
            wav_avg_current=Detections(i).wav_avg{j};
    
            for k=1:length(Detections(i+1).ICI)
                id_k=l_prev+length(Detections(i).ICI)+k;
                Pkk_next=median(Detections(i+1).Pkk{k});
                ref_next=median(Detections(i+1).ref{k});
                ICI_next=median(Detections(i+1).ICI{k});
                wav_avg_next=Detections(i+1).wav_avg{k};

                ToA_prev=sort(Detections(i).ToAs{j});
                ToA_follow=sort(Detections(i+1).ToAs{k});
                ICI_constraint=ToA_follow(1)-ToA_prev(end)>2;
                alpha_Pkk=log(Pkk_next/Pkk_current);
                if ref_next<2 || ref_current<2
                    alpha_ref=1;
                else
                     alpha_ref=log(ref_next/ref_current);
                end
                alpha_ICI=log(ICI_next/ICI_current);
                
                [wa1_aligned,wa2_aligned] = alignsignals(wav_avg_current,wav_avg_next);
                alpha_DTW=dtw(wa1_aligned/max(wa1_aligned),wa2_aligned/max(wa2_aligned));
                alpha_Spectral=spectral_similarity(wa1_aligned,wa2_aligned,F_ds,0);
    
                Features=[alpha_ICI alpha_Pkk alpha_ref alpha_DTW alpha_Spectral];
                L(j,k)=-log(Buffer_pairing_likelihood(Features,Buffer_Params));
                if ICI_constraint
                    L(j,k)=-log(eps);
                end
    
                co_2=co_2+1;
                transition_arcs(co_2,:)=[id_j id_k L(j,k)];
                          
            end        
        end
        l_prev=l_prev+length(Detections(i).ICI);
    
    end
    
    if isempty(transition_arcs)
        c=0;
        for x=1:length(Detections)
            tmp=Detections(x).ToAs;
            if ~isempty(tmp)
            c=c+1;
            trajectories(c)={c};
            end
        end
        trajectories = trajectories';
    else    

        writematrix(detection_arcs,'detection_arcs.txt')
        writematrix(transition_arcs,'transition_arcs.txt')
    
    
        %% load affinity scores
        d_m = dlmread('detection_arcs.txt');
        t_m = dlmread('transition_arcs.txt');
        
        
        %% call mcc4mot

        [trajectories, ~] = mcc4mot(d_m,t_m);      
        
        %% Performance evaluation:
              
        EL=[];
        for j=1:length(trajectories)
           if length(trajectories{j})<2 %3
               EL=[EL j];
           end
        end
        
        trajectories(EL)=[];
        subplot(2,2,4);  
        xlabel('Time [sec]'); ylabel('Slant delay [ms]'); 
        title('Click train formation');
        for j=1:length(trajectories)
            if Tag_flag
                sol=cell2mat(id_j_ToAs(trajectories{j})); 
            else
                sol=cell2mat(id_j_ToAs(trajectories{j})'); 
            end
            [ref_ToAs,~]=determine_reflection_ToAs(test,sol,0,roi,F_ds,0);
            Slant=ref_ToAs;   
            subplot(2,2,4);
            dotH=plot(sol,Slant,'*','LineWidth',2);  hold on;
            pause(0.01);
            set(dotH, 'XData', sol, 'YData', Slant); hold on;
            subplot(2,2,1); hold on;
            Pks=Peaks_extract(test,sol,F_ds);   
            dotH3=plot(sol',Pks,'*','LineWidth',2);  hold on;
            pause(0.01);
            set(dotH3, 'XData', sol, 'YData', Pks); hold on;
        end
        
       
        init=[];
        for j=1:length(trajectories)
            if Tag_flag
                sol=cell2mat(id_j_ToAs(trajectories{j}));
            else
                sol=cell2mat(id_j_ToAs(trajectories{j})');
            end
            init(j)=sol(1);
        end
        [~,I]=sort(init);
        trajectories=trajectories(I);
    end
    subplot(2,2,4);
    xlabel('Time [sec]');
    ylabel('Slant delay [ms]');
    title('Click train formation');

    subplot(2,2,1);
    xlabel('Time [sec]');
    ylabel('Voltage [V]');
end