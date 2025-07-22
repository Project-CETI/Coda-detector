%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%% Detetion and source separation algorithm for sperm whale locomotion signals ("Echolocation clicks")  %%%
                                    %% Main script%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; clc;

%% Arrange folder paths
PF=pwd;
PF=[PF '\functions']; % program folder path
Rec_folder=uigetdir([PF '\', 'Select audio Folder']); % audio folder path
cd(PF)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load model parameters
 load All_objs.mat
 load Buffer_Params
 load F_weights

 %% Default parameters

 Buffer_length=3;    % Analysis buffer length [sec] 
 lone_p=20;          % Lone penalty (used for click train formation)      
 ITI_min=-4;         % minimum allowed time gap (Inter-train interval (ITI)) between click trains [sec]
 ITI_max=30;         % maximum allowed time gap (Inter-train interval (ITI)) between click trains [sec]
 Amplitude_lim=0.5e-3; % minimum allowed signal amplitude
 c_mean=1520;  % mean sound speed [m/s]
 h_depth=20;  % hydrophone depth [m]
 roi=(2*h_depth)/c_mean;    % region of interest (roi)- defines the time window [in sec] around clicks for analyzing their surface echo
 Tag_flag=0;    % flag for adapting parameters to tag recordings (1) | or remote recordings (0)

 %% Run on selected files
 
cd(Rec_folder);   
 if Tag_flag  
     Files=dir('*.flac');
 else
      Files=dir('*.wav');
 end

%% Run the algorithm on each audio file within the selected folder
for fi=1:1%length(Files)
        filename=Files(fi).name;
        cd(Rec_folder);
        [test,F_ds] = audioread(filename); % load signal 
        cd(PF); 
        test=bandpass(test(:,1),[2e3,24e3],F_ds); % apply band pass filter                
        t_test=(0:1/F_ds: (1/F_ds)*(length(test)-1))';

        %% select threshold interactively

        figure('Units','normalized','OuterPosition',[0 0 1 1]);
        plot(t_test,test); hold on; title('Press a random buttom in the keyboard, then select a desired amplitude threshold')
        xlabel('time [sec]'); ylabel('Amplitude'); 

        % Enable zoom
        zoom on;
        pause;  % Wait for user to press a key
        zoom off;
        
        % Get point from user
        [~, Amplitude_lim] = ginput(1);
        
        close;



        %% show signal
        figure('units','normalized','outerposition',[0 0 1 1]) 
        subplot(2,2,1);
        plot(t_test,test); hold on;
        xlabel('time [sec]'); ylabel('Voltage [v]');  
        XL = get(gca, 'XLim'); xlim(XL); 
                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% run separation algorithm %%%%%%%%%%%%%%%%%%%%%%%%
         %%  Click separation within buffers
           
        Detections=run_subtrain_detect(roi,F_weights,Buffer_length,test,F_ds,All_objs,Amplitude_lim,Tag_flag);    
        
        %%  Click train formation (sequences' association between buffers)
        [trajectories,id_j_ToAs]=run_train(Detections,Buffer_Params,test,F_ds,roi,1);
       
        close;


%% Click elimination prior coda detection

% echo_det=[];
% for i=1:length(Detections)
%     for j=1:length(Detections(i).ToAs)
%         echo_det=[echo_det Detections(i).ToAs{j}];
%     end    
% end


echo_det=[];
for i=1:length(trajectories)
    if length(cell2mat(id_j_ToAs(trajectories{i})))>9
        echo_det=[echo_det cell2mat(id_j_ToAs(trajectories{i}))];
    end
end

% echo_det=[];
% for i=1:size(Det,2)
%     echo_det=[echo_det Det(i).ToAs'];
% end

PF_echolocation=PF;
cd ..
CF=pwd; % Current folder
PF=[CF '\Coda detector'] ; % Program folder
cd(PF)

Default_flag=1; % use 1 for setting default parameters and 0 for user-defined. 
[F_low,F_high,FsAnalyze,Plot_Detections_flag,Plot_buffer,U_T,SNR_thresh,MPS_window,T_l,T_sec,rho_corr,rho_IPI,rho_I]=DB_prompt(Default_flag)
creaks_flag=1; % Optional: set 1 to ignore potential creaks | set 0 for normal operation

%% Choose constraints

[fr_flag,P_flag,P_fr_flag,None_flag]=Apply_constraints;

%% Variable initialization 

All_codas={[]}; % An array of cells where each element is a vector contating the arrival time [in seconds] of the clicks of a detected coda
ToA=[];         % A vector containing the arrivel time of all detected coda clicks within the audio
U_max_all=[];   % A vector indicating the likelihood score [ranges between 0 and 1] of each detected group of coda clicks

 %% Run detector

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% run detector over each audio file in the specified folder%%%%%%

 
        if sum(filename(end-3:end)=='flac')==4
            Rec_header=filename(1:end-5);
        elseif sum(filename(end-2:end)=='wav')==3           
            Rec_header=filename(1:end-4);
        elseif sum(filename(end-2:end)=='WAV')==3           
            Rec_header=filename(1:end-4);    
        end
        Fs=F_ds;
        Y=test;                                     % Choose one chanel from the WRU
        
        File_duration=(1/Fs)*(length(Y)-1);           % Calculate duration of the loaded recording
        Y_bpf=Y;
        t_bpf=0:1/F_ds:(1/F_ds)*(length(Y_bpf)-1);
        In=0;       
        Iterations=floor((size(Y,1)/F_ds)/T_sec); % calculate the number of iterations for runing an analysis sliding window of size T_l in steps of T_sec across the audio
        % Iterations=250;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Run detector frame by frame%%%%%%%%%%%%%%%%        
        for Buffer_ind=1:Iterations-1              
           if Buffer_ind==1
               Y_filtered=Y_bpf(1:int32(T_l*F_ds));
           elseif Buffer_ind<Iterations-1      
               Y_filtered=Y_bpf(int32(In*F_ds):int32((In+T_l)*F_ds)); 
           else    
                Y_filtered=Y_bpf(int32(In*F_ds):end);
           end                  
           clc
            fprintf(['\nProcessing: ' num2str(round(100*(Buffer_ind/((Iterations-1))),2)) '%%']);             
            [ey,~]=energyop(Y_filtered,0);          % Apply TKEO (ey is the output signal with enhanced SNR)
            ey_norm=ey/max(ey);                     % Normalize the enhaced signal ey          
            ty=0:1/F_ds:(1/F_ds)*(length(Y_filtered)-1); 
            te=0:1/F_ds:(1/F_ds)*(length(ey_norm)-1); 
            MPS_max=MPS_window*1e-3;     % Set the maximum plausible IPI of sperm whale clicks
            % if Tag_flag
            %    [pks,locs] =findpeaks(ey_norm,F_ds,'MinPeakDistance',5*MPS_max); % For DTag data: Detect transients with an ROI frame of 5*MPS_max
            % else
            %    [pks,locs] =findpeaks(ey_norm,F_ds,'MinPeakDistance',2*MPS_max); % For a recorder deployed from vessel: Detect transients an ROI frame of 2*MPS_max
            % end
            if Tag_flag
               [pks,locs] =findpeaks(Y_filtered,F_ds,'MinPeakDistance',5*MPS_max,'MinPeakHeight',Amplitude_lim); % For DTag data: Detect transients with an ROI frame of 5*MPS_max
            else
               [pks,locs] =findpeaks(Y_filtered,F_ds,'MinPeakDistance',2*MPS_max,'MinPeakHeight',Amplitude_lim); % For a recorder deployed from vessel: Detect transients an ROI frame of 2*MPS_max
            end

            %%%%%%%%%%% Remove clicks close to the buffer edges %%%%%%%%%% 
            El_inds=MPS_max+2e-3;
            pks(locs<El_inds | locs>(T_l-El_inds))=[];
            locs(locs<El_inds | locs>(T_l-El_inds))=[];

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            [Locs,Pks,Amp,SNRs]=Transient_selection(Y_filtered,ey_norm,locs,pks,F_ds,SNR_thresh,20); % Eliminate transients bellow a pre-defined SNR threshold

        %  buffer=Y_filtered;
        % [pks,locs] =findpeaks(buffer,F_ds,'MinPeakDistance',5.5e-3,'MinPeakHeight',snr_lim);
        %     El_inds=MPS_max+2e-3;
        %     pks(locs<El_inds | locs>(T_l-El_inds))=[];
        %     locs(locs<El_inds | locs>(T_l-El_inds))=[];
        % 
        % pulse_len=1e-3;
        % SNR2=[]; fm=[];
        % SNR_min=3; %2.3;
        % for i=1:length(locs)
        %     click=buffer(int32((locs(i)-pulse_len)*F_ds):int32((locs(i)+pulse_len)*F_ds));
        %     noise1=[buffer(int32((locs(i)-2*pulse_len)*F_ds):int32((locs(i)-pulse_len)*F_ds)) ; buffer(int32((locs(i)+pulse_len)*F_ds):int32((locs(i)+2*pulse_len)*F_ds))];
        %     noise2=buffer(int32((locs(i)-3*pulse_len)*F_ds):int32((locs(i)-pulse_len)*F_ds));
        %     noise3=buffer(int32((locs(i)+pulse_len)*F_ds):int32((locs(i)+3*pulse_len)*F_ds));
        %     [~,min_noise]=min([median(abs(noise1)) median(abs(noise2)) median(abs(noise3))]);
        %     Noise_ops={noise1,noise2,noise3};
        %     noise=Noise_ops{min_noise};
        %     [minlen,idx_min]=min([length(click) length(noise)]);
        %     if idx_min==1
        %         noise=noise(1:minlen);
        %     elseif idx_min==2
        %         click=click(1:minlen);
        %     end     
        %     Tz=buffer(int32(locs(i)*F_ds-0.1e-3*F_ds):int32(locs(i)*F_ds+0.1e-3*F_ds));
        %     fm(i)=waveform_features_extraction2(Tz,F_ds);
        %     SNR(i)=snr(click,noise);
        %     SNR2(i)=log(max(abs(click))/median(abs(noise)));                     
        % end
        % 
        % 
        % pks(SNR2<SNR_min | fm>15e3)=[];
        % locs(SNR2<SNR_min | fm>15e3)=[];
        % Locs=locs'; Pks=pks';


            echo_det_ROI=echo_det(echo_det>In & echo_det<In+T_l);
            Tr_toas=In+Locs;
            rem_ind=[];
            for Cand_inds=1:length(Tr_toas)
                 Ind_min=val_return(echo_det_ROI,Tr_toas(Cand_inds));
                 if ~isempty(Ind_min)
                     rem_ind=[rem_ind Cand_inds];
                 end
            end

            Pks(rem_ind)=[];
            Locs(rem_ind)=[];

            NOT=30;
           if length(Pks)>NOT
               [Pks2,I] = maxk(Pks,NOT);
               Locs2=Locs(I);
               Pks=Pks2;
               Locs=Locs2;
           end
           [Locs,LI]=sort(Locs);
           Pks=Pks(LI);

            %% Optional (for reducing runtime): ignore click sequences who are likely to be associated with creaks (frequent clicks used in foraging)
            if creaks_flag
                if length(Locs)>2
                   Creaks=detect_creaks(Locs,Pks);                        
                   Locs(Creaks)=[]; Pks(Creaks)=[];
                   % Amp(Creaks)=[]; SNRs(Creaks)=[];
                   % if isempty(Creaks)
                   %     [Locs,Pks,Amp,SNRs]=Transient_selection(Y_filtered,ey_norm,locs,pks,F_ds,SNR_thresh,20); % Eliminate transients bellow a pre-defined SNR threshold
                   % end
                end
            end
            %% Evaluation of feature similarities
           if length(Locs)>2
               IPI=MPS_extract(MPS_max,F_ds,Y_filtered,Locs);  % Evaluate the IPI of each click candidate
               [Ya,Yp]=Extract_signal(Y_filtered,ey_norm,Locs,F_ds); % Extract clicks' waveform (Ya) and peak amplitude (Yp)

               cluster={};
               s1=[]; beta=[]; gamma=[]; sim=[];
               %% Calculate similarity matrix
                for i=1:size(Ya,2)-1
                    for j=i+1:size(Ya,2)
                          s1(i,j)=max(abs(crosscorr(Ya(:,i),Ya(:,j),'NumLags',150)));    % Clicks similarity in terms of waveform correlation
                          beta(i,j)=1-(abs(Yp(i)-Yp(j))/max([Yp(i) Yp(j)]));             % Similarity of amplitudes
                          gamma(i,j)=1-(abs(IPI(i)-IPI(j))/max([IPI(i) IPI(j)]));        % IPI similarity
                          sim(i,j)=rho_corr*s1(i,j)+rho_I*beta(i,j)+rho_IPI*gamma(i,j);  % Weighted sum of all the above simiarity measures.
                    end
                end

               %% Cluster and detect codas
               
              [L_hat_inds,L_hat,U_max]=Cluster_codas(sim,IPI,Locs,Pks); % Detecet and seperate codas
               
             %% Remove clusters whose clicks' features do not staisfy the constraints imposed on the resonant requency, fr, and the number of interpulses, P.
             
             %Initialize variables%
             P_score=[]; % A vector. Each element contains the median of the number of interpulses of each click in a detected coda
             fr_score=[]; % A vector. Each element contains the average of the resonant frequency of each click in a detected coda
             TOA_current=[]; % A vector of arival times of clciks in a detected coda

               if ~None_flag
                   for D_inds=find(U_max>U_T)
                   Score=[]; fr_coda=[];
                    for q=L_hat_inds{D_inds}
                          P=PSF_extract(IPI(q),Ya(:,q),F_ds,Yp(q));
                          fr=waveform_features_extraction(Y_filtered,F_ds,Locs(q)*F_ds);
                          Score=[Score P];
                          fr_coda=[fr_coda 1e-3*fr];
                    end
                    P_score(D_inds)=mean(Score);
                    fr_score(D_inds)=median(fr_coda);
                   end
                   if P_flag
                      El=find(P_score>2.7);
                   elseif fr_flag
                         El=find(fr_score<12);
                   elseif P_fr_flag
                         El=find(P_score>2.7 & fr_score<12);
                   end
               else
                     El=find(U_max>U_T);
               end
              
             if ~isempty(El)
                 for p_inds=El
                       C_plot=L_hat{p_inds};
                       TOA_current=[TOA_current C_plot(1,:)];
                       U_max_all=[U_max_all U_max(p_inds)];
                       All_codas=[{In+C_plot(1,:)'} All_codas(:)'];                       
                 end
                 ToA=[ToA In+TOA_current];
             end
           end
        %%
           In=In+T_sec; % Move the sliding window in one step of size T_sec
        end       
        
        Det=unique(ToA); 
        U_max_all=fliplr(U_max_all);
        
     
        DF=Rec_folder;

       
        cd(PF)
        [Detected_codas,Coda_U]=Plot_Detections(F_ds,U_max_all,U_T,All_codas,Y_bpf,t_bpf,Rec_header,DF,Plot_Detections_flag); % Plot coda detections
        close;

%% Final stage - remove duplicate detections (in cases where main and multi-pulses detected as two seperate codas)
          Output_name=[Rec_header '.xls'];
          T_haifa=readcell(Output_name);
          Inds=2:size(T_haifa,1);

          if Plot_Detections_flag
              figure;
              plot(t_bpf,Y_bpf); hold on;
          end
        if ~isempty(Inds)
            cd(PF)
            Inds=Windowing_el(Inds,T_haifa);
            cd(Rec_folder)
            delete(Output_name);
            writecell({'Coda number','ToA [sec]'},Output_name,'WriteMode','append');            
            for j=1:length(Inds)
                 click_init=T_haifa{Inds(j),2};
                 ToAs=T_haifa(Inds(j),2:end);
                 Coda_rythm=0;
                 click=ToAs{1};
                 for k=1:length(ToAs)-1
                     if ~ismissing(ToAs{k+1})
                         click(k+1)=ToAs{k+1};
                     end
                 end
                 cd(PF_echolocation)
                 Pk=Peaks_extract(test,click,Fs);
                 if Plot_Detections_flag
                    plot(click,Pk,'*','LineWidth',2);
                 end
                 cd(Rec_folder)
                 writecell({j,click},Output_name,'WriteMode','append');
            end
        end
         
end