%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%% Detetion and source separation algorithm for sperm whale communication signals ("Codas")  %%%
                                    %% Main script%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; clc;

%% Folder settns 

CF=pwd; % Current folder
PF=[CF '\Functions'] ; % Program folder
cd(PF)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameters settings:

Default_flag=1; % use 1 for setting default parameters and 0 for user-defined. 
[files,DF,F_low,F_high,FsAnalyze,Tag_flag,Plot_Detections_flag,Plot_buffer,U_T,SNR_thresh,MPS_window,T_l,T_sec,rho_corr,rho_IPI,rho_I]=DB_prompt(PF,CF,Default_flag);
creaks_flag=0; % Optional: set 1 to ignore potential creaks | set 0 for normal operation

%% Choose constraints

[fr_flag,P_flag,P_fr_flag,None_flag]=Apply_constraints;


%% Variable initialization 

All_codas={[]}; % An array of cells where each element is a vector contating the arrival time [in seconds] of the clicks of a detected coda
ToA=[];         % A vector containing the arrivel time of all detected coda clicks within the audio
U_max_all=[];   % A vector indicating the likelihood score [ranges between 0 and 1] of each detected group of coda clicks

 %% Run detector

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% run detector over each audio file in the specified folder%%%%%%
for file_ind=1:length(files)  

        cd(DF)
        filename=files(file_ind).name;
        if sum(filename(end-3:end)=='flac')==4
            Rec_header=filename(1:end-5);
        elseif sum(filename(end-2:end)=='wav')==3           
            Rec_header=filename(1:end-4);
        elseif sum(filename(end-2:end)=='WAV')==3           
            Rec_header=filename(1:end-4);    
        end
        [y,Fs] = audioread(filename);                 % load recordings
        Y=y(:,1);                                     % Choose one chanel from the WRU
        cd(PF)
        if Fs > FsAnalyze
            S_factor=floor(Fs/FsAnalyze);             % Define factor for resampling towards 48khz  
        else
            S_factor=1;
        end
        File_duration=(1/Fs)*(length(Y)-1);           % Calculate duration of the loaded recording
        Y_decimated = decimate(Y,S_factor);           % Resample recording 
        F_ds=Fs/S_factor;                             % Sample frequency of the decimated recording (48khz by default)       
        Y_bpf=bandpass(Y_decimated,[F_low, F_high],F_ds); % Apply bandpass filter within the frequency range specified by the parameters F_low and F_high
        t_bpf=0:1/F_ds:(1/F_ds)*(length(Y_bpf)-1);
        In=0;       
        Iterations=floor((size(Y_decimated,1)/F_ds)/T_sec); % calculate the number of iterations for runing an analysis sliding window of size T_l in steps of T_sec across the audio

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
            if Tag_flag
               [pks,locs] =findpeaks(ey_norm,F_ds,'MinPeakDistance',5*MPS_max); % For DTag data: Detect transients with an ROI frame of 5*MPS_max
            else
               [pks,locs] =findpeaks(ey_norm,F_ds,'MinPeakDistance',2*MPS_max); % For a recorder deployed from vessel: Detect transients an ROI frame of 2*MPS_max
            end

            %%%%%%%%%%% Remove clicks close to the buffer edges %%%%%%%%%% 
            El_inds=MPS_max+2e-3;
            pks(locs<El_inds | locs>(T_l-El_inds))=[];
            locs(locs<El_inds | locs>(T_l-El_inds))=[];
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            [Locs,Pks,Amp,SNRs]=Transient_selection(Y_filtered,ey_norm,locs,pks,F_ds,SNR_thresh,20); % Eliminate transients bellow a pre-defined SNR threshold
     
            %% Optional (for reducing runtime): ignore click sequences who are likely to be associated with creaks (frequent clicks used in foraging)
            if creaks_flag
                if length(Locs)>2
                   Creaks=detect_creaks(Locs,Pks);                        
                   Locs(Creaks)=[]; Pks(Creaks)=[]; Amp(Creaks)=[]; SNRs(Creaks)=[];
                   if isempty(Creaks)
                       [Locs,Pks,Amp,SNRs]=Transient_selection(Y_filtered,ey_norm,locs,pks,F_ds,SNR_thresh,20); % Eliminate transients bellow a pre-defined SNR threshold
                   end
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
        
       [Detected_codas,Coda_U]=Plot_Detections(F_ds,U_max_all,U_T,All_codas,Y_bpf,t_bpf,Rec_header,DF,Plot_Detections_flag); % Plot coda detections

end 






