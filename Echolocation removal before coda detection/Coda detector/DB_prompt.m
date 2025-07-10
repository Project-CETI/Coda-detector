function [F_low,F_high,FsAnalyze,Plot_Detections_flag,Plot_buffer,U_T,SNR_thresh,MPS_window,T_l,T_sec,rho_corr,rho_IPI,rho_I]=DB_prompt(Default_flag)

    
%% Fixed parameters
    F_low = 2e3;      % BPF: lower pass-band frequency bound
    F_high = 24e3;      % BPF: higher pass-band frequency bound
    FsAnalyze = 48e3;   % Min frequency for decimation purposes
    Plot_buffer=0;          % 1- show plots for the detections of the individual buffers
    MPS_window=5;           % analysis window (in ms) for IPI estimation.  
    rho_corr=0.4;           % normalizing weight for correelation similarity factor 
    rho_IPI=0.3;            % normalizing weight for IPI similarity factor 
    rho_I=0.3;              % normalizing weight for intensity similarity factor 
 

%% Default parameters    
    if Default_flag
        Plot_Detections_flag=1; % 1- show plot for the detections of the entire record 
        U_T=0.1;                % Detection threshold
        SNR_thresh=26;        % max amplitude allowed 
        T_l=3;                  % Set duration for the analysis buffer [in seconds]
        T_sec=1;                % Set overlap of succesive frames [in seconds]

    else
%% Tune parameters
        prompt = {'Set a threshold value for the minimum permissible SNR [in dB]: '};
        dlgtitle = 'Input'; 
        dims = [1 35];
        definput = {'26'};
        answer_scr = inputdlg(prompt,dlgtitle,dims,definput);
        if isempty(answer_scr{1})
            SNR_thresh=26;
        else
            SNR_thresh=str2num(answer_scr{1});
        end


        prompt = {'Set a time buffer length [sec]: '};
        dlgtitle = 'Input'; 
        dims = [1 35];
        definput = {'3'};
        answer_scr = inputdlg(prompt,dlgtitle,dims,definput);
        if isempty(answer_scr{1})
            T_l=3;
        else
            T_l=str2num(answer_scr{1});
        end

        prompt = {'Set overlap of succesive frames [in seconds]:'};
        dlgtitle = 'Input'; 
        dims = [1 35];
        definput = {'1'};
        answer_scr = inputdlg(prompt,dlgtitle,dims,definput);
        if isempty(answer_scr{1})
            T_sec=1;
        else
            T_sec=str2num(answer_scr{1});
        end


        prompt = {'Set a detection threshold U_T: '};
        dlgtitle = 'Input'; 
        dims = [1 35];
        definput = {'0.7'};
        answer_scr = inputdlg(prompt,dlgtitle,dims,definput);

        if isempty(answer_scr{1})
            U_T=0.1;
        else
           U_T=str2num(answer_scr{1});
        end

        answer_plot = questdlg('Visualize results?:','Plot activation','Yes','No','Other');
        % Handle response
        switch answer_plot
            case 'Yes'
                Plot_Detections_flag = 1;
            case 'No'
                Plot_Detections_flag = 0;
        end
    end

end