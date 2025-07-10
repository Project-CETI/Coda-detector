function Detections=run_subtrain_detect(roi,F_weights,Buffer_length,test,F_ds,All_objs,snr_lim,Tag_flag)

%   AUTHOR:         Guy Gubnitsky
%   DATE:           February 2025
%   DESCRIPTION:

% This function provides click separation within buffers. The function gets
% a long audio and divides it into short time buffers and perform source assignment of identified clicks within each buffer.
% This separation is performed by clustering the clicks using a linear assinment (LA) approach.

% The output of this function is a sturct containing the clicks' attributes
% of each source in each buffer

%   INPUT:
%   > roi                - Scalar: region of interest (roi)- defines the time window [in sec] around clicks for analyzing their surface echo
%   > F_weights          - Vector of 1X4 with the weights given to the attributes of classes 1-4 based on their relative information gain.
%   > Buffer_length      - Scalar: Analysis buffer length [sec]
%   > test               - Vector reprsenting a band-passed signal
%   > F_ds               - Scalar: smpale rate of the signal (test)
%   > All_objs           - struct of 1X5 containing the GMM parameters of the clicks' similarity attributes
%   > snr_lim            - Scalar: minimum allowed snr

%   OUTPUT:
%   > Detections         - Sturct containing the clicks' attributes of each source in each buffer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    t_index=0;
    t1=-Buffer_length;
    t2=0;   

    while(1)               
        t1=t1+Buffer_length;
        t2=t2+Buffer_length;
        t_index=t_index+1;
        if t2*F_ds>length(test)
            break;
        end
        buffer=test(t1*F_ds+1:t2*F_ds);

        [pks,locs] =findpeaks(buffer,F_ds,'MinPeakDistance',5.5e-3,'MinPeakHeight',snr_lim);
                        
        El_inds=1.1*roi;
        pks(locs<El_inds | locs>(Buffer_length-El_inds))=[];
        locs(locs<El_inds | locs>(Buffer_length-El_inds))=[];
   
        pulse_len=1e-3;
        SNR2=[]; fm=[];
        SNR_min=2.55; %2.3;
        for i=1:length(locs)
            click=buffer(int32((locs(i)-pulse_len)*F_ds):int32((locs(i)+pulse_len)*F_ds));
            noise1=[buffer(int32((locs(i)-2*pulse_len)*F_ds):int32((locs(i)-pulse_len)*F_ds)) ; buffer(int32((locs(i)+pulse_len)*F_ds):int32((locs(i)+2*pulse_len)*F_ds))];
            noise2=buffer(int32((locs(i)-3*pulse_len)*F_ds):int32((locs(i)-pulse_len)*F_ds));
            noise3=buffer(int32((locs(i)+pulse_len)*F_ds):int32((locs(i)+3*pulse_len)*F_ds));
            [~,min_noise]=min([median(abs(noise1)) median(abs(noise2)) median(abs(noise3))]);
            Noise_ops={noise1,noise2,noise3};
            noise=Noise_ops{min_noise};
            [minlen,idx_min]=min([length(click) length(noise)]);
            if idx_min==1
                noise=noise(1:minlen);
            elseif idx_min==2
                click=click(1:minlen);
            end     
            Tz=buffer(int32(locs(i)*F_ds-0.1e-3*F_ds):int32(locs(i)*F_ds+0.1e-3*F_ds));
            fm(i)=waveform_features_extraction(Tz,F_ds);
            SNR(i)=snr(click,noise);
            SNR2(i)=log(max(abs(click))/median(abs(noise)));                     
        end
        pks(SNR2<SNR_min | fm>16e3)=[];
        locs(SNR2<SNR_min | fm>16e3)=[];


        % if length(locs)>3 
        %     [ref_ToAs,locs,pks]=remove_multipass(F_ds,buffer,locs,pks,t1,roi);
        % else
        %     [ref_ToAs,~]=determine_reflection_ToAs(buffer,locs,0,roi,F_ds,0);
        % end

        [ref_ToAs,~]=determine_reflection_ToAs(buffer,locs,0,roi,F_ds,0);
                
        mode=0;
        if length(locs)>2
            Detected_subtrains=subtrain_detect(roi,ref_ToAs,F_weights,locs,buffer,F_ds,All_objs,mode,Tag_flag); 
            Detected_subtrains=Merge_chains(Detected_subtrains,locs);
        else
            Detected_subtrains={};
        end

        
        subplot(2,2,1); hold on;
        for i=1:length(Detected_subtrains)
            ind=Detected_subtrains{i};              
            dotH=plot(t1+locs(ind),pks(ind),'*','LineWidth',2); hold on; 
               pause(0.05);  % calls DRAWNOW implicitly
               set(dotH, 'XData', t1+locs(ind), 'YData', pks(ind)); hold on;
        end
             
        subplot(2,2,2);
        for i=1:length(Detected_subtrains)
            ind=Detected_subtrains{i};              
            dotH=plot(t1+locs(ind),ref_ToAs(ind),'*','LineWidth',2); hold on; 
               pause(0.05);  % calls DRAWNOW implicitly
               set(dotH, 'XData', t1+locs(ind), 'YData', ref_ToAs(ind)); hold on;
        end
        xlabel('Time [sec]'); ylabel('Slant delay [ms]'); 
        grid on; title('Click separation within buffers');
                              
        for i=1:length(Detected_subtrains)
            ind=Detected_subtrains{i}; 
            Detections(t_index).ToAs(i)={t1+locs(ind)'};
            Detections(t_index).Pkk(i)={Power_estimates(test,t1+locs(ind),F_ds)};
            Detections(t_index).ref(i)={ref_ToAs(ind)};
            Detections(t_index).ICI(i)={diff(sort(t1+locs(ind)))};
            Detections(t_index).wav_avg(i)={Waveform_average(t1+locs(ind),F_ds,test,0)};
        end
        if ~isempty(Detected_subtrains)
            Detections(t_index).Confidence=Mean_JF(Detections(t_index),buffer,F_ds,Buffer_length,t_index);
        end
         
    end
    
    Detections(t_index).Confidence=[];
    Detections(t_index).ToAs={};
    Detections(t_index).Pkk={};
    Detections(t_index).ref={};
    Detections(t_index).ICI={};
    Detections(t_index).wav_avg={};

end



