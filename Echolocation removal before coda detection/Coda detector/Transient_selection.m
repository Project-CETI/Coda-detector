function [Locs,Pks,Amp,SNRs]=Transient_selection(Y_filtered,ey_norm,locs,pks,F_ds,SNR_thresh,NOT)

%   AUTHOR:         Guy Gubnitsky
%   DATE:           February 2024
%   DESCRIPTION:

%   This function receives a set of detected transients and output the most
%   intense ones that exceeds a certain SNR threshold

%% Parameters settings
        locs_samples=locs*F_ds;
        SNR_window=F_ds*0.1;%25e-3;
        SNR_window_l=F_ds*70e-3;
        crop=F_ds*3e-3;
%         NOT=30;

%% Discard transients bellow a pre-defined SNR

        Locs=[]; Pks=[]; c=0;
        Amp=[]; SNRs=[];
        for i=1:length(locs)
            lc=[];
            if locs_samples(i)>SNR_window && (locs_samples(i)+SNR_window)<length(ey_norm)
%                 tmp=Y_filtered(int32(locs_samples(i)-0.1*SNR_window):int32(locs_samples(i)+SNR_window));
%                 tmp=Y_filtered(int32(locs_samples(i)-1*SNR_window):int32(locs_samples(i))-0.2*SNR_window);
                tmp=Y_filtered(int32(locs_samples(i)-SNR_window):int32(locs_samples(i))+SNR_window);

                tmp_crop=Y_filtered(int32(locs_samples(i)-crop):int32(locs_samples(i)+crop));                            
                SNR(i)=10*log(max(abs(tmp_crop))/median(abs(tmp)));
%                 SNR(i)=max(abs(tmp_crop));
                Yf = fft(tmp_crop);
                L=length(Yf);
                P2 = abs(Yf/L);
                P1 = P2(1:int32(L/2));
                P1(2:end-1) = 2*P1(2:end-1);
                f = F_ds*(0:(L/2))/L;
                fm=1e-3*f(find(P1==max(P1)));
%                 if locs_samples(i)>SNR_window_l & locs_samples(i)+SNR_window_l<length(Y_filtered)
%                     tmp_l=Y_filtered(int32(locs_samples(i)-SNR_window_l):int32(locs_samples(i)+SNR_window_l));
%                     [pc,lc] =findpeaks(tmp_l,F_ds,'MinPeakDistance',5e-3,'MinPeakHeight',max(tmp_l)*0.4);
%                 end
                if SNR(i)>SNR_thresh %& fm<18 %& length(lc)<4
                % if max(abs(tmp_crop))>20e-3 %0.025
                    c=c+1;
                    Locs(c)=locs(i);
                    Pks(c)=pks(i);
                    Amp(c)=max(abs(tmp_crop));
                    SNRs(c)=SNR(i);
                end
            end
        end

%% Pick most the NOT (30 by default) most intense transients

       if length(Pks)>NOT
           [Pks2,I] = maxk(Pks,NOT);
           Locs2=Locs(I);
           Amp2=Amp(I);
           SNRs2=SNRs(I);
           Pks=Pks2;
           Locs=Locs2;
           Amp=Amp2;
           SNRs=SNRs2;
       end
       
%% Sort detections by time of arrival
       
       [Locs,LI]=sort(Locs);
       Pks=Pks(LI);
       Amp=Amp(LI);

end