function P=PSF_extract(IPI,Y_bank,F_ds,Inten)

%   AUTHOR:         Guy Gubnitsky
%   DATE:           February 2024
%   DESCRIPTION:

%   This function calculates the number of identified pulses in a click by analyzing the phase slope function (PSF) of a click.

%   INPUT:
%   > IPI                -  1XM vector of IPI measurements of M clicks
%   > F_ds               -  A scalar representing the sample rate of the measured signal
%   > Y_bank            -   1XN vector containing the pressure samples of a measured click

%   OUTPUT:
%   > P                 - A scalar representing the number of identified pulses in a click

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IPI=IPI(q);Y_bank=Ya(:,q);Inten=Yp(q)
exit_flag=0;
PSF_filt_norm=[]; M=[]; PSF_filt2=[]; PSF_filt=[]; xxx=[]; PSF=[];
MP_thresh=0.46;
window=round(F_ds*(IPI)*1e-3);
xxx=Y_bank; %Y_bank/max(Y_bank);
[e_Ya,~]=energyop(xxx,0);
check=length(xxx)-window-2;
if Inten>0
    if (length(xxx)-window-2)>4*window
        Rep=4*window;
    else
        Rep=(length(xxx)-window-2);
    end
else
    Rep=3*window;
end
 if check>100
    %      for n=1:(length(xxx)-window-2)
        for n=1:Rep
%             if window+n<length(e_Ya)
                xn=e_Ya(int32(n):int32(window+n));
                [Tau_w, ~]=grpdelay(xn,F_ds); 
                Tau_filt= medfilt1(Tau_w,3);
                PSF(n)=-mean(Tau_filt); 
%             else
%                 exit_flag=1;
%                 break;
%             end
         end

    if ~exit_flag
         PSF_filt=lowpass(PSF,20,F_ds);
         M=mean(PSF_filt);
         PSF_filt2=PSF_filt+abs(M)*ones(1,length(M));
         PSF_filt_norm=PSF_filt2/max(PSF_filt2);
         
         %% 
% sampling_interval = 1/F_ds; 
% desired_time_delay = -0.5*IPI; % 100 microseconds
% % Calculate the number of samples to shift
% samples_to_shift = round(desired_time_delay / sampling_interval);
% % Shift the signal to the right
% shifted_signal = circshift(PSF_filt_norm, samples_to_shift);
% % Plot the original and shifted signals

         
         
         %%

        time_x=[0:1/F_ds:(1/F_ds)*(length(xxx)-1)];
        time_psf=[0:1/F_ds:(1/F_ds)*(length(PSF_filt_norm)-1)];
%         figure; yyaxis left; plot(1e3*time_x,xxx); ylabel('Normalized signal'); set(gca,'FontSize', 12);
%         hold on; yyaxis right; plot(1e3*time_psf- desired_time_delay,shifted_signal,'-.','Linewidth',1.5); ylabel('PSF');  set(gca,'FontSize', 12);
%         xlabel('t [ms]'); set(gca,'FontSize', 12);
       [pks_F,locs_F] =findpeaks(PSF_filt_norm,F_ds,'MinPeakHeight',MP_thresh,'MinPeakDistance',(0.7*IPI)*1e-3);  % Apply instantaneous energy detector (find peaks)
%        figure; plot(1e3*time_x,xxx); ylabel('Normalized signal'); set(gca,'FontSize', 12);
%        hold on; plot(1e3*time_psf,PSF_filt_norm,'-.','Linewidth',1.5); ylabel('PSF');  set(gca,'FontSize', 12);
%         xlabel('t [ms]'); set(gca,'FontSize', 12);

%         hold on;  plot(1e3*locs_F,pks_F,'*','Linewidth',2);
        if length(pks_F)>1
            if 1e3*(locs_F(2)-locs_F(1))>1.5*IPI
                P=0;
            else
                P=length(pks_F);
            end
        else
           P=length(pks_F);
        end
    else
        P=0;
    end
end

end