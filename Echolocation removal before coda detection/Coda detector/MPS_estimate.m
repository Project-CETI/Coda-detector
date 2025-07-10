function [MPS,MPS_second]=MPS_estimate(Y_bank,F_ds)


    width=0.45e-3;
    width_samples=round(F_ds*width,1);
    [ey,~]=energyop(Y_bank,0);             % Apply TKEO (ey is the output signal with enhanced SNR)
    ey_norm=ey/max(ey);                     % Normalize the enhaced signal ey          
    time=[0:1/F_ds:(1/F_ds)*(length(ey_norm)-1)];   
    [Pk,Lo] =findpeaks(ey_norm,F_ds,'MinPeakDistance',1.8e-3);  % Apply instantaneous energy detector (find peaks)          

    Pk(Lo<2*width)=[]; Lo(Lo<2*width)=[];
    Pk(Lo>time(end)-2*width)=[]; Lo(Lo>time(end)-2*width)=[];

%        figure;
%        plot(time,ey_norm); hold on; 
%        xlabel('time [sec]'); ylabel('TKEO'); ylim([0 1]);
%        hold on; plot(Lo,Pk,'x');
       
%        [~, I_sorted]=sort(Pk);
%        MPS=abs(Lo(I_sorted(end))-Lo(I_sorted(end-1)));
       
       
        [pk1,I1]=sort(Pk);
        Main_peak=Lo(I1(end));
        Lz=length(Pk);
                tk1=Lo(I1);
                if Lz>2
                    RP=(pk1(end-1)-pk1(end-2))/pk1(end-1);                                      
                    if RP<0.4 & abs(Lo(I1(end-2))-Main_peak)>abs(Lo(I1(end-1))-Main_peak) 
                        Second_pulse=find(Pk==pk1(end-2));
                    else
                        Second_pulse=find(Pk==pk1(end-1));
                    end
                else                   
                   Second_pulse=find(Pk==pk1(end-1));
                end
                W=abs(Main_peak-Lo(Second_pulse));
                if W>1.8e-3
                    MPS=W;
                    MPS_second=Lo(Second_pulse);
                else
                    MPS=1.8e-3;        %set min window size
                    MPS_second=Lo(I1(end));
                end

       
end