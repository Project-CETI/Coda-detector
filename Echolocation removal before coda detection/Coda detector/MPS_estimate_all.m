function MPS=MPS_estimate_all(Y_bank,F_ds)

R=4;
    width=0.45e-3;
    width_samples=round(F_ds*width,1);
    [ey,~]=energyop(Y_bank,0);             % Apply TKEO (ey is the output signal with enhanced SNR)
    ey_norm=ey/max(ey);                     % Normalize the enhaced signal ey          
    time=[0:1/F_ds:(1/F_ds)*(length(ey_norm)-1)];   
    [Pk,Lo] =findpeaks(ey_norm,F_ds,'MinPeakDistance',1.8e-3);  % Apply instantaneous energy detector (find peaks)          

    Pk(find(Lo<2*width))=[]; Lo(find(Lo<2*width))=[];
    Pk(find(Lo>time(end)-2*width))=[]; Lo(find(Lo>time(end)-2*width))=[];
    
    [pk1,I1]=sort(Pk,'descend');
    if length(pk1)>R
       Pk=pk1(1:R);
       Lo=Lo(I1(1:R));
    else
       Pk=pk1;
       Lo=Lo(I1);
    end

%        figure;
%        plot(time,ey_norm); hold on; 
%        xlabel('time [sec]'); ylabel('TKEO'); ylim([0 1]);
%        hold on; plot(Lo,Pk,'x');
       
       
       for i=2:length(Lo)
           MPS(i-1)=abs(Lo(1)-Lo(i));
       end
        
        
        
end