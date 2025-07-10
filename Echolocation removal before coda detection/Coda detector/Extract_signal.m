function [Y,Yp]=Extract_signal(Y_zoom,ey_norm,locs,F_ds)

%   AUTHOR:         Guy Gubnitsky
%   DATE:           February 2024
%   DESCRIPTION:

%   This function gets a measured signal and a vector of arrival time of
%   clicks and extract the clicks' waveform (Y) and peak amplitude (Yp)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    W_seg=16e-3; 
    p=0.12;
    seg_ds=round(W_seg*F_ds); 
    locs_samples=locs*F_ds; 
    for ind=1:length(locs)

        Second_pulse=[];
        if  (locs_samples(ind)+seg_ds)<length(ey_norm) && locs_samples(ind)-round(seg_ds)>0
           % Clicks_bank= ey_norm(int32(locs_samples(ind)-p*round(seg_ds)):int32(locs_samples(ind)+seg_ds));   % pick region of analysis
           Y_bank= Y_zoom(int32(locs_samples(ind)-p*round(seg_ds)):int32(locs_samples(ind)+seg_ds)); 
        elseif  (locs_samples(ind)+seg_ds)<length(ey_norm)
           % Clicks_bank= ey_norm(1:int32(locs_samples(ind)+seg_ds));   % pick region of analysis
           Y_bank= Y_zoom(1:int32(locs_samples(ind)+seg_ds));            
        else           
           % Clicks_bank= ey_norm(int32(locs_samples(ind)-p*round(seg_ds)):int32(length(ey_norm)));   % pick region of analysis
           Y_bank= Y_zoom(int32(locs_samples(ind)-p*round(seg_ds)):int32(length(ey_norm)));   % pick region of analysis
        end 

        if ind==1
            Y(:,ind)=Y_bank/max(Y_bank);
            Y_raw(:,ind)=Y_bank;
            Yp(ind)=mean([max(Y_bank) abs(min(Y_bank))]);
        else
            if length(Y_bank)<size(Y,1)
                Y(1:length(Y_bank),ind)=Y_bank/max(Y_bank);
                Y_raw(1:length(Y_bank),ind)=Y_bank;
                Yp(ind)=mean([max(Y_bank) abs(min(Y_bank))]);                 
            else
                Y(:,ind)=Y_bank/max(Y_bank);
                Y_raw(:,ind)=Y_bank;
                Yp(ind)=mean([max(Y_bank) abs(min(Y_bank))]);                
            end
        end


    end

end
