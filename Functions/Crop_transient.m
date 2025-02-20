function [Transient_clip,SNR_test]=Crop_transient(Y_zoom,Fs,crop_inds)

     Transient_clip=Y_zoom(crop_inds(1)-crop_inds(2):crop_inds(1)+crop_inds(2));
     if crop_inds(1)+crop_inds(3)<length(Y_zoom) & crop_inds(1)-crop_inds(3)>0
        Transient_env=Y_zoom(crop_inds(1)-crop_inds(3):crop_inds(1)+crop_inds(3));
     elseif crop_inds(1)+crop_inds(3)>length(Y_zoom) & crop_inds(1)-crop_inds(3)>0
        Transient_env=Y_zoom(crop_inds(1)-crop_inds(3):end);
     elseif crop_inds(1)-crop_inds(3)<0
         Transient_env=Y_zoom(1:crop_inds(1)+crop_inds(3));
     end
     Noise_level = median(abs(Transient_env));
     Magnitude_level = movmean(abs(Transient_clip),10e-3*Fs);
     SNR_test=max(Magnitude_level)/Noise_level;
     
end