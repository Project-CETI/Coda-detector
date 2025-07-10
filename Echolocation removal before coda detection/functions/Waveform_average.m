function wa =Waveform_average(locs,F_ds,test,plot_flag)
        
      if plot_flag
         figure;
      end
      y2=test(int32(F_ds*(locs(1)-2e-3)):int32(F_ds*(locs(1)+2e-3)));
      if length(locs)>2
          for i=1:length(locs)-1
               e1=y2;
               e2=test(int32(F_ds*(locs(i+1)-2e-3)):int32(F_ds*(locs(i+1)+2e-3)));
               [y1,y2] = alignsignals(e1,e2);
               if i==1
                  if plot_flag
                     plot(y1); hold on; plot(y2);
                  end
                  wave(i)={y1}; wave(i+1)={y2};
                  L(i)=length(y1); L(i+1)=length(y2);
               else
                   if plot_flag
                       plot(y2); hold on;
                   end
                   wave(i+1)={y2}; L(i+1)=length(y2);
               end           
          end
          Lmin=min(L);
          for i=1:length(locs)
              tmp=wave{i};
              waves_cut(i,:)=tmp(1:Lmin);
          end
          if plot_flag
             figure; plot(mean(waves_cut,1));
          end
          wa=mean(waves_cut,1);
      else
          wa=y2;
      end
end