function [Detected_codas,Coda_U]=Plot_Detections(F_ds,U_max_all,U_T,All_codas,Y_bpf,t_bpf,Rec_header,DF,Plot_flag)

   Det_inds=find(U_max_all>U_T);
   U_dets=U_max_all(Det_inds);
   All_codas_dets=All_codas(Det_inds);
   [Detected_codas,Coda_U]=Coda_clusters_annotation(F_ds,U_dets,All_codas_dets,Y_bpf,t_bpf,Rec_header,DF,Plot_flag);
   Coda_U=unique(round(Coda_U,4));
   Detected_codas=fliplr(Detected_codas);
   
   xlabel('time [sec]'); ylabel('Amplitude');
   set(gca,'FontSize', 12);
   
end
