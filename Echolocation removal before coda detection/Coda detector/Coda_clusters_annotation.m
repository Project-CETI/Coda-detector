
% function Detected_codas=Coda_clusters_annotation(All_codas,Y_buffer,ttt,Name,RecF)
function [Detected_codas,Coda_U]=Coda_clusters_annotation(F_ds,U_max_all,All_codas,Y_buffer,ttt,Name,RecF,Plot_flag)

%          Y_buffer=Y_bpf; ttt=t_bpf; Name=Rec_header; RecF=DF;
%          Plot_flag=1; U_max_all=U_dets;  All_codas=All_codas_dets;
                  
All_codas_save=All_codas;
c_ind=0; Coda={}; L_intersect=[]; Coda_U=[];
All_codas=All_codas_save;

All_codas=fliplr(All_codas);
U_max_all=fliplr(U_max_all);

for i=1:size(All_codas,2)
    for j=1:size(All_codas,2)
       L_intersect(j)=length(intersect(All_codas{i},All_codas{j}));
       if j==i
           L_intersect(j)=0;
       end
    end
    
    Full_coda_ind=find(L_intersect==max(L_intersect));
    if L_intersect(Full_coda_ind(1))>1
        c_ind=c_ind+1;
        Coda(c_ind)={round(All_codas{Full_coda_ind(1)},3)'};
        Coda_U(c_ind)=U_max_all(Full_coda_ind(1));
    elseif ~isempty(All_codas{i}) && sum(L_intersect)==0
        c_ind=c_ind+1;
        Coda(c_ind)={round(All_codas{i},3)'};
        Coda_U(c_ind)=U_max_all(i);
    end
    
    
        All_codas([i Full_coda_ind(1)])={[]};

%     All_codas([i find(L_intersect>1)])={[]};
    
end

% Coda{:}
Detected_codas=uniquearray(Coda);
Coda_U=unique(Coda_U);
% Detected_codas{:}

Detected_codas=fliplr(Detected_codas);
% Coda_U=fliplr(Coda_U);

% 
% Signal=Y_buffer;
% % ttt=[0:1/F_ds:(1/F_ds)*(length(Signal)-1)];
% 

if Plot_flag
    figure;
    plot(ttt,Y_buffer); hold on;
end
cd(RecF)
writecell({'Coda number','ToA [sec]'},[num2str(Name) '.xls'],'WriteMode','append');
for i=1:size(Detected_codas,2)
    Det=Detected_codas{i};
    if Plot_flag
        Y_pks=[];
        for q=1:length(Det)
           Y_pks(q)=max(Y_buffer(int32(F_ds*(Det(q)-4e-3)):int32(F_ds*(Det(q)+4e-3))));
        end
%         plot(Det,zeros(1,length(Det)),'x','Linewidth',3)
        plot(Det,Y_pks,'x','Linewidth',3)

    end
    writecell({i,Det},[num2str(Name) '.xls'],'WriteMode','append');
end


% if Plot_flag
%     savefig([Name '.fig'])   
% end
end




