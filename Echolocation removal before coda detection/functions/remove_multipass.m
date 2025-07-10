
function [ref_ToAs,locs,pks]=remove_multipass(F_ds,buffer,locs,pks,t1,roi)

Thresh=4;
Locs=t1+locs;
locs=Locs'-t1;
G_mat=[];
G_ind=find(diff(locs)>roi)+1;
if ~isempty(G_ind)
    if G_ind(1)>2
        G_ind=[1  G_ind  length(locs)+1];
    end
    for i=1:length(G_ind)-1
     G_mat(i,:)=[G_ind(i) G_ind(i+1)-1];
    end
    single=find(diff(G_ind)==1);
    G_mat(single,:)=[];
    
    
     El_ind=[]; Ref=[];
     ms_val=0;
    ln=1e-3;
    for i=1:size(G_mat,1)
        Lo=locs(G_mat(i,1):G_mat(i,2));
        if length(Lo)==2
          c1=buffer(int32((Lo(1)-ln)*F_ds):int32((Lo(1)+ln)*F_ds));
          c2=buffer(int32((Lo(2)-ln)*F_ds):int32((Lo(2)+ln)*F_ds));
          score=spectral_similarity(c1,c2,F_ds,0);
          if score<Thresh
             El_ind=[El_ind G_mat(i,2)];
             Ref(G_mat(i,1))=Lo(2)-Lo(1);
          end
        else
              for n=1:length(Lo)-1
                  if n~=ms_val
                      c1=buffer(int32((Lo(n)-ln)*F_ds):int32((Lo(n)+ln)*F_ds));
                      score=[];
                      for m=n+1:length(Lo)
                          c2=buffer(int32((Lo(m)-ln)*F_ds):int32((Lo(m)+ln)*F_ds));
                          score(m)=spectral_similarity(c1,c2,F_ds,0);
                      end
                      score(score==0)=100;
                      [ms,ms_ind]=min(score);
                      if ms<Thresh
                          El_ind=[El_ind G_mat(i,1)-1+ms_ind];
                          Ref(G_mat(i,1)-1+n)=Lo(ms_ind)-Lo(n);
                          ms_val=ms_ind;
                      else
                          ms_val=0;
                      end
                  end
                  
              end
    
    
        end
    end
    
    Ref=Ref*1e3;
    [ref_ToAs,~]=determine_reflection_ToAs(buffer,locs,0,roi,F_ds,0);
    
    Ref_inds=find(Ref>0);
    for i=1:length(Ref_inds)
        if abs(ref_ToAs(Ref_inds(i))-Ref(Ref_inds(i)))>2 && Ref(Ref_inds(i))>6e3 && Ref(Ref_inds(i))<1e3*roi
            ref_ToAs(Ref_inds(i))=Ref(Ref_inds(i));
        end
    end
    
    ref_ToAs(El_ind)=[];
    locs(El_ind)=[];
    pks(El_ind)=[];
else
    [ref_ToAs,~]=determine_reflection_ToAs(buffer,locs,0,roi,F_ds,0);
end

        % figure;
        % subplot(2,1,1);
        %      plot(t_buffer,buffer); hold on; grid on;    
        %     % plot(GT_1,zeros(1,length(GT_1)),'*','LineWidth',2);
        %     % plot(GT_2,zeros(1,length(GT_2)),'*','LineWidth',2);
        %     % plot(GT_3,zeros(1,length(GT_3)),'*','LineWidth',2);
        %     % plot(GT_4,zeros(1,length(GT_4)),'*','LineWidth',2);
        % 
        %     plot(t1+locs,pks,'kx','LineWidth',2);
        %     subplot(2,1,2); hold off;
        %     plot(t1+locs,ref_ToAs,'kx','LineWidth',2); grid on; hold on;
        %     % plot(GT_1',GT_1_ref,'o','LineWidth',2);
        %     % plot(GT_2,GT_2_ref,'o','LineWidth',2);
        %     % plot(GT_3,GT_3_ref,'o','LineWidth',2);
        %     % plot(GT_4,GT_4_ref,'o','LineWidth',2);
end