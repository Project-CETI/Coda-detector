function El_inds=Eliminate_ref(Detected_subtrains,locs)

    El=[]; El_inds=[];
    comb=nchoosek(1:length(Detected_subtrains),2);
    for q=1:size(comb,1)
        a=locs(Detected_subtrains{comb(q,1)})';
        b=locs(Detected_subtrains{comb(q,2)})';
           
         ix = abs(a(:) - b(:).');
         UN=unique(ix);
         creteria=sum(1e3*unique(ix)<30)/max([length(a) length(b)]);
        if  creteria>0.49
            [ind,~]=ismember(ix,UN(1e3*unique(ix)<30));
            [x,y]=find(ind==1);
            if sum(a(x)-b(y))<0
                El=[El b];
            else
                El=[El a];
            end
        end
    end
    El=unique(El);
    [~,El_inds]=intersect(locs,El);

end