function Ind_min=val_return(x,val)

    Dif=abs(x-val);
    if min(Dif)<0.03
       Ind_min=find(Dif==min(Dif));
    else
        Ind_min=[];
    end
end
