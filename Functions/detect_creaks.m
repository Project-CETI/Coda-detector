function Creaks=detect_creaks(Locs,Pks)
%   AUTHOR:         Guy Gubnitsky
%   DATE:           February 2024
%   DESCRIPTION:

%   This function gets a vector of arrival times and avector of peak amplitudes 
%   of M transients and outputs a vector containing the indices of suspected sperm whale creaks 
%   Note: this is an optional function based on a hueristic approach which is not described
%   in the associated paper. This approach searches a sequence of high
%   rhyrhm clicks of roughly similar peak amplitudes.

%   INPUT:
%   > Locs                  - Vector of 1XM with time of arrival in seconds for M identified transients
%   > Pks                   - Vector of 1XM containing the peaks of M identified transients

%   OUTPUT:
%   > Creaks            - Vector of 1XW with indices of suspected sperm whale coda clicks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Creaks=[];  

    ICI_Max=0.06; % Set max ICI for a sequence of creaks
    v = 1:length(Locs);
    rank=10; % Set a max size for creaks candidate clusters (Note that ranks higher than 10 may cause complexity issues)

    while(1)
        rank=rank-1;
        AllCombs = nchoosek(v,rank);
        ICI_g_max=max(diff(Locs(AllCombs),1,2),[],2);
        G_ind=find(ICI_g_max<ICI_Max);
        if ~isempty(G_ind) & rank>8 
            break
        end
        if rank<9
            break
        end
    end
        
    if ~isempty(G_ind)
        [Yamp,R_Inds]=rmoutliers(Pks);
        if mean(Pks(R_Inds))-mean(Yamp)>0.25 
            Creaks=find(R_Inds==0);
        else
            Creaks=[1:length(Locs)];
        end
    end
         
end
