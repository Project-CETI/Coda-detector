
function IPI_vec=MPS_extract(MPS_max,Fs,Y_zoom,Locs)

%   AUTHOR:         Guy Gubnitsky
%   DATE:           February 2024
%   DESCRIPTION:

%   This function gets a set of clicks and output their IPI estimation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    W=MPS_max*Fs;
    for j=1:length(Locs)
        Y_bank=Y_zoom(int32(Locs(j)*Fs-W):int32(Locs(j)*Fs+W));
        [MPS_vec(j),~]=MPS_estimate(Y_bank,Fs);
    end

    IPI_vec=1e3*MPS_vec;
            
end



