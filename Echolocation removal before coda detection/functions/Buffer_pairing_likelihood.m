function L=Buffer_pairing_likelihood(Features,Buffer_Params)

    %% input Features: [ICI, Orientation, Spatial, DTW, Spectral]
    constraint_flag=0;
    % L_prev=0;
    Co=[0.18 0.18 0.28 0.18 0.18];
    for j=1:5
        x=Features(j);
        if j==2
            objA=Buffer_Params.Orientation.object;
            N_val=Buffer_Params.Orientation.normalize;
        elseif j==1
            objA=Buffer_Params.ICI.object;
            N_val=Buffer_Params.ICI.normalize;
        elseif j==3
            objA=Buffer_Params.Spatial.object;
            N_val=Buffer_Params.Spatial.normalize; 
            constraint_flag=abs(x)>0.08; %0.25; %0.08;
        elseif j==4 
            objA=Buffer_Params.DTW.object;
            N_val=Buffer_Params.DTW.normalize;  
        elseif j==5
            objA=Buffer_Params.Spectral.object;
            N_val=Buffer_Params.Spectral.normalize;  
        end
        GMM(j)=Co(j)*pdf(objA,x)/N_val;

        % L=L_prev+GMM;
        % L_prev=L;
    end
    % L=L/5;
    L=sum(GMM);
    if constraint_flag
        L=eps;
    end
end