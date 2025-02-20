function fr=waveform_features_extraction(s,Fs,Ind_min) 

%   AUTHOR:         Guy Gubnitsky
%   DATE:           February 2024
%   DESCRIPTION:

%   This function calculates the resonant frequency of a set of input clicks.

%   INPUT:
%   > s                  -  1XM vector representing a measured signal
%   > Fs                 -  A scalar representing the sample rate of the measured signal
%   > Ind_min            -  1XN vector of indices representing the sampled arrival time of N clicks

%   OUTPUT:
%   > fr                 - Vector of MX1 representing the pdf of all input data points.

    frame=2e-3;
    W_size=int32(Fs*frame);
    Analysis_w=s(Ind_min-W_size:Ind_min+W_size);
    t_w=(0:1/Fs: (1/Fs)*(length(Analysis_w)-1))';
    fois    = linspace(2e3,24e3,100);
    srord   = [1, 3];
    Pf=aslt(Analysis_w, Fs, fois, 3, srord, 0); 
    Pfn=Pf/max(max(Pf));
     
    [pi,~]=find(Pfn==1); 
    fr=fois(pi);
end