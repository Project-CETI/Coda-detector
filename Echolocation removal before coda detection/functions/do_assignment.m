function [C,cost] = do_assignment (D, p, u)

% p=30; u=1e3;

    Nd = size(D, 1);
    D2 = eye(Nd).*p;
    % D2 = eye(Nd)*p;
    % % D2 = eye(Nd)*log(p(:,2));
    D2 (find(D2(:)==0)) = u;
    A = [D D2];
    [rowc, cost] = lapjv (A, 0.00001);
    
    C=zeros(size(A));
    for i=1:length(rowc)
       C(i,rowc(i))=1;
    end
return