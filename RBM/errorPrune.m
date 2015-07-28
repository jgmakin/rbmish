function indexvec = errorPrune(eL,eN,s,params)


Nmods = params.Nmods;
Ndims = params.Ndims;
Nexamples = size(eL,3);
n = size(eL,2)/Nmods;
thresh = 0.1; % .05; % .2;      
%%%%%%%%%%%%%%%%%%
% you just made this number up (10% discrepancy)
% (but see notes at [131])
%%%%%%%%%%%%%%%%%%


indexvec = [];
for iExample = 1:Nexamples
    
    ss = s(iExample,:);
    GOOD = 1;
    for i = 1:n
        
        for iMode = 1:Nmods
            
            ind = (i-1)*Nmods+iMode;
            J = ntrlJacobian(ss,iMode,params);
            eNapprox = J*eL(:,ind,iExample);
            eNtrue = eN(:,ind,iExample);
            if sum( abs((eNapprox - eNtrue)./eNtrue) > thresh*ones(Ndims,1) );
                GOOD = 0;
            end

        end
    end
    if GOOD
        indexvec = [indexvec; iExample];
    end

end


end


% eL = zeros(Ndims,n*Nmods,Nexamples);