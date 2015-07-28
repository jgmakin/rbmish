function InfoTU = setInfoMatrix(SigmaX0,SigmaV0,m)
% InfoTU = setInfoMatrix(SigmaX0,SigmaV0,m)

%-------------------------------------------------------------------------%
% Cribbed: 10/22/13
%   from KF4PPC
%   by JGM
%-------------------------------------------------------------------------%


if any(isinf(SigmaX0(:)))
    InfoX0 = zeros(size(SigmaX0));
else
    InfoX0 = inv(SigmaX0);
end

if any(isinf(SigmaV0(:)))
    InfoV0 = zeros(size(SigmaV0));
else
    InfoV0 = inv(SigmaV0);
end

InfoTU = zeros(2*m);
InfoTU(1:m,1:m) = InfoX0;
InfoTU(m+1:2*m,m+1:2*m) = InfoV0;

end

