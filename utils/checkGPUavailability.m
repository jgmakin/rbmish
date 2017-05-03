function GPUAVAILABLE = checkGPUavailability

%-------------------------------------------------------------------------%
% Cribbed: 01/02/17
%   from Matt J. at matlabcentral
%-------------------------------------------------------------------------%

try    
    gpuArray(1);
    GPUAVAILABLE = true; 
catch
    GPUAVAILABLE = false;
end

end