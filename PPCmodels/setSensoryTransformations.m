function [Finv,Jfxn,InfoNeutralize] = setSensoryTransformations(...
    inmod,outmod,roboparams)
% setSensoryTransformations
%
% USAGES:
%
%   [Finv,Jfxn,InfoNeutralize] =...
%       setSensoryTransformations(inmod,outmod,roboparams);
%
%   F  = setSensoryTransformations(outmod,inmod,roboparams);
%
% NB that, conceptually, the first output is Finv, to which Jfxn and
% InfoNeutralize are companions.  But one can request F instead, simply by
% inverting the order of the input arguments inmod and outmod.  Thus, this
% function can be used both in computing optimal posterior distributions
% for multisensory integration (with or w/o coordinate transformations),
% and for simply retrieving the transformation from one modality to
% another.
%
% For how to use the outputs in the former case, see the supplemental 
% material in [Makin, Fellows, Sabes; PLoS C.B., 2013].  Briefly, suppose 
%   
%   sNonntrl = F(sNeutral),
%  
% and that errors are small compared with curvature of F.  Then
%
%   F(sNeutral) - shatNonntrl \approx J*(sNeutral - Finv(shatNonntrl)).
%
% This function returns Finv and a functional version of J, Jfxn.
%
% Now, mathematically, J is most simply treated by absorbing it into the 
% information matrix via Inew = J'*I*J, which is what congruenceTransform
% does.  For coordinate transformations, something more is needed: the
% linear change of coordinates (sum) decreases one information matrix as a
% harmonic sum with the information of the additive variable (usually gaze
% angle).  That's what harmonicSum does.  The order of application of these
% two functions depends on the "input" and "output" modalities (as does the
% transform Finv and the Jacobian function Jfxn); the composition is set in
% the function InfoNeutralize, whose handle is also returned.
%
% See cumulantNeutralize.m for examples of this usage.
%
% NB: the functions Finv and Jfxn assume a single input X consisting of the
% inmod and the *gaze angle*, concatenated along the third dimension:
%
%   X:  (Nexamples x Nstates x 2)

%-------------------------------------------------------------------------%
% Created: 01/06/17
%   by JGM
%-------------------------------------------------------------------------%

%%% TO DO:
% (1) Accommodate MCDexperiment?
% 
% (2) Accommodate LTI-PPC?  The idea is to allow (in theory) for
% multisensory integration even in this case.  The outputs of this function
% would then become inputs to a Kalman filter.  See KF4PPC.m....
%
% 


switch inmod
    case 'Hand-Position'
        switch outmod
            case {'Joint-Angle','Joint-Angle-Left','Joint-Angle-Right'}
                Finv = @(S,E)(FK2link(S,roboparams,0)+E);
                Jfxn = @(S,E)(IKjacobianFast(S-E,roboparams));
                InfoNeutralize = @(T1,T2,JJ)(harmonicSum(...
                    congruenceTransform(T1,JJ),T2));
            case 'Hand-Position'
                Finv = @(S,E)(S);
                Jfxn = @(S,E)(permute(eye(size(S,2)),[3,1,2]));
                InfoNeutralize = @(T1,T2,JJ)(congruenceTransform(T1,JJ));
            otherwise
                fprintf('No transformation for pairing %s = F(%s)',inmod,outmod);
                fprintf(', just returning identity transforms\n');
                Finv = @(S,E)(S);
                Jfxn = @(S,E)(permute(eye(size(S,2)),[3,1,2]));
                InfoNeutralize = @(T1,T2,JJ)(T1);
        end
        
    case {'Joint-Angle','Joint-Angle-Left','Joint-Angle-Right'}
        switch outmod
            case {'Joint-Angle','Joint-Angle-Left','Joint-Angle-Right'}
                Finv = @(S,E)(S);
                Jfxn = @(S,E)(permute(eye(size(S,2)),[3,1,2]));
                InfoNeutralize = @(T1,T2,JJ)(congruenceTransform(T1,JJ));
            case 'Hand-Position'
                Finv = @(S,E)(IK2link(S-E,roboparams,0));
                Jfxn = @(S,E)(FKjacobianFast(S,roboparams));
                InfoNeutralize = @(T1,T2,JJ)(congruenceTransform(...
                    harmonicSum(T1,T2),JJ));
            otherwise
                fprintf('No transformation for pairing %s = F(%s)',inmod,outmod);
                fprintf(', just returning identity transforms\n');
                Finv = @(S,E)(S);
                Jfxn = @(S,E)(permute(eye(size(S,2)),[3,1,2]));
                InfoNeutralize = @(T1,T2,JJ)(T1);
        end
        
    case 'Gaze-Angle'
        %%% Gaze angle is special
        Finv = @(S,E)(E);
        Jfxn = @(S,E)(eye(size(E,2)));
        InfoNeutralize = @(T1,T2,JJ)(congruenceTransform(T2,JJ));
        
        
        
        %%%%% these cases need work
        %%%% what about E??  Decide which it should add to (not that it
        %%%% matters, as long as you're consistent)....
    case 'Motion-Dots'
        switch outmod
            case 'ICMSpolar'
                Finv = @(S,E)([atan2(S(:,2),S(:,1)), sqrt(S(:,1).^2 + S(:,2).^2)]);

                %%%
                % Jfxn = @(S,E)(IKjacobianFast(S-E,roboparams));
                % InfoNeutralize = @(T1,T2,JJ)(harmonicSum(...
                %    congruenceTransform(T1,JJ),T2));
                %%%
            case 'ICMS'
                Finv = @(S,E)(S);
            otherwise
                error('unexpected pairing: %s = F(%s) -- jgm',inmod,outmod);
        end
        
    case 'ICMSpolar'
        switch outmod
            case 'Motion-Dots'
                Finv = @(S,E)([S(:,2).*cos(S(:,1)), S(:,2).*sin(S(:,1))] + E);
            otherwise
                error('unexpected pairing: %s = F(%s) -- jgm',inmod,outmod);
        end
    case 'ICMS'
        switch outmod
            case 'Motion-Dots'
                Finv = @(S,E)(S);
            otherwise
                error('unexpected pairing: %s = F(%s) -- jgm',inmod,outmod);
        end
        
                
    otherwise
        fprintf('\nNo transformation for mod %s--just ',inmod);
        fprintf('returning identity transform\n\n');
        Finv = @(S,E)(S);
end

end
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function M3 = harmonicSum(M1,M2)
% M3 = inv(inv(M1) + inv(M2)), for *tensors* of matrices, (Nexamples x
% Ndims x Ndims), and with allowances for special cases.

if isempty(M1)
    M3 = M2;
elseif isempty(M2)
    M3 = M1;
else
    if (length(M1) == numel(M1))&&(length(M2) == numel(M2))
        M3 = 1./(1./M1 + 1./M2);
    end
    
end

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function BtrAB = congruenceTransform(TA,B)
% Carry out the operation B'*A*B for a collection of matrices A in the
% tensor TA of size (Nexamples x Ndims x Ndims).  The output has size
% (Nexamples x Mdims x Mdims), with Mdims the number of columns in B.

AB = tensorOp(permute(TA,[2,3,1]),permute(B,[2,3,1]));
BtrAB = tensorOp(permute(B,[3,2,1]),AB);
BtrAB = permute(BtrAB,[3,1,2]);

end
%-------------------------------------------------------------------------%