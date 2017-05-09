function [R,X,Q] = generateData(Nexamples,getLatents,getData,dataclass,varargin)
% generateData  Generate input vectors for the EFH/DBN
%
%   [R,S] = generateData(Nexamples,getStimuli,getData,'gpuArray')
%
%   [R,S,Q] = generateData(Nexamples,getStimuli,getData,'double')
%
% generateData generates training/testing data R for EFHs.  
%
%   The output dimensions are:
%
%       R: (Nexamples x Nmods*N^Ndims)
%       S: (Nexamples x Ndims x Nmods)
%       G: (Nexamples x Nmods)
%
%   Q is a structure with additional information depending on the params.
%
%
%   VARIABLE ARGUMENTS:
%        string |   var   | range  | fxn                     | default val.
%
%   'dbndepth'  | iRBM    |  Z+    | current layer of DBN    | void
%    'dbnwts'   |  wts    |  wts   | produce training data   | void
%               |         |        |  for deep layers of DBN |
% 'deadneurons' | deadInds| vector | zero out responses of   | void
%               |         | of Z+  |  these units            |


%-------------------------------------------------------------------------%
% Revised: 01/02/17
%   -stripped out everything from the case statements into separate fxns.
%   Now this fxn is blind to 'datatype' as well.  Instead, it just receives
%   as input handles to functions that generate "stimuli" and data
%   ("responses").
% Revised: 12/07/16
%   -changed the logic of the whole function: it no longer knows (or cares)
%   about params.Ncases or Nbatches.  I.e., no information about training
%   should be in this function.  Instead it takes Nexamples as an input.
% Revised: 01/26/16
%   -renamed from DATAGEN.m to generateData.m
% Revised: 01/06/15 (JGM)
%   -added 'gains' as a varargin, and then used this hook in the toroidal
%   encoding (eliminating the need to set params.swing to 0 and then back
%   to its old value).
%   -updated the help section for the variable arguments
% Revised: 07/02/14 (JGM)
%   -re-wrote references to FK2link and IK2link to make use of their
%   vectorized forms
%   -this eliminated the last parfor loop!
% Revised: 06/13/14 (JGM)
%   -rewrote encodeStimuli in terms of case statements for different models
%   -added cases for 'MCDexperiment' throughout
% Revised: 05/06/14 (JGM)
%   -changed encodeStims to use vectorized respfxn.m, rather than a parfor
%   loop---so the only parallel loop left is in getStimuliCore.m.
% Revised: 12/18/13 (JGM)
%   -added case for controlled trajectories (tracking of Lissajous curves)
% Revised: 12/17/13 (JGM)
%   -added case for hierarchical data
%   -added internal data structure Q
% Revised: 12/10/13 (JGM)
%   -radically re-worked, especially getStimuli and encodeStimuli (formerly
%   getsourcevars and encodestim), in preparation for the addition of data
%   generated according to a controlled dynamical system
%   -output variable S (formerly X) now has dimensions Ncases x Ndims x
%   Nmods x Nbatches, rather than Ncases x Ndims*Nmods x Nbatches.  This
%   matches the change in estStatsCorePP.m.
% Revised: 07/01/13 (JGM)
%   -changed getsourcevars to return a structure "datainfo" for the case of
%   a dynamic stimulus.
% Revised: 06/27/13 -BKD
%   -for dynamics, src returns full state (including vels etc.)
%   -added "flag," which marks when a restart occurs
% Revised: 02/16/12
%   -heavy revision: functionalized, rationalized, etc.
% Revised: 02/03/11
%   -added varargin for biases and gains on all modalities
%   -made the appropriate other changes to accomodate these
% Revised: 12/20/10
%   -changed to accomodate 3-modality scheme
% Revised: 12/06/10 (happy b'day)
%   -changed to accomodate variable-dimensional stimuli
% Revised: 11/29/10
%   -changed to accomodate x in *true* coordinates
% Revised: 08/23/10
%   -added varargin for funny params
% Revised: 07/05/10
%   -consolidated params into setParams
%   -eliminated randpos fxn :-|
% Revised: 06/07/10
%   -accounted for edge effects, normalized the grid, etc.
% Created: 06/04/10
%   by JGM
%-------------------------------------------------------------------------%

%%% This function can be eliminated entirely.  It's here for legacy
%%% purposes.

% get stimuli and their responses
[X,Q] = getLatents(Nexamples,dataclass);
[R,Q] = getData(X,Q);


end







