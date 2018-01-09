# rbmish
matlab code for exponential family harmoniums, RBMs, DBNs, and relata

## How to run your own EFH
To run an EFH on your own data, you basically do three things: 
  1. Write four functions, one for each of: 
     - generating the _latent_ variables of the training data (for real-world data, this will probably be null)
     - generating the observed variables of the training data
     - generating the latent and observed variables of the testing data
     - computing some kind of scalar-valued learning measure (for the main code to plot; this can be omitted) 
  2. Add a new case to `setParams.m` specifying these functions and the hyperparameters for the EFH you'd like to train; 
  3. Change the call to `setParams` in `DBNtrain.m` (~line 23) to specify your case:
  
        ```params = setParams('datatype','my_new_case');```
        
     Then you can just call `DBNtrain` from the command line.

What exactly has to be specified in your new case in `setParams.m`?  A good case to copy is `'polyphonicmusic'`.  You don't have to set everything set in that case, but you will certainly want to set the following:
 - `numsUnits`, a cell array of integers/vectors of integers
 - `typeUnits`, a cell array of cell arrays of strings
 - `Ncases`, cases/(mini)batch
 - `Nbatches`, batches/epoch
 - `NepochsMax`,
 - `setLearningSchedules`, a function that sets the learning rates and how they change over time
 - `getLatents`: a function that generates the _data_ latent variables
   - inputs: 
     - `Nexamples`, an integer specifying the number of examples per epoch
     - `dataclass`, a string specifying either `'double'` or `'gpuArray'` (in which case the code will execute the gpu-version)
   - outputs: 
     - `X` a matrix of _data_ latent variables of size `[Nexamples, ...]`
     - `Q` a structure for holding other miscellaneous parameters required to generate the training data
     
   These are really just suggested outputs, since they will only ever be passed to `params.getData`, which the user must also write (see next entry).  For real (non-synthetic) data, you probably won't set these outputs to anything at all--just set this to a null function in `setParams.m` like this:
 
        params.getLatents = @(Nexamples,yrclass,varargin)(getLatentsNull(Nexamples,yrclass,T,varargin{:}));
    
    This will set `X` to an empty tensor whose first dimension has size `Nexamples`, allowing that parameter to be passed on to `getData`.
 - `getData`: a function generating the data you will train on.
   - inputs: the outputs of `getLatents`, sc.,
     - `X`, a tensor of latent variables (i.e. of the true generative process, not the hidden vector of the EFH),
     - `Q`, a structure with whatever else you need to generate the data.
   - outputs: 
     - `R`, a tensor of data (observed) variables, of size `[Nexamples, sum(params.numsUnits{1})]`
     - `Q`, the input structure of whatever else was needed to generate the data, possibly modified to reflect something about the data generation
   During training, `X` and `Q` will be returned as outputs from `getLatents` and passed as inputs to `getData`.  So, if you're going to train on synthetic data, `getLatents` basically corresponds to the prior, and `getData` to the emission.  If you're training on real data stored somewhere, you can just write the function you assign to `getData` so that it ignores its inputs.
 - `testEFH`: a function for evaluating your EFH
   - inputs: 
     - `R`, the (observed) testing data,
     - `X`, the _latent_ variables of the testing data, if available
     - `Q`, the catch-all data-related structure (see above)
     - `wts`, the weights for the EFH as produced by `EFH.m`, a cell array of matrices
     - `params`, the master parameters returned by `setParams.m`.  
   - outputs:
     - `e`, some scalar measuring how well or badly your model is doing.
     Don't set this to be the n-step reconstruction error, because that will be reported independently in any case.  Instead, e.g., if you're using synthetic data, this function could infer the latent variables of the EFH and compare (somehow) to the "true" latent variables `X`.  Or in sequential data (synthetic or not), one could infer the hidden vector of the EFH and then use this to predict the next observable data in the sequence (see `testEFHNextFrameError.m`, e.g.), and then report the mean error.  Or etc.
 - `getTestData`: A function for generating the data on which the EFH will be tested by `params.testEFH`.  
   - inputs:
     - `dataclass`: either `'double'` or `'gpuArray'`.  The latter will instruct this function to place data on the GPU.
   - outputs: 
     - `R`, the observed variables, 
     - `X`, the latent variables (if they're available)
     - `Q`, a catch-all structure for anything else.  
     
   But since the outputs of this function _only_ ever get used by `params.testEFH`, the user is really free to set them to whatever he likes, as long as they are compatible with the function he has assigned to `params.testEFH`.
   
   
## How to run your own recurrent EFH
...
