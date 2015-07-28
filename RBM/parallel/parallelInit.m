function [pool,HADBEENCLOSED] = parallelInit



HADBEENCLOSED = isempty(gcp('nocreate'));
if HADBEENCLOSED

    [~, PCname] = system('hostname');
    switch  strtrim(PCname)
        case 'kobayashi-maru'
            numWorkers = 1; %%% will this work??
        case 'CUPCAKE'
            numWorkers = 4;
        case 'Themistocles'
            numWorkers = 2;
        otherwise
            numWorkers = 4;
            fprintf('defaulting to 4 parallel workers -- jgm\n');
    end
    pool = parpool(numWorkers); 

else
    pool = []; %%% kind of hacky
end


% change the workers' default seeds!!
spmd
    RandStream.setGlobalStream ...
        (RandStream('mt19937ar','seed',sum(100*(1 + labindex/numlabs)*clock)));
end


end