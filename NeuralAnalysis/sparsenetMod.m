% sparsenet.m - simulates the sparse coding algorithm
%
% Before running you must first define A and load X.
% See the README file for further instructions.

%-------------------------------------------------------------------------%
% Revised: 03/13/11
%   by JGM
% Created: ??/??/??
%   by B.A. Olshausen
%-------------------------------------------------------------------------%


% params
Ntrials = 10000;
batch_size = 1;



eta = 0.3; % 1.0;
noise_var = 0.01;
beta = 2.2;
sigma = 0.316;
tol = .01;

VAR_GOAL = 0.1;
S_var = VAR_GOAL*ones(Ncoeff,1);
var_eta = .001;
alpha = .02;
gain = sqrt(sum(A.*A))';

display_every = 100; % 10;
% h = display_network(A,S_var);


% main loop
for t = 1:Ntrials

    % calculate coefficients for these data via conjugate-gradient routine
    S = cgf_fitS(A,X,noise_var,beta,sigma,tol,0,0,0);

    % calculate residual error
    E = X - A*S;

    % update bases
    dQ = zeros(Nneurons,Ncoeff);
    for i = 1:batch_size
        dQ = dQ + U'*E(:,i)*S(:,i)';
    end
    dQ = dQ/batch_size;
    Q = Q + eta*dQ;
    A = U*Q;

    % normalize bases to match desired output variance
    for i = 1:batch_size
        S_var = (1-var_eta)*S_var + var_eta*S(:,i).*S(:,i);
    end
    gain = gain .* ((S_var/VAR_GOAL).^alpha);
    normA = sqrt(sum(A.*A));
    for i = 1:Ncoeff
        A(:,i) = gain(i)*A(:,i)/normA(i);
    end



    % display
    if (mod(t,display_every)==0)
        figure(1)
        imagesc(reshape(S,5,4));
        figure(2)
        imagesc(Q)
        drawnow;
        % display_network(A,S_var,h);
    end

end
