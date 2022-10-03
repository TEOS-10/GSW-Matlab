function X = gsw_quadprog(H,f,A,B,Aeq,Beq,lb,ub,X0,options)

% gsw_quadprog                                        Quadratic programming 
%==========================================================================
%
% DESCRIPTION:
% This function attempts to solve the quadratic programming 
%   problem:
%
%            min 0.5*x'*H*x + f'*x   subject to:  A*x <= b 
%             x      
%
% This function was adapted from an old version of Matlab's quadprog. 
% It only alows for the active set solving.
%
%==========================================================================

defaultopt = struct( ...
    'Algorithm','active-set', ...
    'Diagnostics','off', ...
    'Display','final', ...
    'HessMult',[], ... 
    'MaxIter',[], ...    
    'MaxPCGIter','max(1,floor(numberOfVariables/2))', ...   
    'PrecondBandWidth',0, ... 
    'TolCon',1e-8, ...
    'TolFun',[], ...
    'TolPCG',0.1, ...    
    'TolX',100*eps, ...
    'TypicalX','ones(numberOfVariables,1)' ...    
    );

% If just 'defaults' passed in, return the default options in X
if nargin == 1 && nargout <= 1 && isequal(H,'defaults')
   X = defaultopt;
   return
end

% Check for non-double inputs
msg = isoptimargdbl('QUADPROG', {'H','f','A','b','Aeq','beq','LB','UB','X0'}, ...
                                  H,  f,  A,  B,  Aeq,  Beq,  lb,  ub,  X0);
if ~isempty(msg)
    error('gsw_quadprog: NonDoubleInput');
end
                     
% Internal Hessian-multiply function
mtxmpy = @hessMult_optimInternal;
usrSuppliedHessMult = false;

% Set the constraints up: defaults and check size
[nineqcstr,numberOfVariablesineq] = size(A);
[neqcstr,numberOfVariableseq] = size(Aeq);
if isa(H,'double') && ~usrSuppliedHessMult
   % H must be square and have the correct size 
   nColsH = size(H,2);
   if nColsH ~= size(H,1)
      error('gsw_quadprog: NonSquareHessian');
   end
else % HessMult in effect, so H can be anything
   nColsH = 0;
end

% Check the number of variables. The check must account for any combination of these cases:
% * User provides HessMult
% * The problem is linear (H = zeros, or H = [])
% * The objective has no linear component (f = [])
% * There are no linear constraints (A,Aeq = [])
% * There are no, or partially specified, bounds 
% * There is no X0
numberOfVariables = max([length(f),nColsH,numberOfVariablesineq,numberOfVariableseq]);

if numberOfVariables == 0
    % If none of the problem quantities indicate the number of variables,
    % check X0, even though some algorithms do not use it.
    if isempty(X0)
        error('gsw_quadprog: EmptyProblem');
    else
        % With all other data empty, use the X0 input to determine
        % the number of variables.
        numberOfVariables = length(X0);
    end
end

ncstr = nineqcstr + neqcstr;

if isempty(f)
    f = zeros(numberOfVariables,1);
else 
    % Make sure that the number of rows/columns in H matches the length of
    % f under the following conditions:
    % * The Hessian is passed in explicitly (no HessMult)
    % * There is a non-empty Hessian
    if ~usrSuppliedHessMult && ~isempty(H)
        if numel(f) ~= nColsH
            error('gsw_quadprog: MismatchObjCoefSize');
        end
    end
end

if isempty(H)
    H = zeros(numberOfVariables,numberOfVariables);
end
if isempty(A)
    A = zeros(0,numberOfVariables);
end
if isempty(B)
    B = zeros(0,1);
end
if isempty(Aeq)
    Aeq = zeros(0,numberOfVariables);
end
if isempty(Beq)
    Beq = zeros(0,1);
end

% Expect vectors
f = f(:);
B = B(:);
Beq = Beq(:);

if ~isequal(length(B),nineqcstr)
    error('gsw_quadprog: InvalidSizesOfAAndB')
elseif ~isequal(length(Beq),neqcstr)
    error('gsw_quadprog: InvalidSizesOfAeqAndBeq')
elseif ~isequal(length(f),numberOfVariablesineq) && ~isempty(A)
    error('gsw_quadprog: InvalidSizesOfAAndF')
elseif ~isequal(length(f),numberOfVariableseq) && ~isempty(Aeq)
    error('gsw_quadprog: InvalidSizesOfAeqAndf')
end

[X0,lb,ub,msg] = checkbounds(X0,lb,ub,numberOfVariables);
if ~isempty(msg)
   X = X0; 
   return
end

% Check that all data is real
if ~(isreal(H) && isreal(A) && isreal(Aeq) && isreal(f) && ...
     isreal(B) && isreal(Beq) && isreal(lb) && isreal(ub) && isreal(X0))
    error('gsw_quadprog: ComplexData')
end

isLPProblem = false;
% Perform checks on H
if isa(H,'double') && ~usrSuppliedHessMult
   if norm(H,'inf') == 0 
      % an LP problem.  
      isLPProblem = true;
   else
      % Make sure Hessian matrix is symmetric
      if norm(H-H',inf) > eps
         H = (H+H')*0.5;
      end
   end
end

if ~isa(H,'double')
    error('gsw_quadprog: NoHessMult')
end

if isempty(X0)
    X0 = zeros(numberOfVariables,1);
end

% Set default value of MaxIter for qpsub
defaultopt.MaxIter = 200;

% Create options structure for qpsub
qpoptions.MaxIter = optimget(options,'MaxIter',defaultopt,'fast');
qpoptions.TolCon = [];

if issparse(H) || issparse(A) || issparse(Aeq) % Passed in sparse matrices
    H = full(H);
    A = full(A);
    Aeq = full(Aeq);
end

verbosity = 0;
caller = 'quadprog';

[X,lambdaqp,exitflag,output,~,~,msg] =  ...
    qpsub(H,f,[Aeq;A],[Beq;B],lb,ub,X0,neqcstr,verbosity,caller,ncstr,numberOfVariables,qpoptions);
    
end
