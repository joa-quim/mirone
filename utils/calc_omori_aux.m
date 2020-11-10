function  varargout = calc_omori_aux(fun,varargin)
% Auxiliary functions to calc_omori.m NOT to be distributed in source form.
% Include this file in calc_omori.m when compiling it, otherwise distribute
% in the pcode form.

	[varargout{1:nargout}] = feval(fun, varargin{:});

% --------------------------------------------------------------------------------------
function [H,pValue,KSstatistic] = kstest2(x1, x2, alpha, tail)
%KSTEST2 Two-sample Kolmogorov-Smirnov goodness-of-fit hypothesis test.
%   H = KSTEST2(X1,X2,ALPHA,TAIL) performs a Kolmogorov-Smirnov (K-S) test 
%   to determine if independent random samples, X1 and X2, are drawn from 
%   the same underlying continuous population. ALPHA and TAIL are optional
%   scalar inputs: ALPHA is the desired significance level (default = 0.05); 
%   TAIL indicates the type of test (default = 0). H indicates the result of
%   the hypothesis test:
%      H = 0 => Do not reject the null hypothesis at significance level ALPHA.
%      H = 1 => Reject the null hypothesis at significance level ALPHA.
% 
%   Let F1(x) and F2(x) be the empirical distribution functions from sample 
%   vectors X1 and X2, respectively. The 2-sample K-S test hypotheses and 
%   test statistic are:
%
%   Null Hypothesis: F1(x) = F2(x) for all x
%      For TAIL =  0 (2-sided test), alternative: F1(x) not equal to F2(x).
%      For TAIL =  1 (1-sided test), alternative: F1(x) > F2(x).
%      For TAIL = -1 (1-sided test), alternative: F1(x) < F2(x).
%
%   For TAIL = 0, 1, and -1, the test statistics are T = max|F1(x) - F2(x)|,
%   T = max[F1(x) - F2(x)], and T = max[F2(x) - F1(x)], respectively.
%
%   The decision to reject the null hypothesis occurs when the significance 
%   level, ALPHA, equals or exceeds the P-value.
%
%   X1 and X2 are row or column vectors of lengths N1 and N2, respectively, 
%   and represent random samples from some underlying distribution(s). 
%   Missing observations, indicated by NaN's (Not-a-Number), are ignored.
%
%   [H,P] = KSTEST2(...) also returns the asymptotic P-value P.
%
%   [H,P,KSSTAT] = KSTEST2(...) also returns the K-S test statistic KSSTAT
%   defined above for the test type indicated by TAIL.
%
%   The asymptotic P-value becomes very accurate for large sample sizes, and
%   is believed to be reasonably accurate for sample sizes N1 and N2 such 
%   that (N1*N2)/(N1 + N2) >= 4.

% Author(s): R.A. Baker, 08/14/98
% Copyright 1993-2002 The MathWorks, Inc. 
% $Revision: 1.5 $   $ Date: 1998/01/30 13:45:34 $

%
% References:
%   (1) Massey, F.J., "The Kolmogorov-Smirnov Test for Goodness of Fit",
%         Journal of the American Statistical Association, 46 (March 1956), 68-77.
%   (2) Miller, L.H., "Table of Percentage Points of Kolmogorov Statistics",
%         Journal of the American Statistical Association, (March 1951), 111-121.
%   (3) Conover, W.J., "Practical Nonparametric Statistics", 
%         John Wiley & Sons, Inc., 1980.
%   (4) Press, W.H., et. al., "Numerical Recipes in C", 
%         Cambridge University Press, 1992.
 
if (nargin < 2)    error(' At least 2 inputs are required.');   end

% Ensure each sample is a VECTOR.
[rows1 , columns1]  =  size(x1);
[rows2 , columns2]  =  size(x2);

if (rows1 ~= 1) && (columns1 ~= 1) 
    error(' Sample ''X1'' must be a vector.');
end

if (rows2 ~= 1) && (columns2 ~= 1) 
    error(' Sample ''X2'' must be a vector.');
end

% Remove missing observations indicated by NaN's, and 
% ensure that valid observations remain.
x1  =  x1(~isnan(x1));
x2  =  x2(~isnan(x2));
x1  =  x1(:);
x2  =  x2(:);

if (isempty(x1))   error(' Sample vector ''X1'' is composed of all NaN''s.');   end
if (isempty(x2))   error(' Sample vector ''X2'' is composed of all NaN''s.');   end

% Ensure the significance level, ALPHA, is a scalar 
% between 0 and 1 and set default if necessary.
if (nargin >= 3) && ~isempty(alpha)
   if numel(alpha) > 1
      error(' Significance level ''Alpha'' must be a scalar.');
   end
   if (alpha <= 0 || alpha >= 1)
      error(' Significance level ''Alpha'' must be between 0 and 1.'); 
   end
else
   alpha  =  0.05;
end

% Ensure the type-of-test indicator, TAIL, is a scalar integer from 
% the allowable set {-1 , 0 , 1}, and set default if necessary.
if (nargin >= 4) && ~isempty(tail)
   if numel(tail) > 1
      error(' Type-of-test indicator ''Tail'' must be a scalar.');
   end
   if (tail ~= -1) && (tail ~= 0) && (tail ~= 1)
      error(' Type-of-test indicator ''Tail'' must be -1, 0, or 1.');
   end
else
   tail  =  0;
end

% Calculate F1(x) and F2(x), the empirical (i.e., sample) CDFs.
binEdges   = [-inf ; sort([x1;x2]) ; inf];

binCounts1 = histc (x1, binEdges);
binCounts2 = histc (x2, binEdges);

sumCounts1 = cumsum(binCounts1)./sum(binCounts1);
sumCounts2 = cumsum(binCounts2)./sum(binCounts2);

sampleCDF1 = sumCounts1(1:end-1);
sampleCDF2 = sumCounts2(1:end-1);

% Compute the test statistic of interest.
switch tail
   case  0      %  2-sided test: T = max|F1(x) - F2(x)|.
      deltaCDF = abs(sampleCDF1 - sampleCDF2);
   case -1      %  1-sided test: T = max[F2(x) - F1(x)].
      deltaCDF = sampleCDF2 - sampleCDF1;
   case  1      %  1-sided test: T = max[F1(x) - F2(x)].
      deltaCDF = sampleCDF1 - sampleCDF2;
end

KSstatistic = max(deltaCDF);

% Compute the asymptotic P-value approximation and accept or
% reject the null hypothesis on the basis of the P-value.
n1     = length(x1);
n2     = length(x2);
n      = n1 * n2 /(n1 + n2);
lambda = max((sqrt(n) + 0.12 + 0.11/sqrt(n)) * KSstatistic , 0);

if (tail ~= 0)        % 1-sided test.
    pValue = exp(-2 * lambda * lambda);
else                % 2-sided test (default).
    %  Use the asymptotic Q-function to approximate the 2-sided P-value.
   j = (1:101)';
   pValue  =  2 * sum((-1).^(j-1).*exp(-2*lambda*lambda*j.^2));
   if (pValue < 0)  pValue = 0; end
   if (pValue > 1)  pValue = 1; end
end

H = (alpha >= pValue);

% -------------------------------------------------------------------------------------------------
function [X,FVAL,EXITFLAG,OUTPUT,LAMBDA,GRAD,HESSIAN] = fmincon_m(FUN,X,A,B,Aeq,Beq,LB,UB,NONLCON,options,varargin)
%FMINCON Finds a constrained minimum of a function of several variables.
%   FMINCON solves problems of the form:
%       min F(X)  subject to:  A*X  <= B, Aeq*X  = Beq (linear constraints)
%        X                       C(X) <= 0, Ceq(X) = 0   (nonlinear constraints)
%                                LB <= X <= UB            
%
%   Examples
%     FUN can be specified using @:
%        X = fmincon(@humps,...)
%     In this case, F = humps(X) returns the scalar function value F of the HUMPS function
%     evaluated at X.
%

%   Copyright 1990-2003 The MathWorks, Inc. 

defaultopt = struct('Display','final','LargeScale','on', ...
   'TolX',1e-6,'TolFun',1e-6,'TolCon',1e-6,'DerivativeCheck','off',...
   'Diagnostics','off',...
   'GradObj','off','GradConstr','off',...
   'HessMult',[],...% HessMult [] by default
   'Hessian','off','HessPattern','sparse(ones(numberOfVariables))',...
   'MaxFunEvals','100*numberOfVariables',...
   'MaxSQPIter',Inf,...
   'DiffMaxChange',1e-1,'DiffMinChange',1e-8,...
   'PrecondBandWidth',0,'TypicalX','ones(numberOfVariables,1)',...
   'MaxPCGIter','max(1,floor(numberOfVariables/2))', ...
   'TolPCG',0.1,'MaxIter',400,'OutputFcn',[]);

medium = 'medium-scale';

if nargin < 10, options=[];
	if nargin < 9, NONLCON=[];
		if nargin < 8, UB = [];
			if nargin < 7, LB = [];
				if nargin < 6, Beq=[];
					if nargin < 5, Aeq =[];	end
				end
			end
		end
	end
end

XOUT=X(:);
numberOfVariables=length(XOUT);

verbosity = 0;

% Set to column vectors
B = B(:);
Beq = Beq(:);
l = LB(:);
u = UB(:);

meritFunctionType = 0;

gradflag = 0;
hessflag = 0;
gradconstflag = 0;

% Convert to inline function as needed
if ~isempty(FUN)  % will detect empty string, empty matrix, empty cell array
   [funfcn, msg] = optimfcnchk1(FUN,'fmincon',length(varargin),gradflag,hessflag);
end

confcn{1} = '';      % ESTA
[rowAeq,colAeq]=size(Aeq);
hessflag = 0;
OUTPUT.algorithm = medium;       % ESTA

lenvlb=length(l);
lenvub=length(u);

% Ensure starting point lies within bounds
i=1:lenvlb;
lindex = XOUT(i)<l(i);
if any(lindex)
    XOUT(lindex)=l(lindex)+1e-4; 
end
i=1:lenvub;
uindex = XOUT(i)>u(i);
if any(uindex)
    XOUT(uindex)=u(uindex);
end
X(:) = XOUT;

% Evaluate function
GRAD=zeros(numberOfVariables,1);
HESS = [];

try
  f = feval(funfcn{3},X,varargin{:});
catch
  errmsg = sprintf('%s\n%s\n\n%s',...
     'FMINCON cannot continue because user supplied objective function', ...
     ' failed with the following error:', lasterr);
  error(errmsg);
end

% Evaluate constraints
c=[]; ceq =[];
cGRAD = zeros(numberOfVariables,length(c));
ceqGRAD = zeros(numberOfVariables,length(ceq));

% call algorithm
[X,FVAL,lambda,EXITFLAG,OUTPUT,GRAD,HESSIAN]=...
  nlconst1(funfcn,X,l,u,full(A),B,full(Aeq),Beq,confcn,options,defaultopt, ...
  verbosity,gradflag,gradconstflag,hessflag,meritFunctionType,...
  f,GRAD,HESS,c,ceq,cGRAD,ceqGRAD,varargin{:});
LAMBDA=lambda;

% ----------------------------------------------------------------------------
function [x,FVAL,lambda_out,EXITFLAG,OUTPUT,GRADIENT,HESS]= ...
    nlconst1(funfcn,x,lb,ub,Ain,Bin,Aeq,Beq,confcn,OPTIONS,defaultopt,...
    verbosity,gradflag,gradconstflag,hessflag,meritFunctionType,...
    fval,gval,Hval,ncineqval,nceqval,gncval,gnceqval,varargin)
%NLCONST Helper function to find the constrained minimum of a function 
%   of several variables. Called by FMINCON, FGOALATTAIN FSEMINF and FMINIMAX.
%
%   Copyright 1990-2003 The MathWorks, Inc.
%   $Revision: 1.41 $  $Date: 2003/02/10 20:48:36 $

%   Calls OPTEVAL.
%
%   meritFunctionType==5 for fseminf
%                    ==1 for fminimax & fgoalattain (can use 0, but 1 is default)
%                    ==0 for fmincon

% Initialize some parameters
lambda_out=[]; OUTPUT=[];

iter = 0;
XOUT = x(:);
% numberOfVariables must be the name of this variable
numberOfVariables = length(XOUT);
SD = ones(numberOfVariables,1); 
Nlconst = 'nlconst';
bestf = Inf; 
if isempty(confcn{1})
    constflag = 0;
else
    constflag = 1;
end
stepsize = 1;
HESS=eye(numberOfVariables,numberOfVariables); % initial Hessian approximation.
done = false; 
EXITFLAG = 1;

% Get options
tolX = optimgetfast(OPTIONS,'TolX',defaultopt);
tolFun = optimgetfast(OPTIONS,'TolFun',defaultopt);
tolCon = optimgetfast(OPTIONS,'TolCon',defaultopt);
DiffMinChange = optimgetfast(OPTIONS,'DiffMinChange',defaultopt);
DiffMaxChange = optimgetfast(OPTIONS,'DiffMaxChange',defaultopt);
if DiffMinChange >= DiffMaxChange
    errmsg =sprintf('%s %0.5g %s %0.5g%s\n%s\n',...
    'DiffMinChange options parameter is',DiffMinChange,'and DiffMaxChange is',...
    DiffMaxChange,'.','DiffMinChange must be strictly less than DiffMaxChange.');
    error(errmsg)
end
DerivativeCheck = strcmp(optimgetfast(OPTIONS,'DerivativeCheck',defaultopt),'on');
typicalx = optimgetfast(OPTIONS,'TypicalX',defaultopt) ;
if ischar(typicalx)
   if isequal(lower(typicalx),'ones(numberofvariables,1)')
      typicalx = ones(numberOfVariables,1);
   else
      error('Option ''TypicalX'' must be a numeric value if not the default.')
   end
end
maxFunEvals = optimgetfast(OPTIONS,'MaxFunEvals',defaultopt);
maxIter = optimgetfast(OPTIONS,'MaxIter',defaultopt);
% In case the defaults were gathered from calling: optimset('fmincon'):
if ischar(maxFunEvals)
    if isequal(lower(maxFunEvals),'100*numberofvariables')
        maxFunEvals = 100*numberOfVariables;
    else
        error('Option ''MaxFunEvals'' must be an integer value if not the default.')
    end
end

% For SQP problem create default sqpopt to pass to qpsub. 
% (Values not needed until we fix default MaxSQPIter in R14.)
sqpopt = []; 

% Handle bounds as linear constraints
arglb = ~isinf(lb);
lenlb=length(lb); % maybe less than numberOfVariables due to old code
argub = ~isinf(ub);
lenub=length(ub);
boundmatrix = eye(max(lenub,lenlb),numberOfVariables);
if nnz(arglb) > 0     
    lbmatrix = -boundmatrix(arglb,1:numberOfVariables);% select non-Inf bounds 
    lbrhs = -lb(arglb);
else
    lbmatrix = []; lbrhs = [];
end
if nnz(argub) > 0
    ubmatrix = boundmatrix(argub,1:numberOfVariables);
    ubrhs=ub(argub);
else
    ubmatrix = []; ubrhs=[];
end 

% Update constraint matrix and right hand side vector with bound constraints.
A = [lbmatrix;ubmatrix;Ain];
B = [lbrhs;ubrhs;Bin];
if isempty(A)
    A = zeros(0,numberOfVariables); B=zeros(0,1);
end
if isempty(Aeq)
    Aeq = zeros(0,numberOfVariables); Beq=zeros(0,1);
end

% Used for semi-infinite optimization:
s = nan; NEWLAMBDA =[]; NPOINT =[];
OLDLAMBDA = [];

x(:) = XOUT;  % Set x to have user expected size
% Compute the objective function and constraints
f = fval;
nceq = nceqval; ncineq = ncineqval;  % nonlinear constraints only
nc = [nceq; ncineq];
c = [ Aeq*XOUT-Beq; nceq; A*XOUT-B; ncineq];

% Get information on the number and type of constraints.
non_eq = length(nceq);
non_ineq = length(ncineq);
[lin_eq,Aeqcol] = size(Aeq);
[lin_ineq,Acol] = size(A);  % includes upper and lower bounds
eq = non_eq + lin_eq;
ineq = non_ineq + lin_ineq;
ncstr = ineq + eq;
% Boolean inequalitiesExist = true if and only if there exist either
% finite bounds or linear inequalities or nonlinear inequalities. 
% Used only for printing indices of active inequalities at the solution
inequalitiesExist = any(arglb) || any(argub) || size(Ain,1) > 0 || non_ineq > 0;

% Compute the initial constraint violation.
ga=[abs(c( (1:eq)' )) ; c( (eq+1:ncstr)' ) ];
if (~isempty(c))    mg=max(ga);
else                mg = 0;
end

if (isempty(f))
    error('FUN must return a non-empty objective function.')
end

OLDX=XOUT;
OLDC=c;
OLDgf=zeros(numberOfVariables,1);
gf=zeros(numberOfVariables,1);
OLDAN=zeros(ncstr,numberOfVariables);
LAMBDA=zeros(ncstr,1);
lambdaNLP = zeros(ncstr,1);
numFunEvals=1;
numGradEvals=1;
% Just for not having anoying warnings of false errors from the compiler
OLDF = [];  MATX = [];  ACTIND = [];    infeasIllPosedMaxSQPIter = [];

optimError = []; % In case we have convergence in 0th iteration, this needs a value.
%---------------------------------Main Loop-----------------------------
while ~done 
   %----------------GRADIENTS----------------
   
   if ~gradconstflag || ~gradflag || DerivativeCheck
      % Finite Difference gradients (even if just checking analytical)
      POINT = NPOINT; 
      len_nc = length(nc);
      ncstr =  lin_eq + lin_ineq + len_nc;     
      FLAG = 0; % For semi-infinite

      % Compute finite difference gradients
      %
      if (DerivativeCheck || ~gradconstflag)            % ESTA
         [gf,gnc,NEWLAMBDA,OLDLAMBDA,s]=finiteDifferences(XOUT,x,funfcn,confcn,lb,ub,f,nc,...
                               [],DiffMinChange,DiffMaxChange,typicalx,[],'all',...
                               LAMBDA,NEWLAMBDA,OLDLAMBDA,POINT,FLAG,s,varargin{:});
         gnc = gnc'; % nlconst requires the transpose of the Jacobian
      else
         % no derivative check and user-supplied constraint gradients.             
         % Estimate objective gradient only.
         gf=finiteDifferences(XOUT,x,funfcn,confcn,lb,ub,f,[],...   
                         DiffMinChange,DiffMaxChange,typicalx,[],...
                         [],[],[],[],[],[],varargin{:});
      end
      
      numFunEvals = numFunEvals + numberOfVariables;
   end  
   
   % Now add in Aeq, and A
   if ~isempty(gnc)
      gc = [Aeq', gnc(:,1:non_eq), A', gnc(:,non_eq+1:non_ineq+non_eq)];
   elseif ~isempty(Aeq) || ~isempty(A)
      gc = [Aeq',A'];
   else
      gc = zeros(numberOfVariables,0);
   end
   AN=gc';
   
   % Iteration 0 is handled separately below
   if (iter > 0) % All but 0th iteration ----------------------------------------
       % Compute the first order KKT conditions.
       if strcmp(funfcn{2},'fseminf'), lambdaNLP = NEWLAMBDA; end
       normgradLag = norm(gf + AN'*lambdaNLP,inf);
       normcomp = norm(lambdaNLP(eq+1:ncstr).*c(eq+1:ncstr),inf);
       if isfinite(normgradLag) && isfinite(normcomp)
           optimError = max(normgradLag, normcomp);
       else
           optimError = inf;
       end
       feasError  = mg;
       optimScal = 1; feasScal = 1; 
              
       %-------------TEST CONVERGENCE---------------
       
       if (optimError < tolFun*optimScal && feasError < tolCon*feasScal)
           EXITFLAG = 1;
           done = true;

           if inequalitiesExist
              % Report active inequalities
              [activeLb,activeUb,activeIneqLin,activeIneqNonlin] = ...
                  activeInequalities(c,tolCon,arglb,argub,lin_eq,non_eq,size(Ain));           
           end   
       elseif (max(abs(SD)) < 2*tolX || abs(gf'*SD) < 2*tolFun ) && (mg < tolCon || infeasIllPosedMaxSQPIter) 
           % The algorithm can make no more progress.  If feasible, compute 
           % the new up-to-date Lagrange multipliers (with new gradients) 
           % and recompute the KKT error.  Then output appropriate termination
           % message.
           if (mg < tolCon)
               lambdaNLP(:,1) = 0;
               [Q,R] = qr(AN(ACTIND,:)');
               ws = warning('off');
               lambdaNLP(ACTIND) = -R\Q'*gf;
               warning(ws);
               lambdaNLP(eq+1:ncstr) = max(0,lambdaNLP(eq+1:ncstr));
               if strcmp(funfcn{2},'fseminf'), lambdaNLP = NEWLAMBDA; end
               normgradLag = norm(gf + AN'*lambdaNLP,inf);
               normcomp = norm(lambdaNLP(eq+1:ncstr).*c(eq+1:ncstr),inf);
               if isfinite(normgradLag) && isfinite(normcomp)
                   optimError = max(normgradLag, normcomp);
               else
                   optimError = inf;
               end
               optimScal = 1;
               if optimError < tolFun*optimScal
                   EXITFLAG = 1;
               elseif max(abs(SD)) < 2*tolX
                   EXITFLAG = 1;
               else 
                   EXITFLAG = 1;
               end 
               
               if inequalitiesExist
                  % Report active inequalities
                  [activeLb,activeUb,activeIneqLin,activeIneqNonlin] = ...
                      activeInequalities(c,tolCon,arglb,argub,lin_eq,non_eq,size(Ain));  
               end
           else                         % if mg < tolCon
               EXITFLAG = -1;
           end                          % if mg < tolCon
           done = true;
       else % continue
           % NEED=[LAMBDA>0] | G>0
           if numFunEvals > maxFunEvals
               XOUT = MATX;
               f = OLDF;
               EXITFLAG = 0;
               done = true;
           end
           if iter > maxIter
               XOUT = MATX;
               f = OLDF;
               EXITFLAG = 0;
               done = true;
           end
       end 
   else % ------------------------0th Iteration----------------------------------
       
   end % if iter > 0
   
   % Continue if termination criteria do not hold or it is the 0th iteration-------------------------------------------
   if ~done 
      how=''; 
      iter = iter + 1;

      %-------------SEARCH DIRECTION---------------
      % For equality constraints make gradient face in 
      % opposite direction to function gradient.
      for i=1:eq 
         schg=AN(i,:)*gf;
         if schg>0
            AN(i,:)=-AN(i,:);
            c(i)=-c(i);
         end
      end
   
      if numGradEvals>1  % Check for first call    
         if meritFunctionType~=5,   
            NEWLAMBDA=LAMBDA; 
         end
         [ma,na] = size(AN);
         GNEW=gf+AN'*NEWLAMBDA;
         GOLD=OLDgf+OLDAN'*LAMBDA;
         YL=GNEW-GOLD;
         sdiff=XOUT-OLDX;

         % Make sure Hessian is positive definite in update.
         if YL'*sdiff<stepsize^2*1e-3
            while YL'*sdiff<-1e-5
               [YMAX,YIND]=min(YL.*sdiff);
               YL(YIND)=YL(YIND)/2;
            end
            if YL'*sdiff < (eps*norm(HESS,'fro'));
               how=' Hessian modified twice';
               FACTOR=AN'*c - OLDAN'*OLDC;
               FACTOR=FACTOR.*(sdiff.*FACTOR>0).*(YL.*sdiff<=eps);
               WT=1e-2;
               if max(abs(FACTOR))==0; FACTOR=1e-5*sign(sdiff); end
               while YL'*sdiff < (eps*norm(HESS,'fro')) && WT < 1/eps
                  YL=YL+WT*FACTOR;
                  WT=WT*2;
               end
            else
               how=' Hessian modified';
            end
         end
         
         %----------Perform BFGS Update If YL'S Is Positive---------
         if YL'*sdiff>eps
             HESS=HESS +(YL*YL')/(YL'*sdiff)-((HESS*sdiff)*(sdiff'*HESS'))/(sdiff'*HESS*sdiff);
             % BFGS Update using Cholesky factorization  of Gill, Murray and Wright.
             % In practice this was less robust than above method and slower. 
             %   R=chol(HESS); 
             %   s2=R*S; y=R'\YL; 
             %   W=eye(numberOfVariables,numberOfVariables)-(s2'*s2)\(s2*s2') + (y'*s2)\(y*y');
             %   HESS=R'*W*R;
         else
            how=' Hessian not updated';
         end
      else % First call
         OLDLAMBDA=repmat(eps+gf'*gf,ncstr,1)./(sum(AN'.*AN')'+eps);
         ACTIND = 1:eq;     
      end % if numGradEvals>1
      numGradEvals=numGradEvals+1;
   
      LOLD=LAMBDA;
      OLDAN=AN;
      OLDgf=gf;
      OLDC=c;
      OLDF=f;
      OLDX=XOUT;
      XN=zeros(numberOfVariables,1);
      GT =c;
      HESS = (HESS + HESS')*0.5;
   
      [SD,lambda,exitflagqp,outputqp,howqp,ACTIND] ...
         = qpsub(HESS,gf,AN,-GT,[],[],XN,eq,-1,Nlconst,size(AN,1),numberOfVariables,OPTIONS,sqpopt,ACTIND);
    
      lambdaNLP(:,1) = 0;
      lambdaNLP(ACTIND) = lambda(ACTIND);
      lambda((1:eq)') = abs(lambda( (1:eq)' ));
      ga=[abs(c( (1:eq)' )) ; c( (eq+1:ncstr)' ) ];
      if (~isempty(c))  mg = max(ga);
      else              mg = 0;
      end

      if strncmp(howqp,'ok',2); 
          howqp =''; 
      end
      if ~isempty(how) && ~isempty(howqp) 
          how = [how,'; '];
      end

      LAMBDA=lambda((1:ncstr)');
      OLDLAMBDA=max([LAMBDA';0.5*(LAMBDA+OLDLAMBDA)'])' ;

      %---------------LINESEARCH--------------------
      MATX=XOUT;
      MATL = f+sum(OLDLAMBDA.*(ga>0).*ga) + 1e-30;

      infeasIllPosedMaxSQPIter = strcmp(howqp,'infeasible') || strcmp(howqp,'ill posed') || strcmp(howqp,'MaxSQPIter');
     % This merit function looks for improvement in either the constraint
     % or the objective function unless the sub-problem is infeasible in which
     % case only a reduction in the maximum constraint is tolerated.
     % This less "stringent" merit function has produced faster convergence in
     % a large number of problems.
     if (mg > 0),		MATL2 = mg;
     elseif (f >= 0),	MATL2 = -1/(f+1);
	 else				MATL2 = 0;
     end
     if ~infeasIllPosedMaxSQPIter && f < 0
        MATL2 = MATL2 + f - 1;
     end
      if mg < eps && f < bestf
         bestf = f;
         bestx = XOUT;
         bestHess = HESS;
         bestgrad = gf;
         bestOptimError = optimError;
      end
      MERIT = MATL + 1;
      MERIT2 = MATL2 + 1; 
      stepsize=2;
      while ((MERIT2 > MATL2) && (MERIT > MATL) && numFunEvals < maxFunEvals)
         stepsize=stepsize/2;
         if stepsize < 1e-4,  
            stepsize = -stepsize;          
         end
         XOUT = MATX + stepsize*SD;
         x(:)=XOUT; 
      
         f = feval(funfcn{3},x,varargin{:});
         if constflag
            [nctmp,nceqtmp] = feval(confcn{3},x,varargin{:});
            nctmp = nctmp(:); nceqtmp = nceqtmp(:);
         else
            nctmp = []; nceqtmp=[];
         end
            
         nc = [nceqtmp(:); nctmp(:)];
         c = [Aeq*XOUT-Beq; nceqtmp(:); A*XOUT-B; nctmp(:)];  

         numFunEvals = numFunEvals + 1;
         ga=[abs(c( (1:eq)' )) ; c( (eq+1:length(c))' )];
         if (~isempty(c)),
            mg = max(ga);
         else
            mg = 0;
         end

         MERIT = f+sum(OLDLAMBDA.*(ga>0).*ga);
         if (mg > 0),		MERIT2 = mg;
         elseif (f >= 0),	MERIT2 = -1/(f+1);
		 else				MERIT2 = 0;
         end
         if (~infeasIllPosedMaxSQPIter && f < 0)
            MERIT2 = MERIT2 + f - 1;
         end
                                                                                                                                                                                                                            end  % line search loop
      %------------Finished Line Search-------------
   
      mf=abs(stepsize);
      LAMBDA=mf*LAMBDA+(1-mf)*LOLD;

      x(:) = XOUT;
      numFunEvals=numFunEvals+1;
   
      % Evaluate constraint gradients
      switch confcn{1}
      case 'fun'
         gnceq=[]; gncineq=[];
      case ''
         nctmp=[]; nceqtmp =[];
         gncineq = zeros(numberOfVariables,length(nctmp));
         gnceq = zeros(numberOfVariables,length(nceqtmp));
      otherwise
         error('Undefined calltype in FMINCON');
      end
      gnc_user = [gnceq, gncineq];
      gc = [Aeq', gnceq, A', gncineq];
   
   end % if ~done   
end % while ~done

% Gradient is in the variable gf
GRADIENT = gf;

% If a better solution was found earlier, use it:
if (f > bestf )
   XOUT = bestx;
   f = bestf;
   HESS = bestHess;
   GRADIENT = bestgrad;
   optimError = bestOptimError;
end

FVAL = f;
x(:) = XOUT;

OUTPUT.iterations = iter;
OUTPUT.funcCount = numFunEvals;
OUTPUT.stepsize = stepsize;
OUTPUT.algorithm = 'medium-scale: SQP, Quasi-Newton, line-search';
OUTPUT.firstorderopt = optimError;
OUTPUT.cgiterations = [];

[lin_ineq,Acol] = size(Ain);  % excludes upper and lower

lambda_out.lower=zeros(lenlb,1);
lambda_out.upper=zeros(lenub,1);

lambda_out.eqlin = lambdaNLP(1:lin_eq);
ii = lin_eq ;
lambda_out.eqnonlin = lambdaNLP(ii+1: ii+ non_eq);
ii = ii+non_eq;
lambda_out.lower(arglb) = lambdaNLP(ii+1 :ii+nnz(arglb));
ii = ii + nnz(arglb) ;
lambda_out.upper(argub) = lambdaNLP(ii+1 :ii+nnz(argub));
ii = ii + nnz(argub);
lambda_out.ineqlin = lambdaNLP(ii+1: ii + lin_ineq);
ii = ii + lin_ineq ;
lambda_out.ineqnonlin = lambdaNLP(ii+1 : end);

%--------------------------------------------------------------------------
function [activeLb,activeUb,activeIneqLin,activeIneqNonlin] = ...
    activeInequalities(c,tol,arglb,argub,linEq,nonlinEq,linIneq)
% ACTIVEINEQUALITIES returns the indices of the active inequalities
% and bounds.
% INPUT:
% c                 vector of constraints and bounds (see nlconst main code)
% tol               tolerance to determine when an inequality is active
% arglb, argub      boolean vectors indicating finite bounds (see nlconst
%                   main code)
% linEq             number of linear equalities
% nonlinEq          number of nonlinear equalities
% linIneq           number of linear inequalities
%
% OUTPUT
% activeLB          indices of active lower bounds
% activeUb          indices of active upper bounds  
% activeIneqLin     indices of active linear inequalities
% activeIneqNonlin  indices of active nonlinear inequalities
%

% We check wether a constraint is active or not using '< tol'
% instead of '<= tol' to be onsistent with nlconst main code, 
% where feasibility is checked using '<'.
finiteLb = nnz(arglb);      % number of finite lower bounds
finiteUb = nnz(argub);      % number of finite upper bounds

indexFiniteLb = find(arglb); % indices of variables with LB
indexFiniteUb = find(argub); % indices of variables with UB

% lower bounds
i = linEq + nonlinEq; % skip equalities

% Boolean vector that indicates which among the finite bounds is active
activeFiniteLb = abs(c(i + 1 : i + finiteLb)) < tol;

% indices of the finite bounds that are active
activeLb = indexFiniteLb(activeFiniteLb);

% upper bounds
i = i + finiteLb;

% Boolean vector that indicates which among the finite bounds is active
activeFiniteUb = abs(c(i + 1 : i + finiteUb)) < tol;

% indices of the finite bounds that are active
activeUb = indexFiniteUb(activeFiniteUb);

% linear inequalities
i = i + finiteUb;
activeIneqLin = find(abs(c(i + 1 : i + linIneq)) < tol); 
% nonlinear inequalities
i = i + linIneq;
activeIneqNonlin = find(abs(c(i + 1 : end)) < tol);   

% --------------------------------------------------------------------------------
function [gradf,cJac,NEWLAMBDA,OLDLAMBDA,s] = finiteDifferences(xCurrent,...
                  xOriginalShape,funfcn,confcn,lb,ub,fCurrent,cCurrent,...
                  YDATA,DiffMinChange,DiffMaxChange,typicalx,finDiffType,...
                  variables,LAMBDA,NEWLAMBDA,OLDLAMBDA,POINT,FLAG,s,varargin)
%
%  [gradf,cJac,NEWLAMBDA,OLDLAMBDA,s] = FINITEDIFFERENCES(xCurrent,...
%                  xOriginalShape,funfcn,confcn,lb,ub,fCurrent,cCurrent,...
%                  YDATA,DiffMinChange,DiffMaxChange,typicalx,finDiffType,...
%                  variables,LAMBDA,NEWLAMBDA,OLDLAMBDA,POINT,FLAG,s,varargin)
%
% computes the finite-difference gradients of the objective and
% constraint functions.
%
%  gradf = FINITEDIFFERENCES(xCurrent,xOriginalShape,funfcn,[],lb,ub,fCurrent,...
%                  [],YDATA,DiffMinChange,DiffMaxChange,typicalx,finDiffType,...
%                  variables,[],[],[],[],[],[],varargin)
%
% computes the finite-difference gradients of the objective function.
%
%
% INPUT:
% xCurrent              Point where gradient is desired
% xOriginalShape        Shape of the vector of variables supplied by the user
%                       (The value of xOriginalShape is NOT used)
% funfcn, confcn        Cell arrays containing info about objective and 
%                       constraints
% lb, ub                Lower and upper bounds
% fCurrent, cCurrent    Values at xCurrent of the function and the constraints 
%                       to be differentiated. Note that fCurrent can be a scalar 
%                       or a vector. 
% DiffMinChange, 
% DiffMaxChange         Minimum and maximum values of perturbation of xCurrent 
% finDiffType           Type of finite difference desired (only forward 
%                       differences implemented so far)
% variables             Variables w.r.t which we want to differentiate. Possible 
%
% LAMBDA,NEWLAMBDA,
% OLDLAMBDA,POINT,
% FLAG,s                Parameters for semi-infinite constraints
% varargin              Problem-dependent parameters passed to the objective and 
%                       constraint functions
%
% OUTPUT:
% gradf                 If fCurrent is a scalar, gradf is the finite-difference 
%                       gradient of the objective; if fCurrent is a vector,
%                       gradf is the finite-difference Jacobian  
% cJac                  Finite-difference Jacobian of the constraints
% NEWLAMBDA,
% OLDLAMBDA,s           Parameters for semi-infinite constraints
%
%   Copyright 1990-2003 The MathWorks, Inc.

[fNumberOfRows,fNumberOfColumns] = size(fCurrent);

% Make sure that fCurrent is either a scalar or a column vector
% (to ensure that the given fCurrent and the computed fplus 
% will have the same shape)
if fNumberOfColumns ~= 1
   fCurrent(:) = fCurrent;
else
   functionIsScalar = (fNumberOfRows == 1);
end

numberOfVariables = length(xCurrent); 

% nonEmptyLowerBounds = true if lb is not empty, false if it's empty;
% analogoulsy for nonEmptyUpperBound
nonEmptyLowerBounds = ~isempty(lb);
nonEmptyUpperBounds = ~isempty(ub);

% Make sure xCurrent and typicalx are column vectors so that the 
% operation max(abs(xCurrent),abs(typicalx)) won't error
xCurrent = xCurrent(:); typicalx = typicalx(:);
% Value of stepsize suggested in Trust Region Methods, Conn-Gould-Toint, section 8.4.3
CHG = sqrt(eps)*sign(xCurrent).*max(abs(xCurrent),abs(typicalx));
%
% Make sure step size lies within DiffminChange and DiffMaxChange
%
CHG = sign(CHG+eps).*min(max(abs(CHG),DiffMinChange),DiffMaxChange);
len_cCurrent = length(cCurrent);
cJac = zeros(len_cCurrent,numberOfVariables);  % For semi-infinite

if nargout < 3
   NEWLAMBDA=[]; OLDLAMBDA=[]; s=[];
   if (nargout == 1),	cJac=[];   end
end

% allVariables = true/false if finite-differencing wrt to one/all variables
allVariables = false;
if ischar(variables)
   if strcmp(variables,'all')
      variables = 1:numberOfVariables;
      allVariables = true;
   else
      error('Unknown value of input ''variables''.')
   end
end

% Preallocate gradf for speed 
if functionIsScalar
   gradf = zeros(numberOfVariables,1);
elseif allVariables % vector-function and gradf estimates full Jacobian 
   gradf = zeros(fNumberOfRows,numberOfVariables); 
else % vector-function and gradf estimates one column of Jacobian
   gradf = zeros(fNumberOfRows,1);
end
   
for gcnt=variables
   temp = xCurrent(gcnt);
   xCurrent(gcnt)= temp + CHG(gcnt);
         
   if (nonEmptyLowerBounds && isfinite(lb(gcnt))) || (nonEmptyUpperBounds && isfinite(ub(gcnt)))
      % Enforce bounds while finite-differencing.
      % Need lb(gcnt) ~= ub(gcnt), and lb(gcnt) <= temp <= ub(gcnt) to enforce bounds.
      % (If the last qpsub problem was 'infeasible', the bounds could be currently violated.)
      if (lb(gcnt) ~= ub(gcnt)) && (temp >= lb(gcnt)) && (temp <= ub(gcnt)) 
          if  ((xCurrent(gcnt) > ub(gcnt)) || (xCurrent(gcnt) < lb(gcnt))) % outside bound ?
              CHG(gcnt) = -CHG(gcnt);
              xCurrent(gcnt)= temp + CHG(gcnt);
              if (xCurrent(gcnt) > ub(gcnt)) || (xCurrent(gcnt) < lb(gcnt)) % outside other bound ?
                  [newchg,indsign] = max([temp-lb(gcnt), ub(gcnt)-temp]);  % largest distance to bound
                  if newchg >= DiffMinChange
                      CHG(gcnt) = ((-1)^indsign)*newchg;  % make sure sign is correct
                      xCurrent(gcnt)= temp + CHG(gcnt);
                  else
                      errmsg = sprintf('%s %d %s\n%s\n%s %0.5g%s\n',...
                          'Distance between lower and upper bounds, in dimension',gcnt,', is too small to compute', ...
                          'finite-difference approximation of derivative. Increase distance between these', ...
                          'bounds to be at least',2*DiffMinChange,'.');
                      error(errmsg)
                  end          
              end
          end
      end
   end % of 'if isfinite(lb(gcnt)) || isfinite(ub(gcnt))'
   
   xOriginalShape(:) = xCurrent;
   fplus = feval(funfcn{3},xOriginalShape,varargin{:});
   % YDATA: Only used by lsqcurvefit, which has no nonlinear constraints
   % (the only type of constraints we do finite differences on: bounds 
   % and linear constraints do not require finite differences) and thus 
   % no needed after evaluation of constraints
   if ~isempty(YDATA)    
      fplus = fplus - YDATA;
   end
   % Make sure it's in column form
   fplus = fplus(:);

   if functionIsScalar
      gradf(gcnt,1) =  (fplus-fCurrent)/CHG(gcnt);
   elseif allVariables % vector-function and gradf estimates full Jacobian 
      gradf(:,gcnt) = (fplus-fCurrent)/CHG(gcnt);
   else % vector-function and gradf estimates only one column of Jacobian
      gradf = (fplus-fCurrent)/CHG(gcnt);
   end

   if ~isempty(cJac) % Constraint gradient required
      % This is necessary in case confcn is empty, then in the comparison 
      % below confcn{2} would be out of range
      if ~isempty(confcn{1}) 
         [ctmp,ceqtmp] = feval(confcn{3},xOriginalShape,varargin{:});
         cplus = [ceqtmp(:); ctmp(:)];
      end
      % Next line used for problems with varying number of constraints
      if ~isempty(cplus)
         cJac(:,gcnt) = (cplus - cCurrent)/CHG(gcnt); 
      end           
   end
    xCurrent(gcnt) = temp;
end % for 

% --------------------------------------------------------------------------------
function [X,lambda,exitflag,output,how,ACTIND] = ...
    qpsub(H,f,A,B,lb,ub,X,neqcstr,verbosity,caller,ncstr,numberOfVariables,options,defaultopt,ACTIND)
%QP Quadratic programming subproblem. Handles qp and constrained
%   linear least-squares as well as subproblems generated from NLCONST.
%
%   X=QP(H,f,A,b) solves the quadratic programming problem:
%
%            min 0.5*x'Hx + f'x   subject to:  Ax <= b 
%             x    
%

%   Copyright 1990-2003 The MathWorks, Inc. 

% Define constant strings
NewtonStep = 'Newton';
NegCurv = 'negative curvature chol';   
ZeroStep = 'zero step';
SteepDescent = 'steepest descent';
Conls = 'lsqlin';
Lp = 'linprog';
Qp = 'quadprog';
Qpsub = 'qpsub';
how = 'ok'; 

exitflag = 1;
output = [];
iterations = 0;
if nargin < 15, ACTIND = [];
    if nargin < 13, options = []; 
    end
end

lb=lb(:); ub = ub(:);

if isempty(neqcstr), neqcstr = 0; end

LLS = 0;
if strcmp(caller, Conls)
    LLS = 1;
    [rowH,colH]=size(H);
    numberOfVariables = colH;
end
if strcmp(caller, Qpsub)
    normalize = -1;
else
    normalize = 1;
end

simplex_iter = 0;
if  norm(H,'inf')==0 || isempty(H), is_qp=0; else, is_qp=1; end

if LLS==1
    is_qp=0;
end

normf = 1;
if (normalize > 0)			% Check for lp
    if (~is_qp && ~LLS)
        normf = norm(f);
        if (normf > 0)  f = f./normf;   end
    end
end

% Handle bounds as linear constraints
arglb = ~eq(lb,-inf);
lenlb=length(lb); % maybe less than numberOfVariables due to old code
if nnz(arglb) > 0     
    lbmatrix = -eye(lenlb,numberOfVariables);
    A=[A; lbmatrix(arglb,1:numberOfVariables)]; % select non-Inf bounds
    B=[B;-lb(arglb)];
end

argub = ~eq(ub,inf);
lenub=length(ub);
if (nnz(argub) > 0)
    ubmatrix = eye(lenub,numberOfVariables);
    A=[A; ubmatrix(argub,1:numberOfVariables)];
    B=[B; ub(argub)];
end 

% Bounds are treated as constraints: Reset ncstr accordingly
ncstr=ncstr + nnz(arglb) + nnz(argub);

% Figure out max iteration count
% For linprog/quadprog/lsqlin/qpsub problems, use 'MaxIter' for this.
% For nlconst (fmincon, etc) problems, use 'MaxSQPIter' for this.

if (isempty(options) || isempty(options.MaxSQPIter))
  maxiter = 10*max(numberOfVariables,ncstr-neqcstr);
else 
  maxiter = options.MaxSQPIter;
end 

% Used for determining threshold for whether a direction will violate a constraint.
normA = ones(ncstr,1);
if (normalize > 0)
    for i=1:ncstr
        n = norm(A(i,:));
        if (n ~= 0)
            A(i,:) = A(i,:)/n;
            B(i) = B(i)/n;
            normA(i,1) = n;
        end
    end
else 
    normA = ones(ncstr,1);
end
errnorm = 0.01*sqrt(eps); 

tolDep = 100*numberOfVariables*eps;      
lambda = zeros(ncstr,1);
eqix = 1:neqcstr;

% Modifications for warm-start.
% Do some error checking on the incoming working set indices.
ACTCNT = length(ACTIND);
if isempty(ACTIND)
    ACTIND = eqix;
elseif neqcstr > 0
    i = max(find(ACTIND<=neqcstr));
    if isempty(i) || i > neqcstr % safeguard which should not occur
        ACTIND = eqix;
    elseif i < neqcstr
        % A redundant equality constraint was removed on the last
        % SQP iteration.  We will go ahead and reinsert it here.
        numremoved = neqcstr - i;
        ACTIND(neqcstr+1:ACTCNT+numremoved) = ACTIND(i+1:ACTCNT);
        ACTIND(1:neqcstr) = eqix;
    end
end
aix = zeros(ncstr,1);
aix(ACTIND) = 1;
ACTCNT = length(ACTIND);
ACTSET = A(ACTIND,:);

% Check that the constraints in the initial working set are not
% dependent and find an initial point which satisfies the initial working set.
indepInd = 1:ncstr;
if ACTCNT > 0 && normalize ~= -1
    % call constraint solver
    [Q,R,A,B,X,Z,how,ACTSET,ACTIND,ACTCNT,aix,eqix,neqcstr,ncstr, ...
            remove,exitflag]= ...
        eqnsolv(A,B,eqix,neqcstr,ncstr,numberOfVariables,LLS,H,X,f, ...
        normf,normA,verbosity,aix,ACTSET,ACTIND,ACTCNT,how,exitflag); 
    
    if ~isempty(remove)
        indepInd(remove)=[];
        normA = normA(indepInd);
    end
    
    if strcmp(how,'infeasible')
        % Equalities are inconsistent, so X and lambda have no valid values
        % Return original X and zeros for lambda.
        ACTIND = indepInd(ACTIND);
        output.iterations = iterations;
        exitflag = -1;
        return
    end
    
    if neqcstr >= numberOfVariables
        err = max(abs(A(eqix,:)*X-B(eqix)));
        if (err > 1e-8)  % Equalities not met
            how='infeasible';
            exitflag = -1;
            % Equalities are inconsistent, X and lambda have no valid values
            % Return original X and zeros for lambda.
            ACTIND = indepInd(ACTIND);
            output.iterations = iterations;
            return
        else % Check inequalities
            if (max(A*X-B) > 1e-8)
                how = 'infeasible';
                % was exitflag = 8; 
                exitflag = -1;
            end
        end
        if is_qp
            actlambda = -R\(Q'*(H*X+f));
        elseif LLS
            actlambda = -R\(Q'*(H'*(H*X-f)));
        else
            actlambda = -R\(Q'*f);
        end
        lambda(indepInd(ACTIND)) = normf * (actlambda ./normA(ACTIND));
        ACTIND = indepInd(ACTIND);
        output.iterations = iterations;
        return
    end
    
    % Check whether in Phase 1 of feasibility point finding. 
    if (verbosity == -2)
        cstr = A*X-B; 
        mc=max(cstr(neqcstr+1:ncstr));
        if (mc > 0)
            X(numberOfVariables) = mc + 1;
        end
    end
else 
    if ACTCNT == 0 % initial working set is empty 
        Q = eye(numberOfVariables,numberOfVariables);
        R = [];        Z = 1;
    else           % in Phase I and working set not empty
        [Q,R] = qr(ACTSET');
        Z = Q(:,ACTCNT+1:numberOfVariables);
    end   
end

% Find Initial Feasible Solution 
cstr = A*X-B;
if ncstr > neqcstr
    mc = max(cstr(neqcstr+1:ncstr));
else
    mc = 0;
end
if mc > eps
    quiet = -2;
    optionsPhase1 = []; % Use default options in phase 1
    ACTIND2 = 1:neqcstr;
    A2=[[A;zeros(1,numberOfVariables)],[zeros(neqcstr,1);-ones(ncstr+1-neqcstr,1)]];
    [XS,lambdaS,exitflagS,outputS,howS,ACTIND2] = ...
        qpsub([],[zeros(numberOfVariables,1);1],A2,[B;1e-5], ...
        [],[],[X;mc+1],neqcstr,quiet,Qpsub,size(A2,1),numberOfVariables+1, ...
        optionsPhase1,defaultopt,ACTIND2);
    slack = XS(numberOfVariables+1);
    X=XS(1:numberOfVariables);
    cstr=A*X-B;
    if slack > eps 
        if slack > 1e-8 
            how='infeasible';
            exitflag = -1;
        else
            how = 'overly constrained';
            exitflag = -1;
        end
        lambda(indepInd) = normf * (lambdaS((1:ncstr)')./normA);
        ACTIND = 1:neqcstr;
        ACTIND = indepInd(ACTIND);
        output.iterations = iterations;
        return
    else
        % Initialize active set info based on solution of Phase I.
        %      ACTIND = ACTIND2(find(ACTIND2<=ncstr));
        ACTIND = 1:neqcstr;
        ACTSET = A(ACTIND,:);
        ACTCNT = length(ACTIND);
        aix = zeros(ncstr,1);
        aix(ACTIND) = 1;
        if ACTCNT == 0
            Q = zeros(numberOfVariables,numberOfVariables);
            R = [];
            Z = 1;
        else
            [Q,R] = qr(ACTSET');
            Z = Q(:,ACTCNT+1:numberOfVariables);
        end
    end
end

if ACTCNT >= numberOfVariables - 1  
    simplex_iter = 1; 
end
[m,n]=size(ACTSET);

if (is_qp)
    gf=H*X+f;
    [SD, dirType] = compdir(Z,H,gf,numberOfVariables,f);
elseif (LLS)
    HXf=H*X-f;
    gf=H'*(HXf);
    HZ= H*Z;
    [mm,nn]=size(HZ);
    if mm >= nn
        [QHZ, RHZ] =  qr(HZ,0);
        Pd = QHZ'*HXf;
        % Now need to check which is dependent
        if min(size(RHZ))==1 % Make sure RHZ isn't a vector
            depInd = find( abs(RHZ(1,1)) < tolDep);
        else
            depInd = find( abs(diag(RHZ)) < tolDep );
        end  
    end
    if mm >= nn && isempty(depInd) % Newton step
        SD = - Z*(RHZ(1:nn, 1:nn) \ Pd(1:nn,:));
        dirType = NewtonStep;
    else % steepest descent direction
        SD = -Z*(Z'*gf);
        dirType = SteepDescent;
    end
else % lp
    gf = f;
    SD=-Z*Z'*gf;
    dirType = SteepDescent; 
    if norm(SD) < 1e-10 && neqcstr
        % This happens when equality constraint is perpendicular to objective function f.x.
        actlambda = -R\(Q'*(gf));
        lambda(indepInd(ACTIND)) = normf * (actlambda ./ normA(ACTIND));
        ACTIND = indepInd(ACTIND);
        output.iterations = iterations;
        return;
    end
end

%--------------Main Routine-------------------

while iterations < maxiter
    iterations = iterations + 1;
    if isinf(verbosity)
      curr_out = sprintf('Iter: %5.0f, Active: %5.0f, step: %s, proc: %s',iterations,ACTCNT,dirType,how);
        disp(curr_out); 
    end
    
    % Find distance we can move in search direction SD before a 
    % constraint is violated.
    % Gradient with respect to search direction.
    GSD=A*SD;
    
    % Note: we consider only constraints whose gradients are greater
    % than some threshold. If we considered all gradients greater than 
    % zero then it might be possible to add a constraint which would lead to
    % a singular (rank deficient) working set. The gradient (GSD) of such
    % a constraint in the direction of search would be very close to zero.
    indf = find((GSD > errnorm * norm(SD))  &  ~aix);
    
    if isempty(indf) % No constraints to hit
        STEPMIN=1e16;
        ind=[];
    else % Find distance to the nearest constraint
        dist = abs(cstr(indf)./GSD(indf));
        [STEPMIN,ind2] =  min(dist);
        ind2 = find(dist == STEPMIN);
        % Bland's rule for anti-cycling: if there is more than one 
        % blocking constraint then add the one with the smallest index.
        ind=indf(min(ind2));
    end
    
    %----------------Update X---------------------
    
    % Assume we do not delete a constraint
    delete_constr = 0;   
    
    if ~isempty(indf) && isfinite(STEPMIN) % Hit a constraint
        if strcmp(dirType, NewtonStep)
            % Newton step and hit a constraint: LLS or is_qp
            if STEPMIN > 1  % Overstepped minimum; reset STEPMIN
                STEPMIN = 1;
                delete_constr = 1;
            end
            X = X+STEPMIN*SD;
        else
            % Not a Newton step and hit a constraint: is_qp or LLS or maybe lp
            X = X+STEPMIN*SD;  
        end              
    else %  isempty(indf) | ~isfinite(STEPMIN)
        % did not hit a constraint
        if strcmp(dirType, NewtonStep)
            % Newton step and no constraint hit: LLS or maybe is_qp
            STEPMIN = 1;   % Exact distance to the solution. Now delete constr.
            X = X + SD;
            delete_constr = 1;
        else % Not a Newton step: is_qp or lp or LLS
            
            if (~is_qp && ~LLS) || strcmp(dirType, NegCurv) % LP or neg def (implies is_qp)
                % neg def -- unbounded
                if norm(SD) > errnorm
                    if normalize < 0
                        STEPMIN=abs((X(numberOfVariables)+1e-5)/(SD(numberOfVariables)+eps));
                    else 
                        STEPMIN = 1e16;
                    end
                    X=X+STEPMIN*SD;
                    how='unbounded'; 
                    % was exitflag = 5; 
                    exitflag = -1;
                else % norm(SD) <= errnorm
                    how = 'ill posed';
                    % was exitflag = 6; 
                    exitflag = -1;
                end
                ACTIND = indepInd(ACTIND);
                output.iterations = iterations;
                return
            else % singular: solve compatible system for a solution: is_qp or LLS
                if is_qp
                    projH = Z'*H*Z; 
                    Zgf = Z'*gf;
                    projSD = pinv(projH)*(-Zgf);
                else % LLS
                    projH = HZ'*HZ; 
                    Zgf = Z'*gf;
                    projSD = pinv(projH)*(-Zgf);
                end
                
                % Check if compatible
                if norm(projH*projSD+Zgf) > 10*eps*(norm(projH) + norm(Zgf))
                    % system is incompatible --> it's a "chute": use SD from compdir
                    % unbounded in SD direction
                    if norm(SD) > errnorm
                        if normalize < 0
                            STEPMIN=abs((X(numberOfVariables)+1e-5)/(SD(numberOfVariables)+eps));
                        else 
                            STEPMIN = 1e16;
                        end
                        X=X+STEPMIN*SD;
                        how='unbounded'; 
                        % was exitflag = 5;
                        exitflag = -1;
                    else % norm(SD) <= errnorm
                        how = 'ill posed';
                        %was exitflag = 6;
                        exitflag = -1;
                    end
                    ACTIND = indepInd(ACTIND);
                    output.iterations = iterations;
                    return
                else % Convex -- move to the minimum (compatible system)
                    SD = Z*projSD;
                    if gf'*SD > 0
                        SD = -SD;
                    end
                    dirType = 'singular';
                    % First check if constraint is violated.
                    GSD=A*SD;
                    indf = find((GSD > errnorm * norm(SD))  &  ~aix);
                    if isempty(indf) % No constraints to hit
                        STEPMIN=1;
                        delete_constr = 1;
                        dist=[]; ind2=[]; ind=[];
                    else % Find distance to the nearest constraint
                        dist = abs(cstr(indf)./GSD(indf));
                        [STEPMIN,ind2] =  min(dist);
                        ind2 = find(dist == STEPMIN);
                        % Bland's rule for anti-cycling: if there is more than one 
                        % blocking constraint then add the one with the smallest index.
                        ind=indf(min(ind2));
                    end
                    if STEPMIN > 1  % Overstepped minimum; reset STEPMIN
                        STEPMIN = 1;
                        delete_constr = 1;
                    end
                    X = X + STEPMIN*SD; 
                end
            end % if ~is_qp | smallRealEig < -eps
        end % if strcmp(dirType, NewtonStep)
    end % if ~isempty(indf)& isfinite(STEPMIN) % Hit a constraint
    
    %----Check if reached minimum in current subspace-----
    
    if delete_constr
        % Note: only reach here if a minimum in the current subspace found
        %       LP's do not enter here.
        if ACTCNT>0
            if is_qp
                rlambda = -R\(Q'*(H*X+f));
            elseif LLS
                rlambda = -R\(Q'*(H'*(H*X-f)));
                % else: lp does not reach this point
            end
            actlambda = rlambda;
            actlambda(eqix) = abs(rlambda(eqix));
            indlam = find(actlambda < 0);
            if (~~isempty(indlam)) 
                lambda(indepInd(ACTIND)) = normf * (rlambda./normA(ACTIND));
                ACTIND = indepInd(ACTIND);
                output.iterations = iterations;
                return
            end
            % Remove constraint
            lind = find(ACTIND == min(ACTIND(indlam)));
            lind = lind(1);
            ACTSET(lind,:) = [];
            aix(ACTIND(lind)) = 0;
            [Q,R]=qrdelete(Q,R,lind);
            ACTIND(lind) = [];
            ACTCNT = length(ACTIND);
            simplex_iter = 0;
            ind = 0;
        else % ACTCNT == 0
            output.iterations = iterations;
            return
        end
    end
    
    % If we are in the Phase-1 procedure check if the slack variable
    % is zero indicating we have found a feasible starting point.
    if normalize < 0
        if X(numberOfVariables,1) < eps
            ACTIND = indepInd(ACTIND);
            output.iterations = iterations;
            return;
        end
    end   
    
    % Calculate gradient w.r.t objective at this point
    if is_qp
        gf=H*X+f;
    elseif LLS % LLS
        gf=H'*(H*X-f);
        % else gf=f still true.
    end
    
    % Update constraints
    cstr = A*X-B;
    cstr(eqix) = abs(cstr(eqix));
    if max(cstr) > 1e5 * errnorm
        if max(cstr) > norm(X) * errnorm 
            how='unreliable'; 
            exitflag = -1;
        end
    end
    
    %----Add blocking constraint to working set----
    if ind % Hit a constraint
        aix(ind)=1;
        CIND = length(ACTIND) + 1;
        ACTSET(CIND,:)=A(ind,:);
        ACTIND(CIND)=ind;
        [m,n]=size(ACTSET);
        [Q,R] = qrinsert(Q,R,CIND,A(ind,:)');
        ACTCNT = length(ACTIND);
    end
    if ~simplex_iter
        [m,n]=size(ACTSET);
        Z = Q(:,m+1:n);
        if ACTCNT == numberOfVariables - 1, simplex_iter = 1; end
        oldind = 0; 
    else
        
        %---If Simplex Alg. choose leaving constraint---
        rlambda = -R\(Q'*gf);
        if isinf(rlambda(1)) && rlambda(1) < 0 
            fprintf('         Working set is singular; results may still be reliable.\n');
            [m,n] = size(ACTSET);
            rlambda = -(ACTSET + sqrt(eps)*randn(m,n))'\gf;
        end
        actlambda = rlambda;
        actlambda(eqix)=abs(actlambda(eqix));
        indlam = find(actlambda<0);
        if ~isempty(indlam)
            if STEPMIN > errnorm
                % If there is no chance of cycling then pick the constraint 
                % which causes the biggest reduction in the cost function. 
                % i.e the constraint with the most negative Lagrangian 
                % multiplier. Since the constraints are normalized this may 
                % result in less iterations.
                [minl,lind] = min(actlambda);
            else
                % Bland's rule for anti-cycling: if there is more than one 
                % negative Lagrangian multiplier then delete the constraint
                % with the smallest index in the active set.
                lind = find(ACTIND == min(ACTIND(indlam)));
            end
            lind = lind(1);
            ACTSET(lind,:) = [];
            aix(ACTIND(lind)) = 0;
            [Q,R]=qrdelete(Q,R,lind);
            Z = Q(:,numberOfVariables);
            oldind = ACTIND(lind);
            ACTIND(lind) = [];
            ACTCNT = length(ACTIND);
        else
            lambda(indepInd(ACTIND))= normf * (rlambda./normA(ACTIND));
            ACTIND = indepInd(ACTIND);
            output.iterations = iterations;
            return
        end
    end %if ACTCNT<numberOfVariables
    
    %----------Compute Search Direction-------------      
    if (is_qp)
        Zgf = Z'*gf; 
        if ~isempty(Zgf) && (norm(Zgf) < 1e-15)
            SD = zeros(numberOfVariables,1); 
            dirType = ZeroStep;
        else
            [SD, dirType] = compdir(Z,H,gf,numberOfVariables,f);
        end
    elseif (LLS)
        Zgf = Z'*gf;
        HZ = H*Z;
        if (norm(Zgf) < 1e-15)
            SD = zeros(numberOfVariables,1);
            dirType = ZeroStep;
        else
            HXf=H*X-f;
            gf=H'*(HXf);
            [mm,nn]=size(HZ);
            if mm >= nn
                [QHZ, RHZ] =  qr(HZ,0);
                Pd = QHZ'*HXf;
                % SD = - Z*(RHZ(1:nn, 1:nn) \ Pd(1:nn,:));
                % Now need to check which is dependent
                if min(size(RHZ))==1 % Make sure RHZ isn't a vector
                    depInd = find( abs(RHZ(1,1)) < tolDep);
                else
                    depInd = find( abs(diag(RHZ)) < tolDep );
                end  
            end
            if mm >= nn && isempty(depInd) % Newton step
                SD = - Z*(RHZ(1:nn, 1:nn) \ Pd(1:nn,:));
                dirType = NewtonStep;
            else % steepest descent direction
                SD = -Z*(Z'*gf);
                dirType = SteepDescent;
            end
        end
    else % LP
        if ~simplex_iter
            SD = -Z*(Z'*gf);
            gradsd = norm(SD);
        else
            gradsd = Z'*gf;
            if  (gradsd > 0)    SD = -Z;
            else                SD = Z;
            end
        end
        if abs(gradsd) < 1e-10 % Search direction null
            % Check whether any constraints can be deleted from active set.
            % rlambda = -ACTSET'\gf;
            if ~oldind
                rlambda = -R\(Q'*gf);
                ACTINDtmp = ACTIND; Qtmp = Q; Rtmp = R;
            else
                % Reinsert just deleted constraint.
                ACTINDtmp = ACTIND;
                ACTINDtmp(lind+1:ACTCNT+1) = ACTIND(lind:ACTCNT);
                ACTINDtmp(lind) = oldind;
                [Qtmp,Rtmp] = qrinsert(Q,R,lind,A(oldind,:)');
            end
            actlambda = rlambda;
            actlambda(1:neqcstr) = abs(actlambda(1:neqcstr));
            indlam = find(actlambda < errnorm);
            lambda(indepInd(ACTINDtmp)) = normf * (rlambda./normA(ACTINDtmp));
            if ~(~isempty(indlam))
                ACTIND = indepInd(ACTIND);
                output.iterations = iterations;
                return
            end
            cindmax = length(indlam);
            cindcnt = 0;
            m = length(ACTINDtmp);
            while (abs(gradsd) < 1e-10) && (cindcnt < cindmax)
                cindcnt = cindcnt + 1;
                lind = indlam(cindcnt);
                [Q,R]=qrdelete(Qtmp,Rtmp,lind);
                Z = Q(:,m:numberOfVariables);
                if m ~= numberOfVariables
                    SD = -Z*Z'*gf;
                    gradsd = norm(SD);
                else
                    gradsd = Z'*gf;
                    if  (gradsd > 0)    SD = -Z;
                    else                SD = Z;
                    end
                end
            end
            if abs(gradsd) < 1e-10  % Search direction still null
                ACTIND = indepInd(ACTIND);
                output.iterations = iterations;
                return;
            else
                ACTIND = ACTINDtmp;
                ACTIND(lind) = [];
                aix = zeros(ncstr,1);
                aix(ACTIND) = 1;
                ACTCNT = length(ACTIND);
                ACTSET = A(ACTIND,:);
            end
            lambda = zeros(ncstr,1);
        end
    end % if is_qp
end % while 

if iterations >= maxiter
    exitflag = 0;
    how = 'MaxSQPIter';
end

output.iterations = iterations;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[Q,R,A,B,X,Z,how,ACTSET,ACTIND,ACTCNT,aix,eqix,neqcstr,ncstr,remove,exitflag]= ...
    eqnsolv(A,B,eqix,neqcstr,ncstr,numberOfVariables,LLS,H,X,f,normf, ...
    normA,verbosity,aix,ACTSET,ACTIND,ACTCNT,how,exitflag)
% EQNSOLV Helper function for QPSUB.
%    Checks whether the working set is linearly independent and
%    finds a feasible point with respect to the working set constraints.
%    If the equalities are dependent but not consistent, warning
%    messages are given. If the equalities are dependent but consistent, 
%    the redundant constraints are removed and the corresponding variables 
%    adjusted.

% set tolerances
tolDep = 100*numberOfVariables*eps;      

remove =[];

% First see if the equality constraints form a consistent system.
[Qa,Ra]=qr(A(eqix,:));

% Form vector of dependent indices.
if min(size(Ra))==1 % Make sure Ra isn't a vector
    depInd = find( abs(Ra(1,1)) < tolDep);
else
    depInd = find( abs(diag(Ra)) < tolDep );
end
if neqcstr > numberOfVariables
    depInd = [depInd; ((numberOfVariables+1):neqcstr)'];
end      

if ~isempty(depInd)    % equality constraints are dependent
    how='dependent';
    exitflag = 1;
    bdepInd =  abs(Qa(:,depInd)'*B(eqix)) >= tolDep ;
    
    if any( bdepInd ) % Not consistent
        how='infeasible';   
        exitflag = -1;
    else % the equality constraints are consistent
        % Delete the redundant constraints
        % By QR factoring the transpose, we see which columns of A' (rows of A) move to the end
        [Qat,Rat,Eat]=qr(A(eqix,:)');        
        [i,j] = find(Eat); % Eat permutes the columns of A' (rows of A)
        remove = i(depInd);
        numDepend = nnz(remove);
        A(eqix(remove),:)=[];
        B(eqix(remove))=[];
        neqcstr = neqcstr - numDepend;
        ncstr = ncstr - numDepend;
        eqix = 1:neqcstr;
        aix(remove) = [];
        ACTIND(1:numDepend) = [];
        ACTIND = ACTIND - numDepend;      
        ACTSET = A(ACTIND,:);
        ACTCNT = ACTCNT - numDepend;
    end % consistency check
end % dependency check

% Now that we have done all we can to make the equality constraints
% consistent and independent we will check the inequality constraints
% in the working set.  First we want to make sure that the number of 
% constraints in the working set is only greater than or equal to the
% number of variables if the number of (non-redundant) equality 
% constraints is greater than or equal to the number of variables.
if ACTCNT >= numberOfVariables
    ACTCNT = max(neqcstr, numberOfVariables-1);
    ACTIND = ACTIND(1:ACTCNT);
    ACTSET = A(ACTIND,:);
    aix = zeros(ncstr,1);
    aix(ACTIND) = 1;
end

% Now check to see that all the constraints in the working set are linearly independent.
if ACTCNT > neqcstr
    [Qat,Rat,Eat]=qr(ACTSET');
    
    % Form vector of dependent indices.
    if min(size(Rat))==1            % Make sure Rat isn't a vector
        depInd = find( abs(Rat(1,1)) < tolDep);
    else
        depInd = find( abs(diag(Rat)) < tolDep );
    end
    
    if ~isempty(depInd)
        [i,j] = find(Eat);          % Eat permutes the columns of A' (rows of A)
        remove2 = i(depInd);
        removeEq   = remove2(remove2 <= neqcstr);
        removeIneq = remove2(remove2 > neqcstr);
        
        if (~isempty(removeEq))     % Just take equalities as initial working set.
            ACTIND = 1:neqcstr; 
        else                        % Remove dependent inequality constraints.
            ACTIND(removeIneq) = [];
        end
        aix = zeros(ncstr,1);
        aix(ACTIND) = 1;
        ACTSET = A(ACTIND,:);
        ACTCNT = length(ACTIND);
    end  
end

[Q,R]=qr(ACTSET');
Z = Q(:,ACTCNT+1:numberOfVariables);

if ~strcmp(how,'infeasible') && ACTCNT > 0
    % Find point closest to the given initial X which satisfies working set constraints.
    minnormstep = Q(:,1:ACTCNT) * ...
        ((R(1:ACTCNT,1:ACTCNT)') \ (B(ACTIND) - ACTSET*X));
    X = X + minnormstep; 
    % Sometimes the "basic" solution satisfies Aeq*x= Beq 
    % and A*X < B better than the minnorm solution. Choose the one
    % that the minimizes the max constraint violation.
    err = A*X - B;
    err(eqix) = abs(err(eqix));
    if any(err > eps)
        Xbasic = ACTSET\B(ACTIND);
        errbasic = A*Xbasic - B;
        errbasic(eqix) = abs(errbasic(eqix));
        if (max(errbasic) < max(err))   X = Xbasic;     end
    end
end

% --------------------------------------------------------------------------------
function [SD, dirType] = compdir(Z,H,gf,nvars,f)
% COMPDIR Computes a search direction in a subspace defined by Z.  [SD,
% dirType] = compdir(Z,H,gf,nvars,f) returns a search direction for the
% subproblem 0.5*Z'*H*Z + Z'*gf. Helper function for NLCONST. SD is Newton
% direction if possible. SD is a direction of negative curvature if the
% Cholesky factorization of Z'*H*Z fails. If the negative curvature
% direction isn't negative "enough", SD is the steepest descent direction.
% For singular Z'*H*Z, SD is the steepest descent direction even if small,
% or even zero, magnitude.

%   Copyright 1990-2003 The MathWorks, Inc.

% Define constant strings
Newton = 'Newton';                     % Z'*H*Z positive definite
NegCurv = 'negative curvature chol';   % Z'*H*Z indefinite
SteepDescent = 'steepest descent';     % Z'*H*Z (nearly) singular

% Compute the projected Newton direction if possible
projH = Z'*H*Z;
[R, p] = chol(projH);
if ~p  % positive definite: use Newton direction
    SD = - Z*(R \ ( R'\(Z'*gf)));
    dirType = Newton;
else % not positive definite
    [L,sneg] = choltrap(projH);
    if ~isempty(sneg) && sneg'*projH*sneg < -sqrt(eps) % if negative enough
        SD = Z*sneg;
        dirType = NegCurv;
    else % Not positive definite, not negative definite "enough" so use steepest descent direction
        stpDesc = - Z*(Z'*gf);
        % ||SD|| may be (close to) zero, but qpsub handles that case
        SD = stpDesc;
        dirType = SteepDescent;
    end %   
end % ~p  (positive definite)

% Make sure it is a descent direction
if gf'*SD > 0
    SD = -SD;
end

%-----------------------------------------------
function [L,sneg] = choltrap(A)
% CHOLTRAP Compute Cholesky factor or direction of negative curvature.
%     [L, SNEG] = CHOLTRAP(A) computes the Cholesky factor L, such that
%     L*L'= A, if it exists, or returns a direction of negative curvature
%     SNEG for matrix A when A is not positive definite. If A is positive
%     definite, SNEG will be []. 
%
%     If A is singular, it is possible that SNEG will not be a direction of
%     negative curvature (but will be nonempty). In particular, if A is
%     positive semi-definite, SNEG will be nonempty but not a direction of
%     negative curvature. If A is indefinite but singular, SNEG may or may
%     not be a direction of negative curvature.

sneg = [];
n = size(A,1);
L = eye(n);
tol = 0;    % Dividing by sqrt of small number isn't a problem 
for k=1:n-1
    if A(k,k) <= tol
        elem = zeros(length(A),1); 
        elem(k,1) = 1;
        sneg = L'\elem;
        return;
    else
        L(k,k) = sqrt(A(k,k));
        s = k+1:n;
        L(s,k) = A(s,k)/L(k,k);
        A(k+1:n,k+1:n) =  A(k+1:n,k+1:n)  - tril(L(k+1:n,k)*L(k+1:n,k)');   
    end
end
if A(n,n) <= tol
    elem = zeros(length(A),1); 
    elem(n,1) = 1;
    sneg = L'\elem;
else
    L(n,n) = sqrt(A(n,n));
end

%------------------------------------------------------------------
function value = optimgetfast(options,name,defaultopt)
%OPTIMGETFAST Get OPTIM OPTIONS parameter with no error checking so fast.

if (~isempty(options))      value = options.(name);
else                        value = [];
end

if (isempty(value))    value = defaultopt.(name);   end

% -----------------------------------------------------------------------
function [allfcns,msg] = optimfcnchk1(funstr,caller,lenVarIn,gradflag,hessflag,constrflag)
% OPTIMFCNCHK Pre- and post-process function expression for FUNCHK.

%   Copyright 1990-2002 The MathWorks, Inc. 
%   $Revision: 1.7 $  $Date: 2002/03/12 20:36:19 $

% Initialize
if nargin < 6
	if nargin < 5
		hessflag = 0;
		if nargin < 4
			gradflag = 0;
		end
	end
end

allfcns = {};
gradfcn = [];
hessfcn = [];
calltype = 'fun';

% {fun}
if (isa(funstr, 'cell') && length(funstr)==1)    % take the cellarray apart: we know it is nonempty
    [funfcn, msg] = fcnchk1(funstr{1},lenVarIn);
    if (~isempty(msg))      error(msg);    end
    % {fun,[]}      
elseif isa(funstr, 'cell') && length(funstr)==2 && isempty(funstr{2})
    [funfcn, msg] = fcnchk1(funstr{1},lenVarIn);
    if (~isempty(msg))      error(msg);    end
    % {fun, grad}   
elseif isa(funstr, 'cell') && length(funstr)==2 % and ~isempty(funstr{2})
    
    [funfcn, msg] = fcnchk1(funstr{1},lenVarIn);
    if (~isempty(msg))      error(msg);    end
    [gradfcn, msg] = fcnchk1(funstr{2},lenVarIn);
    if (~isempty(msg))      error(msg);    end
    calltype = 'fun_then_grad';
    if (~gradflag)
        calltype = 'fun';
    end
    % {fun, [], []}   
elseif isa(funstr, 'cell') && length(funstr)==3 && ~isempty(funstr{1}) && isempty(funstr{2}) && isempty(funstr{3})
    [funfcn, msg] = fcnchk1(funstr{1},lenVarIn);
    if (~isempty(msg))      error(msg);    end
    % {fun, grad, hess}   
elseif isa(funstr, 'cell') && length(funstr)==3 && ~isempty(funstr{2}) && ~isempty(funstr{3})
    [funfcn, msg] = fcnchk1(funstr{1},lenVarIn);
    if (~isempty(msg))      error(msg);    end
    [gradfcn, msg] = fcnchk1(funstr{2},lenVarIn);
    if (~isempty(msg))      error(msg);    end
    [hessfcn, msg] = fcnchk1(funstr{3},lenVarIn);
    if (~isempty(msg))      error(msg);    end
    calltype = 'fun_then_grad_then_hess';
    if ~hessflag && ~gradflag
        calltype = 'fun';
    elseif hessflag && ~gradflag
        calltype = 'fun';
    elseif ~hessflag && gradflag
        calltype = 'fun_then_grad';
    end
    
    % {fun, grad, []}   
elseif isa(funstr, 'cell') && length(funstr)==3 && ~isempty(funstr{2}) && isempty(funstr{3})
    [funfcn, msg] = fcnchk1(funstr{1},lenVarIn);
    if (~isempty(msg))      error(msg);    end
    [gradfcn, msg] = fcnchk1(funstr{2},lenVarIn);
    if (~isempty(msg))      error(msg);    end
    calltype = 'fun_then_grad';
    if ~gradflag
        calltype = 'fun';
    end
    
elseif ~isa(funstr, 'cell')  %Not a cell; is a string expression, function name string or inline object
    [funfcn, msg] = fcnchk1(funstr,lenVarIn);
    if (~isempty(msg))      error(msg);    end
else
    errmsg = sprintf('%s\n%s', ...
        'FUN must be a function or an inline object;', ...
        ' or, FUN may be a cell array that contains these type of objects.');
    error(errmsg)
end

allfcns{1} = calltype;      allfcns{2} = caller;
allfcns{3} = funfcn;        allfcns{4} = gradfcn;
allfcns{5} = hessfcn;

% ----------------------------------------------------
function [f,msg] = fcnchk1(fun,varargin)
msg = '';
if isstr(fun)
    f = strtrim(fun);
elseif isa(fun,'function_handle') 
    f = fun; 
else
    f = '';
    msg = ['FUN must be a function, a valid string expression, ', ...
            sprintf('\n'),'or an inline function object.'];
end

%------------------------------------------
function s1 = strtrim(s)
%STRTRIM Trim spaces from string.
if isempty(s)
    s1 = s;
else            % remove leading and trailing blanks (including nulls)
    c = find(s ~= ' ' & s ~= 0);
    s1 = s(min(c):max(c));
end

