% Define global variables to pass to a residual function 
function test2dgau()
global X Y Zdata;

[X,Y] = meshgrid(-3:.05:3 , -3:.05:3);
% Define a specific 2D function with 3 parameters to test "fminsearch"
A = .5;
sigx = .33;
sigy = 1.3;
Zdata = A*exp(-X.^2/(2*sigx^2)).*exp(-Y.^2/(2*sigy^2));
% plot (X,Y,Zdata);

% Attempt to invoke fminsearch
param_guess = [1 1 1];
bestparams = fminsearch(@residual,param_guess);
disp ('bestparms');
bestparams

% Define a residual function
function err = residual(par)
  global X Y Zdata;
  Zmodel = par(1)*exp(-X.^2/(2*par(2)^2)).*exp(-Y.^2/(2*par(3)^2));
  err = sum(sum( (Zmodel-Zdata).^2 ));

