clear
close all
clc

X0 = [0 3]';
X_DFP = DFP(X0, 5);
X_BFGS = BFGS(X0, 5);

function Cost_Value = Cost(X) % Insert Cost Function Here
x1 = X(1); x2 = X(2);

Cost_Value = (x1 - 2)^4 + (x1 - 3*x2)^2;
end

function Jacibian_Value = Jacibian(X) % Insert Jacobian of the Cost Function Here
x1 = X(1); x2 = X(2);

Jacibian_Value = [4*(x1 - 2)^3 + 2*(x1 - 3*x2)
    -6*(x1 - 3*x2)];
end

function x_Optimal = GoldenSection(xl, xu, xm, iteration_type, number_of_iterations, X0, dX)
a = (5^0.5 - 1)/2;

for i = 1:number_of_iterations
    if iteration_type == 0 % first iteration
        xb = xl + a*(xu - xl);
        xa = xl + (1 - a)*(xu - xl);
    elseif iteration_type == 1 % xu is deleted
        xb = xm;
        xa = xl + (1 - a)*(xu - xl);
    else iteration_type == 2 % xl is deleted
        xb = xl + a*(xu - xl);
        xa = xm;
    end
    
    fa = Cost(X0 + xa*dX);
    fb = Cost(X0 + xb*dX);
    
    if fa >= fb
        xl = xa;
        xm = xb;
        iteration_type = 2;
    else
        xu = xb;
        xm = xa;
        iteration_type = 1;
    end
end
x_Optimal = xm;
end

function X_Optimal = DFP(X0, number_of_iterations)
L = length(X0);
H = eye(L);

for i = 1 : number_of_iterations
    C = Jacibian(X0);
    D = -H*C;
    alpha = GoldenSection(0, 1, 0, 0, 10, X0, D);
    X1 = X0 + alpha*D;
    P = alpha*D;
    C1 = Jacibian(X1);
    Q = C1 - C;
    H1 = H + (P*P')/(P'*Q) - (H*Q*Q'*H)/(Q'*H*Q);
    
    H = H1;
    X0 = X1;
end
X_Optimal = X0;
end

function  X_Optimal = BFGS(X0, number_of_iterations)
L = length(X0);
B = eye(L);

for i = 1 : number_of_iterations
    C = Jacibian(X0);
    D = -B\C;
    alpha = GoldenSection(0, 1, 0, 0, 10, X0, D);
    X1 = X0 + alpha*D;
    P = alpha*D;
    C1 = Jacibian(X1);
    Q = C1 - C;
    B1 = B + (Q*Q')/(Q'*P) - (B*P*P'*B)/(P'*B*P);
    
    B = B1;
    X0 = X1;
end
X_Optimal = X0;
end