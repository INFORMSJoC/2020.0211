%PARTIAL R4B X = X1

%------------------------ calculate upper bound ---------------------------
maxtime = 1000;
y = sdpvar(Jvert*K,1); 
W = sdpvar(n,Jvert*K, 'full');
x = sdpvar(n,1);

%constraints c,d,f and y >= 0
constraints = [sum(reshape(y,Jvert,[])) == 1, A*W - b*y' <= 0, y >= 0, W(:) >= 0];

%constraints g
for k = 1:K-1
    constraints = [constraints, sum(W(:,[1+(k-1)*Jvert:Jvert+(k-1)*Jvert]),2) == sum(W(:,[1+k*Jvert:Jvert+k*Jvert]),2)];
end
constraints = [constraints, sum(W(:,[1:Jvert]),2) == x];

obj = -trace(Q*W) - c'*y
ops = sdpsettings('solver','gurobi','savesolveroutput',2,'verbose',1,'debug',1);
sol = optimize(constraints,obj,ops);
optobj = value(-obj)
ubx1 = optobj
time_ubx1 = sol.solveroutput.result.runtime

%------------------------- calculate lower bound --------------------------
x = value(x);
xvalue_ubx1 = x;
y = value(y);
yvalue_ubx1 = y;
W = value(W);
Wvalue_ubx1 = W;

tau = x'*Q'*y + c'*y;
eps = optobj - tau

if value(eps) > 0.0001
    X = [xvalue_ubx1];
for j = 1:Jvert*K
    if y(j) <= 0.0001
        X = X;
    else
        X = [X,W(:,j)/y(j)];
    end
end

X = value(X);
XX = [];
YY = [];
l = 0;
maxobjx1 = [];
time1 = [];
time2 = zeros(size(X,2),1);
for ii = 1:size(X,2)
    y = sdpvar(Jvert*K,1);
    constraints = [y >= 0, sum(reshape(y,Jvert,[])) == 1];
    obj = - X(:,ii)'*Q'*y - c'*y;
    ops = sdpsettings('solver','gurobi','savesolveroutput',2,'verbose',1,'debug',1);
    sol = optimize(constraints,obj,ops);
    a = value(-obj);
    time1(ii) = sol.solveroutput.result.runtime;
    
    x = X(:,ii);
    y = value(y);
    delta = 1;
    val = 0;
    
    while(delta >= 0.0001)
    tau = x'*Q'*y + c'*y
    
    x = sdpvar(n,1);
    t = sdpvar(n,1);
    
    constraints = [A*x <= b, x>= 0];
    obj = -x'*Q'*y - c'*y;
    ops = sdpsettings('solver','gurobi','savesolveroutput',2,'verbose',1,'debug',1);
    sol = optimize(constraints,obj,ops);
    x = value(x);
    time2(ii) = [time2(ii) + sol.solveroutput.result.runtime];
    
    y = sdpvar(Jvert*K,1)
    constraints = [y >= 0, sum(reshape(y,Jvert,[])) == 1];
    obj = - x'*Q'*y - c'*y;
    ops = sdpsettings('solver','gurobi','savesolveroutput',2,'verbose',1,'debug',1);
    sol = optimize(constraints,obj,ops);
    value(-obj)
    y = value(y);   
    delta = value(-obj) - tau
    time2(ii) = [time2(ii) + sol.solveroutput.result.runtime];
    end
    
XX = [XX,x];
YY = [YY,y];
maxobjx1 = [maxobjx1,-obj];
l = l + 1
end
maxobjx1 = value(maxobjx1);
[tau,index] = max(maxobjx1)
time_lbx1 = max(time1) + max(time2)
xvalue_lbx1 = value(XX(:,index));
yvalue_lbx1 = value(YY(:,index));
else
    time_lbx1 = 0;
    xvalue_lbx1 = xvalue_ubx1;
end

lbx1 = tau;
eps = optobj - tau
time_partialR4B_x1 = time_lbx1 + time_ubx1
%%
%Initial solution for MBCO X = X1
thetamax = [];
idxmax = [];
z0 = zeros(K,Jvert);
for k =  1:K
        [thetam,idx] = max(Q([I(:,k)],:)*xvalue_lbx1 + c(I(:,k)));
        thetamax(k) = thetam;
        idxmax(k) = idx;
        z0(k,idx) = 1;
end
thetamax = thetamax';
%%
%MIXED BINARY CONVEX OPTIMIZATION PROBLEM (MBCO) IN GUROBI X = X1
if eps > 0.0001
    maxtime = 1000 - time_partialR4B_x1;

    x = sdpvar(n,1);
    theta = sdpvar(K,1);
    z = binvar(K,Jvert,'full');

    constraints = [x >= 0, A*x <= b];

    for k = 1:K
        for j = 1:Jvert
        constraints = [constraints, theta(k) <= Q([I(j,k)],:)*x + c(I(j,k)) + M*(1 - z(k,j))];
        end
    end

    constraints = [constraints, sum(z,2) == 1];

    constraints = [constraints, -sum(theta) >= -ubx1];
    assign(x,xvalue_lbx1)
    assign(theta,thetamax)
    assign(z,z0)
    obj = -sum(theta);
    ops = sdpsettings('solver','gurobi','gurobi.timelimit',maxtime,'savesolveroutput',2,'verbose',1,'debug',1);
    sol = optimize(constraints,obj,ops);
    sol.solveroutput.result.runtime

    R4B_lbx1 = value(-obj)
else
    R4B_lbx1 = lbx1
end
%%
%PARTIAL R4B X = X2

%--------------------------- calculate upper bound ------------------------
y = sdpvar(Jvert*K,1); 
W = sdpvar(n,Jvert*K,'full');
z = sdpvar(n,Jvert*K, 'full');
u = sdpvar(n,1);
x = sdpvar(n,1);

%constraints c,d,f and y >= 0
constraints = [sum(reshape(y,Jvert,[])) == 1, A*W - b*y' <= 0, y >= 0, W(:) >= 0];

%constraints b,e
constraints = [constraints, sum(u) <= 1, y >= sum(z)'];
for i = 1:n
   constraints = [constraints, expcone([x(i) - Gamma, 1, u(i)])];
   for j = 1:Jvert*K
       constraints = [constraints, expcone([W(i,j)-y(j)*Gamma,y(j),z(i,j)])];
   end
end

%constraints g
for k = 1:K-1
    constraints = [constraints, sum(W(:,[1+(k-1)*Jvert:Jvert+(k-1)*Jvert]),2) == sum(W(:,[1+k*Jvert:Jvert+k*Jvert]),2)];
end
constraints = [constraints, sum(W(:,[1:Jvert]),2) == x];

obj = -trace(Q*W) - c'*y
ops = sdpsettings('solver','mosek','savesolveroutput',2,'verbose',1,'debug',1);
sol = optimize(constraints,obj,ops);
optobj = value(-obj)
ubx2 = optobj

time_ubx2 = sol.solveroutput.res.info.MSK_DINF_OPTIMIZER_TIME
ubx2 = optobj;

%------------------- calculate lower bound --------------------------------

x = value(x);
y = value(y);
W = value(W);
xvalue_ubx2 = x;
yvalue_ubx2 = y;
Wvalue_ubx2 = W;

tau = x'*Q'*y + c'*y;
eps = optobj - tau

if eps > 0.0001
    X = [xvalue_ubx2];
    y = yvalue_ubx2;
    W = Wvalue_ubx2;
for j = 1:Jvert*K
    if y(j) <= 0.0001
        X = X;
    else
        X = [X,W(:,j)/y(j)];
    end
end

X = value(X);
XX = [];
YY = [];
maxobj = [];
time1 = [];
time2 = zeros(size(X,2),1);
l = 0;
for ii = 1:size(X,2)
    y = sdpvar(Jvert*K,1);
    constraints = [y >= 0, sum(reshape(y,Jvert,[])) == 1];
    obj = - X(:,ii)'*Q'*y - c'*y;
    ops = sdpsettings('solver','mosek','savesolveroutput',2,'verbose',1,'debug',1);
    sol = optimize(constraints,obj,ops);
    a = value(-obj);
    time1(ii) = sol.solveroutput.res.info.MSK_DINF_OPTIMIZER_TIME;
    
    x = X(:,ii);
    y = value(y);
    delta = 1;
    val = 0;
    
    while(delta >= 0.0001)
    tau = x'*Q'*y + c'*y;
    
    x = sdpvar(n,1);
    t = sdpvar(n,1);
    
    constraints = [A*x <= b, x>= 0];
    constraints = [constraints, sum(t) <= 1];
    for i = 1:n
   constraints = [constraints, expcone([x(i) - Gamma, 1, t(i)])];
    end

    obj = -x'*Q'*y - c'*y
    ops = sdpsettings('solver','mosek','savesolveroutput',2,'verbose',1,'debug',1);
    sol = optimize(constraints,obj,ops);
    x = value(x);
    time2(ii) = time2(ii) + sol.solveroutput.res.info.MSK_DINF_OPTIMIZER_TIME;
    
    y = sdpvar(Jvert*K,1)
    constraints = [y >= 0, sum(reshape(y,Jvert,[])) == 1];
    obj = - x'*Q'*y - c'*y;
    ops = sdpsettings('solver','mosek','savesolveroutput',2,'verbose',1,'debug',1);
    sol = optimize(constraints,obj,ops);
    value(-obj)
    y = value(y);   
    delta = value(-obj) - tau
    time2(ii) = time2(ii) + sol.solveroutput.res.info.MSK_DINF_OPTIMIZER_TIME;
    end
    
XX = [XX,x];
YY = [YY,y];
maxobj = [maxobj,-obj];
l = l + 1
end
maxobj = value(maxobj);
[tau,index] = max(maxobj)
time_lbx2 = max(time1) + max(time2)
xvalue_lbx2 = value(XX(:,index));
yvalue_lbx2 = value(YY(:,index));
else 
    time_lbx2 = 0;
end

eps = optobj - tau
lbx2 = tau;
time_partialR4B_x2 = time_lbx2 + time_ubx2
%%
%MBCO IN MOSEK X = X2
if eps > 0.0001
    maxtime = 1000 - time_partialR4B_x2
    x = sdpvar(n,1);
    theta = sdpvar(K,1);
    z = binvar(K,Jvert,'full');
    t = sdpvar(n,1);

    constraints = [x >= 0, A*x <= b];

    for k = 1:K
        for j = 1:Jvert
        constraints = [constraints, theta(k) <= Q([I(j,k)],:)*x + c(I(j,k)) + M*(1 - z(k,j))];
        end
    end

    constraints = [constraints, sum(t) <= 1];
    for i = 1:n
       constraints = [constraints, expcone([x(i) - Gamma, 1, t(i)])];
    end

    constraints = [constraints, sum(z,2) == 1];

    constraints = [constraints, -sum(theta) >= -ubx2];
    obj = -sum(theta);
    ops = sdpsettings('solver','mosek','mosek.MSK_DPAR_OPTIMIZER_MAX_TIME',maxtime,'verbose',1,'debug',1);
    sol = optimize(constraints,obj,ops);

    R4B_lbx2 = value(-obj)
else
    R4B_lbx2 = lbx2
end

