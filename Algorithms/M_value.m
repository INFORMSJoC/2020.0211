function [M,time] = M_value(A_y,b_y,Q,x_opt,Lower_bound,Best_solution,time_limit,tolerr)
    M=sdpvar;
    F=[M>=0];
    n_y=size(A_y,2);
    m_y=size(A_y,1);
    y=sdpvar(n_y,1);
    z=sdpvar(m_y,1);
    W=[A_y*y>= b_y, y>=0, A_y'*z<=Q'*x_opt, z>=0, b_y'*z>=Lower_bound, b_y'*z <= Best_solution, uncertain([y' z'])];
    F=F+[A_y*y-b_y<=ones(m_y,1)*M];
    F=F+[y<=ones(n_y,1)*M];
    F=F+[Q'*x_opt-A_y'*z<=ones(n_y,1)*M];
    F=F+[z<=ones(m_y,1)*M];
    sol_M=optimize(F+W,M,sdpsettings('solver','gurobi','verbose', 0,'gurobi.TimeLimit',time_limit,'gurobi.MIPGap',tolerr,'savesolveroutput',1));
    if sol_M.problem==12
        M=1000000;
    else
        M=value(M);
    end
    time=sol_M.solvertime;
end

