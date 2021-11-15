function [val,time,sol_x,sol_y,obj] = R4B(Q,A_x,b_x,A_y,b_y,tolerr,time_limit)
% R4B is a code written by Ahmadreza Marandi and Jianzhe Zhen
% based on the paper: 
% Disjoint bilinear optimization: A tow-stage adjustable robust optimization
% perspective by Jianzhe Zhen, Ahmadreza Marandi, Danique de Moor, Dick den Hertog, and 
% Lieven Vandenberghe. 
% 
% R4B function solves the following bilinear problem 
%           min x'*Q*y
%           s.t. A_x*x>=b_x 
%                A_y*y>=b_y
%                y>=0, x>=0
% NOTICE: b_x and b_y are column vectors!
% tolerr is the allowed error gap between lower and upper bound
% time_limit is the time restriction on the solver to be terminated
%
%
%
%The outputs:
%val:       lower bound on the objective value after termination
%time:      solver time
%sol_x      x solution after termination
%sol_y      y solution after termination
%obj        objective value of the solution after termination
%
%
%Requirements: Yalmip, Gurobi

    limit_constraint=Inf;
    current_time=0;
    sdp=sdpsettings;
    optionGUROBI=sdpsettings('solver','gurobi','verbose', 0);
    val=-inf;  
    obj=+inf;
    sol_x_LDR=[];
    sol_y_LDR=[];
    obj_LDR=+inf;
    x_opt = [];
    y_opt = [];
    time=[];
    time_LDR=0;
    n_x=size(A_x,2);
    n_y=size(A_y,2);
    m_x=size(A_x,1);
    m_y=size(A_y,1);
    fprintf('\n R4B is using LDR to approximate the solution of the bilinear problem: \n\n');
    
      %% McCormick form
    % We first solve the relaxation of the problem as it has the least
    % complexity
    x=sdpvar(n_x,1);
    y=sdpvar(n_y,1);
    gama=sdpvar(n_x,n_y,'full');
    C=[];
     for k=1:m_y
        C=C+[A_y(k,:)*gama'-b_y(k)*x'>=0];
        C+[A_y(k,:)*(gama)'*(A_x)'+b_y(k)*(b_x)'-A_y(k,:)*y*(b_x)'-b_y(k)*(A_x*x)'];
    end
    for k=1:m_x
        C=C+[A_x(k,:)*gama-b_x(k)*y'];
    end
    C=C+[A_x*x>=b_x,(A_y*y>=b_y):'z0',x>=0,y>=0,gama(:)>=0];
    sol=optimize(C,sum(sum(Q.*gama)),sdpsettings('solver','gurobi','verbose', 0,'gurobi.TimeLimit',time_limit,'gurobi.MIPGap',tolerr));
    val=sum(sum(Q.*value(gama)));
    x_opt=value(x);
    x_Mc=x_opt;
    y_Mc=y_opt;
    y_opt=value(y);
    obj=x_opt'*Q*y_opt;
    if sol.problem==2
         lin_val=-Inf;
         val=-Inf;
    else
       lin_val=sum(sum(Q.*value(gama)));
    end 
    time_LDR=time_LDR+sol.solvertime;
    z0=dual(C('z0'));
    z1=zeros(m_y,n_x);
    for k=1:m_y
        z1(k,:)=dual(C(k))';
    end
     %% MCP 
    y_GBD=y_opt;
    x_GBD=x_opt;
    z_opt=z0;
    GBD_val=-Inf;
    count_GBD=0;
    while GBD_val<obj
        y=sdpvar(n_y,1);
        sdp.gurobi.TimeLimit=max([time_limit-time_LDR,0]);
        C=[(A_y*y>=b_y):'z',y>=0];
        sol=optimize(C,x_GBD'*Q*y,sdpsettings('solver','gurobi','verbose', 0,'gurobi.TimeLimit',max([time_limit-time_LDR,0]),'gurobi.MIPGap',tolerr));
        time_LDR=time_LDR+sol.solvertime;
        y_GBD=value(y);
        if x_GBD'*Q*y_GBD<=obj
            x_opt=x_GBD;
            y_opt=y_GBD;
            z_opt=dual(C('z'));
            obj=x_GBD'*Q*y_GBD;
            Best_solution=obj;
            Lower_bound=val;
            fprintf('Time=%f      lower bound=%f      best value=%f\n',time_LDR,Lower_bound,Best_solution)
        end
        if abs(obj-val)>tolerr
                y=sdpvar(n_y,1);
                sdp.gurobi.TimeLimit=max([time_limit-time_LDR,0]);
                sol=optimize([A_y*y>=b_y,y>=0],x_GBD'*Q*y,sdpsettings('solver','gurobi','verbose', 0,'gurobi.TimeLimit',max([time_limit-time_LDR,0])));
                time_LDR=time_LDR+sol.solvertime;
                y_GBD=value(y);
                GBD_val=x_GBD'*Q*y_GBD;
         else
                break;
         end
    end
    Best_solution=obj;
    Lower_bound=val;
    fprintf('Time=%f      lower bound=%f      best value=%f\n',time_LDR,Lower_bound,Best_solution)
    obj_LDR= Lower_bound;           
     if sol.problem ~= 1    % Infeasibility Check
        val=lin_val;        

     %% McCormick
        gama_Mc=value(gama);
        ind_y=find(y_Mc);
        if size(ind_y)~=0
            gama_s=gama_Mc(:,ind_y)./repmat(y_Mc(ind_y)',[n_x,1]);
            size_x=size(ind_y);
            for i=1:size_x
                x=gama_s(:,i);
                if abs(obj-val)>tolerr && sum(A_x*x-b_x>-10^(-5))==size(b_x,1) 
                    y=sdpvar(n_y,1);
                    sdp.gurobi.TimeLimit=max([time_limit-time_LDR,0]);
                    sol=optimize([A_y*y>=b_y,y>=0],x'*Q*y,sdpsettings('solver','gurobi','verbose', 0,'gurobi.TimeLimit',max([time_limit-time_LDR,0])));
                    time_LDR=time_LDR+sol.solvertime;
                    y=value(y);
                    if x'*Q*y<obj
                        y_opt=y;
                        x_opt=x;
                        obj=x'*Q*y;
                        %% GBD 
                        y_GBD=y;
                        GBD_val=-Inf;
                        count_GBD=0;
                        while GBD_val<obj
                            x=sdpvar(n_x,1);
                            sdp.gurobi.TimeLimit=max([time_limit-time_LDR,0]);
                            sol=optimize([A_x*x>=b_x,x>=0],x'*Q*y_GBD,sdpsettings('solver','gurobi','verbose', 0,'gurobi.TimeLimit',max([time_limit-time_LDR,0])));
                            time_LDR=time_LDR+sol.solvertime;
                            x_GBD=value(x);
                            if x_GBD'*Q*y_GBD<obj
                                x_opt=x_GBD;
                                y_opt=y_GBD;
                                obj=x_GBD'*Q*y_opt;
                                Best_solution=obj;
                                Lower_bound=val;
                                fprintf('Time=%f      lower bound=%f      best value=%f\n',time_LDR,Lower_bound,Best_solution)
                            end
                            if abs(obj-val)>tolerr
                                    y=sdpvar(n_y,1);
                %                    sdp.cplex.timelimit=max([time_limit-time_LDR,0]);
                                    sdp.gurobi.TimeLimit=max([time_limit-time_LDR,0]);
                                    sol=optimize([A_y*y>=b_y,y>=0],x_GBD'*Q*y,sdpsettings('solver','gurobi','verbose', 0,'gurobi.TimeLimit',max([time_limit-time_LDR,0])));
                                    time_LDR=time_LDR+sol.solvertime;
                                    y_GBD=value(y);
                                    GBD_val=x_GBD'*Q*y_GBD;
                             else
                                    break;
                             end
                        end
                        Best_solution=obj;
                        Lower_bound=val;
                        fprintf('Time=%f      lower bound=%f      best value=%f\n',time_LDR,Lower_bound,Best_solution)
                    end
                end
            end
        end
     end
    if abs(val-obj)<=tolerr
        Best_solution=obj;
        Lower_bound=val;
        fprintf('Time=%f      lower bound=%f      best value=%f\n',time_LDR,Lower_bound,Best_solution)
        fprintf('\n The problem is solved to OPTIMALITY using LDR. \n');
        time=[time_LDR,0,0,-1];
        sol_x=x_opt;
        sol_y=y_opt;
    else
        sol_x=x_opt;
        sol_y=y_opt;
        time=[time_LDR,0];
        fprintf('\n The problem is APPROXIMATED using LDR. \n');
        fprintf('\n R4B is using MILP reformulation now to solve the solution of the bilinear problem: \n\n');
    %%  MILP 
        obj_opt=obj;
        Best_solution=obj_opt;
        Lower_bound=val;
        MILP_time    =0;    %time for solving the MILP
        n_z=size(b_y,1);          
        current_time=time_LDR+MILP_time;
        fprintf('Time=%f      lower bound=%f      best value=%f\n',current_time,Lower_bound,Best_solution)
        Best_solution=obj_opt;
        Lower_bound=val;
        fprintf('Time=%f      lower bound=%f      best value=%f\n',current_time,Lower_bound,Best_solution)
        %finding the value for M in the reformulation
        fprintf('The value of BigM is being calculated:\n')
        [M,time_M] = M_value(A_y,b_y,Q,x_opt,Lower_bound,Best_solution,time_limit-time_LDR-MILP_time,tolerr);
        %% 
        x=sdpvar(n_x,1);
        y=sdpvar(n_y,1);
        z=sdpvar(m_y,1);
        u=sdpvar(n_y,1);
        v=sdpvar(m_y,1);
        assign(x,sol_x);
        assign(y,sol_y);
        assign(z,z_opt);
        assign(u,sign(round(Q'*x_opt-A_y'*z_opt,3)));
        assign(v,abs(sign(round(b_y-A_y*y_opt,3))));
        M=10*M;              
        constraints=[ binary(u), binary(v), z>= 0-tolerr, y>=0-tolerr, x>=0-tolerr, A_x*x >=b_x-tolerr, A_y*y>=b_y-tolerr, A_y'*z<=Q'*x+tolerr, Q'*x-A_y'*z<=M*u+tolerr, y<=M*(1-u), b_y-A_y*y>=-M*v-tolerr, z<= M*(1-v),b_y'*z>=Lower_bound];
        sol_MILP=optimize(constraints,b_y'*z,sdpsettings('usex0',1,'solver','gurobi','verbose', 0,'gurobi.Cutoff',Best_solution-tolerr,'gurobi.TimeLimit',time_limit-time_LDR,'gurobi.MIPGap',tolerr,'savesolveroutput',1));
        obj=sol_MILP.solveroutput.result.objval;
        val=sol_MILP.solveroutput.result.objbound;
        sol_x=value(x);
        sol_y=value(y);
        MILP_time=time_M+sol_MILP.solveroutput.result.runtime;


    end
    

    %reporting the best solution

    time=time_LDR+MILP_time;

    
    
end