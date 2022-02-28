%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 function [f1,t1]=bound_1(Q0,c0,Q,c,d,A,b)
%    just polytope
%%%%%% min    x'Q0x+2c0'x
%%%%%% s.t.   Ax=b,   A_{p,n}
%%%%%%        x'Q_ix+2c_i'x+d_i<=0 .  i=1 ... m
%%%%%%        0<x<0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n=size(Q0,1);
m=size(Q,3);
p=size(A,1);
%%%%%%%%%%%%%%%%%%%%%%%%
x=sdpvar(n,1);
lam_q=sdpvar(m,1);   %  lagrange multiplier quadratic constarint
lam_el=sdpvar(p,n);  %  linear part lagrange multiplier linear equality constarint
lam_ec=sdpvar(p,1);  %  constant part lagrange multiplier linear equality constarint
lam_ul=sdpvar(n,n);  %  linear part lagrange multiplier upper bound
lam_uc=sdpvar(n,1);  %  constant part lagrange multiplier upper bound
lam_ll=sdpvar(n,n);  %  linear part lagrange multiplier lower bound
lam_lc=sdpvar(n,1);  %  constant part lagrange multiplier lower bound
opt=sdpvar(1);       %  optimal value
mu_e_l=sdpvar(n,p);
mu_e_u=sdpvar(n,p);
mu_u_l=sdpvar(n,n,'full');
mu_u_u=sdpvar(n,n,'full');
mu_l_l=sdpvar(n,n,'full');
mu_l_u=sdpvar(n,n,'full');
S=sdpvar(n+1);
mu_c_l=sdpvar(n,1);
mu_c_u=sdpvar(n,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
con=[];
P=x'*Q0*x+2*c0'*x-opt;

B = reshape(Q,n,n*m);
bB = reshape(x'*B,n,m);
anew = x.'*bB;
anew = anew(:);
P=P+lam_q'*(anew+2*c*x+d);


P=P+(lam_el*x+lam_ec)'*(A*x-b);
P=P+(lam_ul*x+lam_uc)'*(x-ones(n,1));
P=P+(lam_ll*x+lam_lc)'*(-x);
 con=[con, coefficients(lam_ll*x+lam_lc-mu_e_l*(A*x-b)-...
          mu_u_l*(ones(n,1)-x)-mu_l_l*x-mu_c_l,x)==0];
    con=[con, coefficients(lam_ul*x+lam_uc-mu_e_u*(A*x-b)-...
          mu_u_u*(ones(n,1)-x)-mu_l_u*x-mu_c_u,x)==0];

con=[con,lam_q>=0, mu_u_l>=0,mu_u_u>=0, mu_l_l>=0,mu_l_u>=0, mu_c_l>=0, mu_c_u>=0, S>=0];
P=P-[x' 1]*S*[x;1];
con=[con,coefficients(P,x)==0];
ops = sdpsettings('solver','mosek','verbose',0);
g0=optimize(con,-opt,ops);
f1=value(opt);
t1=g0.solvertime;
end

