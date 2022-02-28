%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 function [f1,t1]=bound_14(Q0,c0,Q,c,d,A,b)
%    just polytope
%%%%%% min    x'Q0x+2c0'x
%%%%%% s.t.   Ax=b,   A_{p,n}
%%%%%%        x'Q_ix+2c_i'x+d_i<=0 .  i=1 ... m
%%%%%%        0<x<0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n=size(Q0,1);
m=size(Q,3);
p=size(A,1);
m_c=2;
clear('yalmip');
%%%%%%%%%%%%%%%%%%%%%%%%
x=sdpvar(n,1);
lam_ql=sdpvar(m,n,'full');   %  linear part lagrange multiplier quadratic constarint
lam_qc=sdpvar(m,1);          %  constant part lagrange multiplier quadratic constarint
lam_eq=sdpvar(n,n,p);        %  quadratic part lagrange multiplier linear equality constarint
lam_el=sdpvar(p,n,'full');   %  linear part lagrange multiplier linear equality constarint
lam_ec=sdpvar(p,1);          %  constant part lagrange multiplier linear equality constarint
lam_ul=sdpvar(n,n,'full');   %  linear part lagrange multiplier upper bound
lam_uc=sdpvar(n,1);          %  constant part lagrange multiplier upper bound
lam_uq=sdpvar(n,n,n);        %  quadratic part lagrange multiplier linear equality constarint
lam_ll=sdpvar(n,n,'full');   %  linear part lagrange multiplier lower bound
lam_lc=sdpvar(n,1);          %  constant part lagrange multiplier lower bound
lam_lq=sdpvar(n,n,n);        %  quadratic part lagrange multiplier linear equality constarint
opt=sdpvar(1);               %  optimal value
mu_e_l=sdpvar(n,p);
mu_e_u=sdpvar(n,p);
mu_u_l=sdpvar(n,n,'full');
mu_u_u=sdpvar(n,n,'full');
mu_l_l=sdpvar(n,n,'full');
mu_l_u=sdpvar(n,n,'full');
S=sdpvar(n+1);
mu_q_l=sdpvar(n,m,'full');
mu_q_u=sdpvar(n,m,'full');
Sl=sdpvar(n+1,n+1,n);
Su=sdpvar(n+1,n+1,n);
mu_c_q=sdpvar(m,1);
mu_e_q=sdpvar(m,p);
mu_u_q=sdpvar(m,n,'full');
mu_l_q=sdpvar(m,n,'full');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
con=[];
P=x'*Q0*x+2*c0'*x-opt;

B = reshape(Q,n,n*m);
bB = reshape(x'*B,n,m);
anew = x.'*bB;
anew = anew(:);
P=P+(lam_ql*x+lam_qc)'*(anew+2*c*x+d);
con=[con, coefficients(lam_ql*x+lam_qc-mu_e_q*(A*x-b)-...
          mu_u_q*(ones(n,1)-x)-mu_l_q*x-mu_c_q,x)==0];
     
B = reshape(lam_eq,n,n*p);
bB = reshape(x'*B,n,p);
anew = x.'*bB;
anew = anew(:);
P=P+(anew+lam_el*x+lam_ec)'*(A*x-b);


B = reshape(lam_uq,n,n*n);
bB = reshape(x'*B,n,n);
anew = x.'*bB;
anew2 = anew(:);
P=P+(anew2+lam_ul*x+lam_uc)'*(x-ones(n,1));

B = reshape(lam_lq,n,n*n);
bB = reshape(x'*B,n,n);
anew = x.'*bB;
anew1 = anew(:);
P=P+(anew1+lam_ll*x+lam_lc)'*(-x);

    B = reshape(Q,n,n*m);
    bB = reshape(x'*B,n,m);
    anew = x.'*bB;
    anew = anew(:);
    pl=mu_q_l*(anew+2*c*x+d);
    pu=mu_q_u*(anew+2*c*x+d);
   
    B = reshape(Sl,n+1,(n+1)*n);
    bB = reshape([x ;1]'*B,(n+1),n);
    anew = [x ;1].'*bB;
    anew11 = anew(:);
   
    con=[con, coefficients(anew1+lam_ll*x+lam_lc-mu_e_l*(A*x-b)-...
                           mu_u_l*(ones(n,1)-x)-mu_l_l*x+pl-anew11,x)==0];
     
     
     
    B = reshape(Su,n+1,(n+1)*n);
    bB = reshape([x ;1]'*B,(n+1),n);
    anew = [x ;1].'*bB;
    anew22 = anew(:);
   
    con=[con, coefficients( anew2+lam_ul*x+lam_uc-mu_e_u*(A*x-b)-...
                           mu_u_u*(ones(n,1)-x)-mu_l_u*x+pu-anew22,x)==0];  
for i=1:n
  con=[con,Su(:,:,i)>=0];
  con=[con,Sl(:,:,i)>=0];
end
con=[con,  mu_l_q>=0, mu_u_q>=0, mu_c_q>=0, mu_u_l>=0, mu_u_u>=0, mu_l_l>=0,...
                     mu_l_u>=0, S>=0, mu_q_l>=0, mu_q_u>=0];
P=P-[x' 1]*S*[x;1];
con=[con,coefficients(P,x)>=0];
ops = sdpsettings('solver','mosek','verbose',0);
g0=optimize(con,-opt,ops);
f1=value(opt);
t1=g0.solvertime;
end

