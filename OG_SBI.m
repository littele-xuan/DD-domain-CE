function [h_tilda_est,k_v_est,l_tau_est,mse]=OG_SBI(y,N,M,x,hw,h_tilda,k_v,l_tau)
k=(1:1:N)';
l=(1:1:M);
kk=repmat([k(1:N-1);0],[M,1]);
ll=reshape(repmat([0,l(1:M-1)],[N,1]),[],1);
iter_num=50;
e=0.4;
kk=kk+e;
ll=ll+e;
hw_old=zeros(N,M);
for i=1:1:M*N
w(i,:,:)=w_Sample((repmat(k,[1,M])-(kk(i))),(repmat(l,[N,1])-(ll(i))),N,M);
Phi(:,:,i)=cconv2(reshape(w(i,:,:),N,M),x)  ;
end
phi0=reshape(Phi,[],M*N);
kv_est=zeros(M*N,1);
ltau_est=zeros(M*N,1);
phi_tau=phi_t(kk,ll,N,M,x);
phi_nu=phi_v(kk,ll,N,M,x);
phi=phi0+phi_tau*diag(ltau_est)+phi_nu*diag(kv_est);

alpha=abs(phi'*y);
b0=M*N/norm(y,"fro")^2;

rho=M*N/sum(alpha);
c=1;
d=1/b0;

for t=1:1:iter_num
%Sigma_h=inv(b0*(phi'*phi)+diag(1./alpha));
Sigma_h=diag(alpha)-diag(alpha)*phi'*inv(1/b0*eye(M*N)+phi*diag(alpha)*phi')*phi*diag(alpha);
mu_h=Sigma_h*(b0*phi'*y);
b0=(c-1+M*N)/(d+norm(y-phi*mu_h,"fro")^2+real(trace(phi'*phi*Sigma_h)));
%b0=(c-1+M*N)/(d+norm(y-phi*mu_h,"fro")^2+1/b0*sum(1-real(diag(Sigma_h))./alpha));
d=1/b0;
alpha=(sqrt(1+4*rho*(abs(mu_h).^2+diag(real(Sigma_h))))-1)/(2*rho);
rho=M*N/sum(alpha);

A_v=real(conj(phi_nu'*phi_nu).*(mu_h*mu_h'+Sigma_h.'));
%b_v=real(diag(mu_h.')*phi_nu.'*conj(y)-diag((mu_h*mu_h'+Sigma_h)*(phi0+phi_tau*diag(ltau_est))'*phi_nu));
b_v=real(conj(mu_h).*(phi_nu'*(y-(phi0+phi_tau*diag(ltau_est))*mu_h))-diag(phi_nu'*(phi0+phi_tau*diag(ltau_est))* Sigma_h));
A_tau=real(conj(phi_tau'*phi_tau).*(mu_h*mu_h'+Sigma_h.')) ;
%b_tau=real(diag(mu_h.')*phi_tau.'*conj(y)-diag((mu_h*mu_h'+Sigma_h)*(phi0+phi_nu*diag(kv_est))'*phi_tau));
b_tau=real(conj(mu_h).*(phi_tau'*(y-(phi0+phi_nu*diag(kv_est))*mu_h))-diag(phi_tau'*(phi0+phi_nu*diag(kv_est))* Sigma_h));
kv_est=(A_v)\b_v;
kv_est(kv_est>0.5)=0.5;
kv_est(kv_est<-0.5)=-0.5;
kv_est(diag(A_v)==0)=0;
ltau_est=(A_tau)\b_tau;
ltau_est(ltau_est>0.5)=0.5;
ltau_est(ltau_est<-0.5)=-0.5;
ltau_est(diag(A_tau)==0)=0;

kk=kk+kv_est;
ll=ll+ltau_est;
for i=1:1:M*N
w(i,:,:)=w_Sample((repmat(k,[1,M])-(kk(i))),(repmat(l,[N,1])-(ll(i))),N,M);
Phi(:,:,i)=cconv2(reshape(w(i,:,:),N,M),x)  ;
end
phi0=reshape(Phi,[],M*N);
kv_est=zeros(M*N,1);
ltau_est=zeros(M*N,1);
phi_tau=phi_t(kk,ll,N,M,x);
phi_nu=phi_v(kk,ll,N,M,x);
% norm(y-phi*mu_h,"fro")^2
phi=phi0+phi_tau*diag(ltau_est)+phi_nu*diag(kv_est);



    [h_main,pos]=sort(alpha);
    pos=pos(end-50:end);
    h_main=mu_h(pos);
    Mu_v_post=kk(pos)+kv_est(pos);
    Mu_tau_post=ll(pos)+ltau_est(pos);
    hw_e=zeros(N,M);
    for i=1:1:length(pos)
        w(i,:,:)=w_Sample((repmat(k,[1,M])-Mu_v_post((i))),(repmat(l,[N,1])-Mu_tau_post((i))),N,M);
        hw_e=hw_e+h_main((i))*reshape(w(i,:,:),N,M); 
    end
     mse(t) = norm(hw_e-hw,"fro")^2/norm(hw,"fro")^2
      
     if (norm(hw_e-hw_old,"fro")/norm(hw_old,"fro")<1e-3) || (norm(hw_old)==0&&norm(hw_old-hw_e)==0) || (t >= iter_num)||(t>10&& mse(t)> mse(t-1))
        mse(t-1:end) = mse(t-1);
        break
    end
    hw_old = hw_e;
end
    h_tilda_est=h_main;
    k_v_est=Mu_v_post;
    l_tau_est=Mu_tau_post;


    



end
























function P_tau=phi_t(k_v,l_tau,N,M,x)
k=(1:1:N)';
l=(1:1:M);
for i=1:1:M*N
w(i,:,:)=w_Sample_tau((repmat(k,[1,M])-(k_v(i))),(repmat(l,[N,1])-l_tau(i)),l_tau(i),N,M);
end
for i=1:1:M*N
Phi(:,:,i)=cconv2(reshape(w(i,:,:),N,M),x)  ;
end
P_tau=reshape(Phi,[],M*N);
end
function P_v=phi_v(k_v,l_tau,N,M,x)
k=(1:1:N)';
l=(1:1:M);
for i=1:1:N*M
w(i,:,:)=w_Sample_v((repmat(k,[1,M])-k_v(i)),k_v(i),(repmat(l,[N,1])-l_tau(i)),N,M);
end
for i=1:1:N*M
Phi(:,:,i)=cconv2(reshape(w(i,:,:),N,M),x)  ;
end
P_v=reshape(Phi,[],N*M);
end 


function w=w_Sample_tau(k,l,l_tau,N,M)
w=1/(M*N)*(1j*(M-1)*pi/M*exp(-1j*(M-1)*pi/M.*l).*(sin(pi*l)./sin(pi*l/M))+exp(-1j*(M-1)*pi/M.*l).*((sin(pi*l).*cos(pi*l/M)*pi/M-cos(pi*l).*sin(pi*l/M)*pi)./(sin(pi*l/M).^2)) ).*(exp(-1j*(N-1)*pi/N.*k).*(sin(pi*k)./sin(pi*k/N)));
end
function w=w_Sample_v(k,k_v,l,N,M)
w=1/(M*N)*(1j*(N-1)*pi/N*exp(-1j*(N-1)*pi/N.*k).*(sin(pi*k)./sin(pi*k/N))+exp(-1j*(N-1)*pi/N.*k).*((sin(pi*k).*cos(pi*k/N)*pi/N-cos(pi*k).*sin(pi*k/N)*pi)./(sin(pi*k/N).^2)) ).*(exp(-1j*(M-1)*pi/M.*l).*(sin(pi*l)./sin(pi*l/M)));
end
function w=w_Sample(k,l,N,M)
w=1/(M*N)*(exp(-1j*(N-1)*pi/N.*k).*(sin(pi*k)./sin(pi*k/N))).*(exp(-1j*(M-1)*pi/M.*l).*(sin(pi*l)./sin(pi*l/M)));
end


