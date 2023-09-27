function [h_tilda_est,k_v_est,l_tau_est,mse,Kt,PP,gamma_v,gamma_tau,gamma_noise]=GL_VBI(y,N,M,x,hw,h_tilda,k_v,l_tau)
    Max_iter=50;
    p_tilda=24;
    t=1;
    mse   = zeros(Max_iter+1,1);
    Kt    = zeros(Max_iter+1,1);
    hw_old=zeros(N,M);
   %% parameter 
   for it=1:p_tilda
       if it==1
        y2=y'*y;
        gamma_noise  = 100;
        s_0   = floor(p_tilda/2); 
        rho = s_0/p_tilda;   
        tau_h = 1/(y2/(M*N)-1/gamma_noise)/(rho*p_tilda);
        gamma_v=1e6;
        gamma_tau=1e6;
       end
       
   end
   %%
    Mu_tau_post=zeros(p_tilda,1);
    Mu_tau_post_a=zeros(p_tilda,1);
    Mu_v_post=zeros(p_tilda,1);
    Mu_v_post_a=zeros(p_tilda,1);
    Sigma_tau_post_a=zeros(p_tilda,p_tilda);
    Sigma_v_post_a=zeros(p_tilda,p_tilda);
    Sigma_tau_post=zeros(p_tilda,p_tilda);
    Sigma_v_post=zeros(p_tilda,p_tilda);
    Mu_v_pri=[round(k_v(1:end)')+0.8*(rand(1,size(k_v,2))'-0.5);N*rand(1,p_tilda-size(k_v,2))'];
    Mu_tau_pri=[round(l_tau(1:end)')+0.8*(rand(1,size(k_v,2))'-0.5);N*rand(1,p_tilda-size(k_v,2))'];     
%     Mu_v_pri=[round(k_v(1:end)')+0.5;N*rand(1,p_tilda-size(k_v,2))'];
%     Mu_tau_pri=[round(l_tau(1:end)')+0.5;M*rand(1,p_tilda-size(l_tau,2))'];
%     Mu_v_pri=N*rand(1,p_tilda)';
%     Mu_tau_pri=M*rand(1,p_tilda)';
    P_mu=Phi_mu(Mu_v_pri,Mu_tau_pri,p_tilda,N,M,x); 
    P_tau=Phi_tau(Mu_v_pri,Mu_tau_pri,p_tilda,N,M,x);
    P_v=Phi_v(Mu_v_pri,Mu_tau_pri,p_tilda,N,M,x);
    P_hat=Phi_hat(Mu_v_pri,Mu_tau_pri,p_tilda,N,M,(1/gamma_v)*eye(p_tilda),(1/gamma_tau)*eye(p_tilda),x);%?
    %PP=Phi_Phi(Mu_v_pri,Mu_tau_pri,p_tilda,N,M,x);
    PP=P_mu'*P_mu+(P_tau'*P_tau.*(1/gamma_tau*eye(p_tilda)))+(P_v'*P_v.*(1/gamma_v*eye(p_tilda))); 
   %%
    cont=1;
    while cont
    t = t + 1; 
    [ s_0, s, Mu_h_post, Sigma_h_post ] = maxZ( PP, P_hat'*y, M*N, 1/gamma_noise, rho, 1/tau_h );
    if s_0>0
        if s_0<p_tilda
            rho = s_0/p_tilda;
        else
            rho = (p_tilda-1)/p_tilda; 
        end
    else
        rho = 1/p_tilda; 
    end
    PtPt=P_tau'*P_tau;
    Sigma_tau_post_a(s,s)=real(((Mu_h_post(s)*Mu_h_post(s)'+Sigma_h_post(s,s)).').*(PtPt(s,s)));
    Mu_tau_post_a(s)=real(diag(Mu_h_post(s))'*P_tau(:,s)'*(y-P_mu(:,s)*Mu_h_post(s))-diag(P_tau(:,s)'*P_mu(:,s)*Sigma_h_post(s,s)));
    Sigma_tau_post(s,s)=inv(gamma_noise*Sigma_tau_post_a(s,s)+gamma_tau*eye(s_0));
    Mu_tau_post(s)=Sigma_tau_post(s,s)*(gamma_noise*(Mu_tau_post_a(s)+Sigma_tau_post_a(s,s)*Mu_tau_pri(s))+gamma_tau*Mu_tau_pri(s));
    
    PvPv=P_v'*P_v;
    Sigma_v_post_a(s,s)=real(((Mu_h_post(s)*Mu_h_post(s)'+Sigma_h_post(s,s)).').*(PvPv(s,s)));
    Mu_v_post_a(s)=real(diag(Mu_h_post(s))'*P_v(:,s)'*(y-P_mu(:,s)*Mu_h_post(s))-diag(P_v(:,s)'*P_mu(:,s)*Sigma_h_post(s,s)));
    Sigma_v_post(s,s)=inv(gamma_noise*Sigma_v_post_a(s,s)+gamma_v*eye(s_0));
    Mu_v_post(s)=Sigma_v_post(s,s)*(gamma_noise*(Mu_v_post_a(s)+Sigma_v_post_a(s,s)*Mu_v_pri(s))+gamma_v*Mu_v_pri(s));
    

    gamma_v=s_0/(real(trace(Sigma_v_post(s,s))+(Mu_v_post(s)-Mu_v_pri(s))'*(Mu_v_post(s)-Mu_v_pri(s))));
    gamma_tau=s_0/(real(trace(Sigma_tau_post(s,s))+(Mu_tau_post(s)-Mu_tau_pri(s))'*(Mu_tau_post(s)-Mu_tau_pri(s))));
    Py=P_hat'*y;
    gamma_noise  =M*N/real( y2 - 2*real( Mu_h_post(s)'*Py(s)) + Mu_h_post(s)'*PP(s,s)*Mu_h_post(s) + trace(PP(s,s)*Sigma_h_post(s,s)) );
    tau_h = s_0/real( Mu_h_post(s)'*Mu_h_post(s)+trace(Sigma_h_post(s,s)));

    
    
    %
    P_mu=Phi_mu(Mu_v_post,Mu_tau_post,p_tilda,N,M,x);
    P_tau=Phi_tau(Mu_v_post,Mu_tau_post,p_tilda,N,M,x);
    P_v=Phi_v(Mu_v_post,Mu_tau_post,p_tilda,N,M,x);
    P_hat=Phi_hat(Mu_v_post,Mu_tau_post,p_tilda,N,M,Sigma_v_post,Sigma_tau_post,x);
    PP=P_mu'*P_mu+(P_tau'*P_tau.*conj(Sigma_tau_post))+(P_v'*P_v.*conj(Sigma_v_post)); 
    
    
    Mu_v_pri=Mu_v_post;
    Mu_tau_pri=Mu_tau_post;
    
    
    
    
    hw_e=zeros(N,M);
    ind = 1:p_tilda; 
    ind = ind(s);
    k=(1:1:N)';
    l=(1:1:M);
    for i=1:1:length(ind)
        w(i,:,:)=w_Sample((repmat(k,[1,M])-Mu_v_post(ind(i))),(repmat(l,[N,1])-Mu_tau_post(ind(i))),N,M);
        hw_e=hw_e+Mu_h_post(ind(i))*reshape(w(i,:,:),N,M); 
    end
    
    mse(t) = norm(hw_e-hw,"fro")^2/norm(hw,"fro")^2;
    if mse(t)>0.5
        mse(t)=0;
    end
    Kt(t)  = s_0;
    if (norm(hw_e-hw_old,"fro")/norm(hw_old,"fro")<1e-3) || (norm(hw_old)==0&&norm(hw_old-hw_e)==0) || (t >= Max_iter)
        cont = 0;
        mse(t+1:end) = mse(t);
        Kt(t+1:end)  = Kt(t);
    end
    hw_old = hw_e;
    end
h_tilda_est=Mu_h_post(s);
k_v_est=Mu_v_post(s);
l_tau_est=Mu_tau_post(s);
end

function P_mu=Phi_mu(k_v,l_tau,p_tilda,N,M,x)
k=(1:1:N)';
l=(1:1:M);
for i=1:1:p_tilda
w(i,:,:)=w_Sample((repmat(k,[1,M])-k_v(i)),(repmat(l,[N,1])-l_tau(i)),N,M);
end
for i=1:1:p_tilda
Phi(:,:,i)=cconv2(reshape(w(i,:,:),N,M),x)  ;
end
P_mu=reshape(Phi,[],p_tilda);
end
function P_tau=Phi_tau(k_v,l_tau,p_tilda,N,M,x)
k=(1:1:N)';
l=(1:1:M);
for i=1:1:p_tilda
w(i,:,:)=w_Sample_tau((repmat(k,[1,M])-(k_v(i))),(repmat(l,[N,1])-l_tau(i)),l_tau(i),N,M);
end
for i=1:1:p_tilda
Phi(:,:,i)=cconv2(reshape(w(i,:,:),N,M),x)  ;
end
P_tau=reshape(Phi,[],p_tilda);
end
function P_v=Phi_v(k_v,l_tau,p_tilda,N,M,x)
k=(1:1:N)';
l=(1:1:M);
for i=1:1:p_tilda
w(i,:,:)=w_Sample_v((repmat(k,[1,M])-k_v(i)),k_v(i),(repmat(l,[N,1])-l_tau(i)),N,M);
end
for i=1:1:p_tilda
Phi(:,:,i)=cconv2(reshape(w(i,:,:),N,M),x)  ;
end
P_v=reshape(Phi,[],p_tilda);
end 
function P_hat=Phi_hat(k_v,l_tau,p_tilda,N,M,S2_v,S2_tau,x)
k=(1:1:N)';
l=(1:1:M);
for i=1:1:p_tilda
w(i,:,:)=w_Sample_hat((repmat(k,[1,M])),k_v(i),(repmat(l,[N,1])),(l_tau(i)),N,M,S2_v(i,i),S2_tau(i,i));
end
for i=1:1:p_tilda
Phi(:,:,i)=cconv2(reshape(w(i,:,:),N,M),x)  ;
end
P_hat=reshape(Phi,[],p_tilda);
end 
function w=w_Sample_hat(k,m_v,l,m_tau,N,M,s2v,s2tau)
w=1/(M*N)*((1-exp(N*(0.5*s2v-1j*(2*pi/N.*(k-m_v)))))./(1-exp(0.5*s2v-1j*(2*pi/N.*(k-m_v)))) ).*((1-exp(M*(0.5*s2tau-1j*(2*pi/M.*(l-m_tau)))))./(1-exp(0.5*s2tau-1j*(2*pi/M.*(l-m_tau)))));
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

% function PP=Phi_Phi(Mu_v_pri,Mu_tau_pri,p_tilda,N,M,x,Sigma_tau_post,Sigma_v_post)
% for i=1:1:p_tilda
%     for j=1:1:p_tilda
%         
%     end
% end
% end



function [ K, s, w, C ] = maxZ( J, h, M, nu, rho, tau )
%maxZ maximizes the function Z of the binary vector s, see Appendix A of
%the paper

L = size(h,1);
cnst = log(rho/(1-rho)/tau);

K = 0; % number of components
s = false(L,1); % Initialize s
w = zeros(L,1);
C = zeros(L);
u = zeros(L,1);
v = zeros(L,1);
Delta = zeros(L,1);
if L > 1
    cont = 1;
    while cont
        if K<M-1
            v(~s) = nu ./ ( M + nu/tau - real(sum(J(s,~s).*conj(C(s,s)*J(s,~s)),1))/nu );
            u(~s) = v(~s) .* ( h(~s) - J(s,~s)'*w(s))/nu;
            Delta(~s) = log(v(~s)) + u(~s).*conj(u(~s))./v(~s) + cnst;
        else
            Delta(~s) = -1; % dummy negative assignment to avoid any activation
        end
        if ~isempty(h(s))
            Delta(s) = -log(diag(C(s,s))) - w(s).*conj(w(s))./diag(C(s,s)) - cnst;
        end
        [~, k] = max(Delta);
        if Delta(k)>0
            if s(k)==0 % activate
                w(k) = u(k);
                ctemp = C(s,s)*J(s,k)/nu;
                w(s) = w(s) - ctemp*u(k);
                C(s,s) = C(s,s) + v(k)*(ctemp*ctemp');
                C(s,k) = -v(k)*ctemp;
                C(k,s) = C(s,k)';
                C(k,k) = v(k);
                s(k) = ~s(k); 
                K = K+1;
            else % deactivate
                s(k) = ~s(k); 
                K = K-1;
                w(s) = w(s) - C(s,k)*w(k)/C(k,k);
                C(s,s) = C(s,s) - C(s,k)*C(k,s)/C(k,k);
            end
            C = (C+C')/2; % ensure the diagonal is real
        else
            break
        end
    end
elseif L == 1
    if s == 0
        v = nu ./ ( M + nu/tau );
        u = v * h/nu;
        Delta = log(v) + u*conj(u)/v + cnst;
        if Delta>0
            w = u; 
            C = v; 
            s = 1; 
            K = 1;
        end
    else
        Delta = -log(C) - w*conj(w)/C - cnst;
        if Delta>0
            w = 0; 
            C = 0; 
            s = 0; 
            K = 0;
        end
    end
end
end