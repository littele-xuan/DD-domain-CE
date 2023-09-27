%xiao xuan, 2023
%%
clc
clear
%%
N=32; %Doopler grid i.e.=k
M=40; %Delay grid i.e.=l
SymOrder=2; %bits per symbol
p=6; %path number
snr=30;
sim_num=1;
NMSE=zeros(length(snr),sim_num);
overhead=1;

%%
% nmse_l=load('E:\OTFS\OTFS\reslut0lmmse.mat');
% NMSE=nmse_l.NMSE;
%n2=0;
for ii=1:1:length(snr)
for j=1:1:sim_num
Snr=snr(ii);    
BitSource=randi([0 2^SymOrder-1],N,M);
x=qammod(BitSource,2^SymOrder,'UnitAveragePower',true);
x=[x(:,1:int32(M*overhead)),zeros(N,M-int32(M*overhead))];
h_tilda=1/sqrt(2)*(randn(1,p)+1j*randn(1,p));
k_v=N*rand(1,p);
l_tau=M*rand(1,p);
k=(1:1:N)';
l=(1:1:M);
hw=zeros(N,M);
for i=1:1:p
w(i,:,:)=w_Sample((repmat(k,[1,M])-k_v(i)),(repmat(l,[N,1])-l_tau(i)),N,M);
hw=hw+h_tilda(i)*reshape(w(i,:,:),N,M);
end
hw_norm=(norm(hw,"fro"));
h_tilda=h_tilda/hw_norm;
hw=hw/hw_norm;
y_true=cconv2(hw,x);
y_true_imp=cconv2(hw,x,1);
y_awgn_imp=y_true_imp+ abs(y_true_imp).*10^(-Snr/20).*sqrt(0.5).*(randn(size(y_true_imp))+1i*randn(size(y_true_imp)));
y_awgn=y_true+ abs(y_true).*10^(-Snr/20).*sqrt(0.5).*(randn(size(y_true))+1i*randn(size(y_true)));
%y_awgn=awgn(y_true,Snr,'measured');
%y_awgn_imp=awgn(y_true_imp,Snr,'measured');
y=reshape(y_awgn,[],1);
for i=1:1:p
Phi(:,:,i)=cconv2(reshape(w(i,:,:),N,M),x)  ;
end
phi=reshape(Phi,[],p);
%[h_tilda_est,k_v_est,l_tau_est,nmse,Kt]=GL_VBI(y,N,M,x,hw,h_tilda,k_v,l_tau);
[h_tilda_est_sbl,k_v_est_sbl,l_tau_est_sbl,nmse]=OG_SBI(y,N,M,x,hw,h_tilda,k_v,l_tau);
%[nmse,h_imp]=CE_imp(y_awgn_imp,hw,x);
%[nmse,h_imp]=CE_lmmse(y_awgn_imp,hw,x,Snr);
%n2=n2+norm(y_true,"fro")^2/norm(x,"fro")^2;
NMSE(ii,j)=nmse(end);
str=['the snr=' num2str(Snr) ' the sim_num=' num2str(j) ' NMSE= ' num2str(nmse(end))];
disp(str);


end
end
NMSE(find(isnan(NMSE)==1))=0;
%save('E:\OTFS\OTFS\reslut0sbl.mat','NMSE');
%save('E:\OTFS\OTFS\reslut.mat','NMSE');
%save('E:\OTFS\OTFS\reslut0imp.mat','NMSE');
%save('E:\OTFS\OTFS\reslut0lmmse.mat','NMSE');



function [nmse,h]=CE_lmmse(y,hw,x1,snr)
y1=ifftshift(y/numel(y));
y2=ifft2(y1); %点乘的结果
x_ifft=ifft2(x1);
h_ifft=y2.*conj(x_ifft)./(abs(x_ifft).^2+abs(x_ifft).^2.*10^(-snr/10));
h=fft2(h_ifft);
nmse=norm(hw-h,"fro")^2/norm(hw,"fro")^2;
if nmse>1e6
    nmse=0;
end



end


function [nmse,h]=CE_imp(y,hw,x1)
y1=ifftshift(y/numel(y));
y2=sqrt(numel(y1))*ifft2(y1); %点乘的结果
x_ifft=sqrt(numel(x1))*ifft2(x1);
h_ifft=y2./(x_ifft);
h=fft2(h_ifft);
nmse=norm(hw-h,"fro")^2/norm(hw,"fro")^2;
if nmse>1e6
    nmse=0;
end
end

function w=w_Sample(k,l,N,M)
w=1/(M*N)*(exp(-1j*(N-1)*pi/N.*k).*(sin(pi*k)./sin(pi*k/N))).*(exp(-1j*(M-1)*pi/M.*l).*(sin(pi*l)./sin(pi*l/M)));
end