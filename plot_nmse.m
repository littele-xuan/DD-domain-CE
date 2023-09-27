clc
clear
snr=0:5:30;
m=40;
n=32;
nmse=load('E:\OTFS\OTFS\reslut.mat');
nmse_i=load('E:\OTFS\OTFS\reslut0imp.mat');
nmse_l=load('E:\OTFS\OTFS\reslut0lmmse.mat');
nmse_s=load('E:\OTFS\OTFS\reslut0sbl.mat');
nmse=nmse.NMSE;
nmse_i=nmse_i.NMSE;
nmse_l=nmse_l.NMSE;
nmse_s=nmse_s.NMSE;

mse=sum(nmse,2)/size(nmse,2);
mse1=sum(nmse_i,2)/size(nmse_i,2);
mse2=sum(nmse_l,2)/size(nmse_l,2);
mse3=sum(nmse_s,2)/size(nmse_s,2);
% mse_db=10*log10(mse);   
% mse1_db=10*log10(mse1);
% mse2_db=10*log10(mse2);
semilogy(snr,mse,'-rs');
hold on
semilogy(snr,mse1,'-bs');
semilogy(snr,mse2,'-gs');
semilogy(snr,mse3,'-k*');
%305.1559*
%semilogy(snr,10.^(-snr/10),'-+');
grid on
%plot(snr,10*log10(10.^(-snr/10)/(32*40)),'r')
legend('GL-VBI','LS','LMMSE','OG-SBI');
set(gca,'FontName','Times New Roman','FontSize',13);
xlabel('\fontname{Times New Roman}SNR\fontname{Times New Roman}(dB)','FontSize',18);ylabel('\fontname{Times New Roman}NMSE\fontname{Times New Roman}','FontSize',18);

% x=sqrt(0.5)*(randn(m,n)+1j*randn(m,n));
% x2=fft2(x);
% sum(abs(x2).^2)/sum(abs(x).^2)
% sum(sum(abs(x).^2))/(m*n)



