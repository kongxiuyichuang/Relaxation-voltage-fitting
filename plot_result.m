function [RMSE,Chi_2,U_diff]=plot_result(parameter,exp_data,options)
% input:
% parameter:The model parameters to be fitted,
% parameters = [alpha, beta, OCV],when model type is 'Relaxation'
% parameters = [alpha, beta, gamma, OCV],when model type is 'Relaxation-nolinear';
% parameters = [a1,a2,a3,a4,OCV],when model type is 'Nernst';
% exp_data: Time&voltage data of the relaxation process,size:(N*2)
% options is a struct with fields:
% methode: types of fitting models ,the default is 'Relaxation',also provide
% 'Relaxation-nolinear' and 'Nernst'
% fitness: types of the target loss function, support 'RMSE' and
% 'Chi-2',the default is 'RMSE'.
%------------------------------
% output:
% RMSE:RMSE of the experimental data and the fitting data,size:(1*1)
% Chi_2:-(1-RSS/TSS),size:(1*1)
% U_diff:The difference between the experimental data and the fitting
% data,size:(N*1)

arguments
    parameter (1,:) double
    exp_data (:,2) double
    options.methode string ='Relaxation'
    options.fitness string ='RMSE'
end
methode=options.methode;
t_seri=exp_data(:,1);
U_exp=exp_data(:,2);
switch methode
    case 'Relaxation'
        alph=parameter(1);
        beta=parameter(2);
        OCV=parameter(3);
        Ut=zeros(length(t_seri),1);
        Ut(1)=U_exp(1);
        for i=2:length(t_seri)
            Ut(i)=(Ut(i-1)-OCV).*exp(-(t_seri(i)-t_seri(i-1))./(alph*t_seri(i)+beta))+OCV;
            U_fit=Ut;
        end
    case 'Relaxation-nolinear'
        alph=parameter(1);
        beta=parameter(2);
        gamma=parameter(3);
        OCV=parameter(4);
        Ut=zeros(length(t_seri),1);
        Ut(1)=U_exp(1);
        for i=2:length(t_seri)
            Ut(i)=(Ut(i-1)-OCV).*exp(-(t_seri(i)-t_seri(i-1))./(alph*t_seri(i).^gamma+beta))+OCV;
            U_fit=Ut;
        end
    case 'Nernst'
        a1=parameter(1);
        a2=parameter(2);
        a3=parameter(3);
        a4=parameter(4);
        OCV=parameter(5);
        U_fit=zeros(length(t_seri),1);
        U_fit(1)=U_exp(1);
        U_fit(2:end)=OCV-a3.*power(t_seri(2:end),a4).*log(t_seri(2:end))-a1.*power(t_seri(2:end),a2);
    otherwise
        error('%s is not surpport',methode)
end
U_ave=mean(U_exp);
U_diff=U_exp-U_fit;
Chi_2=(1-sqrt(sum(U_diff.^2))./sqrt((sum((U_exp-U_ave).^2))))*100;
RMSE= sqrt(mean(U_diff.^2));
figure
plot(t_seri,U_exp,'-r',LineWidth=2,DisplayName='exp-data')
hold on;
plot(t_seri,U_fit,'-b',LineWidth=2,DisplayName='fit-data')
xlabel('time(s)');ylabel('voltage(mV)');grid on;legend;

end