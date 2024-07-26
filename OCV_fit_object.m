function object=OCV_fit_object(parameter,exp_data,options)
% Relaxation voltage fitting objective function
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
% object:Value of the loss function

arguments
    parameter (1,:) double
    exp_data (:,2) double
    options.methode string ='Relaxation'
    options.fitness string ='RMSE'
end
methode=options.methode;
fitness=options.fitness;
t_seri=exp_data(:,1);
U_exp=exp_data(:,2);

switch methode
    case 'Relaxation'
        %Ref:
        %Fang, Qiaohua, et al. "A state of health estimation method for 
        % lithium-ion batteries based on voltage relaxation model." 
        % Energies 12.7 (2019): 1349.
        %Qian, Yimin, et al. "Fast open circuit voltage estimation of lithium-ion 
        % batteries using a relaxation model and genetic algorithm." 
        % IEEE Access 10 (2022): 96643-96651.
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
        %Ref:
        %耿陈, et al. "基于弛豫过程特征提取的锂离子电池健康状态估计." 储能科学与技术 12.11 (2023): 3479.
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
        %Ref:
        %Yang, Jie, et al. "Rapid prediction of the open-circuit-voltage of lithium ion batteries 
        % based on an effective voltage relaxation model." Energies 11.12 (2018): 3444.
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

switch fitness
    case 'Chi-2'
        U_ave=mean(U_exp);
        U_diff=U_exp-U_fit;
        object=-(1-sqrt(sum(U_diff.^2))./sqrt(sum((U_exp-U_ave).^2)))*100;
    case 'RMSE'
        diff=U_fit-U_exp;
        RMSE= sqrt(mean(diff.^2));
        object=RMSE;
    otherwise
        error('%s is not surpport',fitness)
end

end