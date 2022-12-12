clear;
%Project Tasks
%Set L=25 and number of noise realizations N_it = 1e5
L = 25;
N_it = 1e5;
%Generate vector x of length L randomly
x = exp(1j*2*pi*rand(L,1));
%xH is the Hermitian transpose and xH*x = L
xH = x';
%SNR vector from -5dB to 25dB, in steps of 2dB
SNR_dB = -5:2:25;
%SNR in linear scale
SNR_L = 10.^(SNR_dB/10);
%Variance vector sigma_square, from variance = xH*x/SNR
sigma_square = (xH*x)./SNR_L;
%Initialize array to store P[||y||^2>=beta], i.e. simulated probability, PFA = 0.01
P_sim_1 = zeros(1,length(SNR_dB));
%Initialize array to store theoretical probability of detection P_D, PFA = 0.01
P_D_1 = zeros(1,length(SNR_dB));
%Starting loop over SNR vector, and PFA = 0.01
for k=1:length(SNR_dB) 
    norm_y_square = zeros(1,N_it); %Array to store ||y||^2
    PFA = 0.01;
    beta = sigma_square(k)*gammaincinv(PFA,L,'upper');
    for m=1:N_it
        %Generating zero-mean white gaussian noise vector
        %variance=sigma_square, therefore we need to take square root 
        n = randn(length(x),1)*sqrt(sigma_square(k));
        y = x + n;
        norm_y_square(m)=norm(y)^2;
    end
    P_sim_1(k)=length(find(norm_y_square>=beta))/N_it;
    P_D_1(k)=marcumq(sqrt(2)*norm(x)/sqrt(sigma_square(k)),sqrt(2*beta)/sqrt(sigma_square(k)),L);
end

%Initialize array to store P[||y||^2>=beta], i.e. simulated probability, PFA = 0.001
P_sim_2 = zeros(1,length(SNR_dB));
%Initialize array to store probability of detection, P_D, PFA = 0.001
P_D_2 = zeros(1,length(SNR_dB));
%Starting loop over SNR vector, and PFA = 0.001
for k=1:length(SNR_dB) 
    norm_y_square = zeros(1,N_it); %Array to store ||y||^2
    PFA_2 = 0.001;
    beta = sigma_square(k)*gammaincinv(PFA_2,L,'upper');
    for m=1:N_it
        %Generating zero-mean white gaussian noise vector
        %variance=sigma_square, therefore we need to take square root 
        n = randn(length(x),1)*sqrt(sigma_square(k));
        y = x + n;
        norm_y_square(m)=norm(y)^2;
    end
    P_sim_2(k)=length(find(norm_y_square>=beta))/N_it;
    P_D_2(k)=marcumq(sqrt(2)*norm(x)/sqrt(sigma_square(k)),sqrt(2*beta)/sqrt(sigma_square(k)),L);
end

%Initialize array to store P[||y||^2>=beta], i.e. simulated probability, PFA = 0.1
P_sim_3 = zeros(1,length(SNR_dB));
%Initialize array theoretical Probability of detection array, P_D, PFA = 0.1
%Starting loop over SNR vector, and PFA = 0.1
P_D_3 = zeros(1,length(SNR_dB));
for k=1:length(SNR_dB) 
    norm_y_square = zeros(1,N_it); %Array to store ||y||^2
    PFA_3 = 0.1;
    beta = sigma_square(k)*gammaincinv(PFA_3,L,'upper');
    for m=1:N_it
        %Generating zero-mean white gaussian noise vector
        %variance=sigma_square, therefore we need to take square root 
        n = randn(length(x),1)*sqrt(sigma_square(k));
        y = x + n;
        norm_y_square(m)=norm(y)^2;
    end
    P_sim_3(k)=length(find(norm_y_square>=beta))/N_it;
    P_D_3(k)=marcumq(sqrt(2)*norm(x)/sqrt(sigma_square(k)),sqrt(2*beta)/sqrt(sigma_square(k)),L);
end
figure(1);
plot(SNR_dB, P_sim_1, 'b'); %Plot of Simulated Probability of Detection, PFA = 0.01
%semilogy(SNR_dB, P_sim_1, 'b'); %Plot of Simulated Probability of Detection, PFA = 0.01
hold on;
plot(SNR_dB, P_D_1, 'r');%Plot of Theoretical Probability of Detection, PFA = 0.01
%semilogy(SNR_dB, P_D_1, 'r');%Plot of Theoretical Probability of Detection, PFA = 0.01
hold on;
plot(SNR_dB, P_sim_2, 'g'); %Plot of Simulated Probability of Detection, PFA = 0.001
%semilogy(SNR_dB, P_sim_2, 'g'); %Plot of Simulated Probability of Detection, PFA = 0.001
hold on;
plot(SNR_dB, P_D_2, 'c');%Plot of Theoretical Probability of Detection, PFA = 0.001
%semilogy(SNR_dB, P_D_2, 'c');%Plot of Theoretical Probability of Detection, PFA = 0.001
hold on;
plot(SNR_dB, P_sim_3, 'm'); %Plot of Simulated Probability of Detection, PFA = 0.1
%semilogy(SNR_dB, P_sim_3, 'm'); %Plot of Simulated Probability of Detection, PFA = 0.1
hold on;
plot(SNR_dB, P_D_3, 'y');%Plot of Theoretical Probability of Detection, PFA = 0.1
%semilogy(SNR_dB, P_D_3, 'y');%Plot of Theoretical Probability of Detection, PFA = 0.1
hold on;
title('Simulated and Theoretical Probability of Detection for different PFA scenarios')
xlabel('SNR (dB)') 
ylabel('Probability of Detection, P_D')
grid on;
legend({'Simulation, PFA = 0.01', 'Theory, PFA = 0.01','Simulation, PFA = 0.001', 'Theory, PFA = 0.001','Simulation, PFA = 0.1', 'Theory, PFA = 0.1'},"Location","northeastoutside");
set(gcf,'position',[300,200,800,400])
hold off;

figure(2);
subplot(2,1,1)
plot(SNR_dB, P_sim_1, 'b'); %Plot of Simulated Probability of Detection, PFA = 0.01
hold on;
plot(SNR_dB, P_sim_2, 'g'); %Plot of Simulated Probability of Detection, PFA = 0.001
hold on;
plot(SNR_dB, P_sim_3, 'm'); %Plot of Simulated Probability of Detection, PFA = 0.1
legend('Simulation, PFA = 0.01','Simulation, PFA = 0.001','Simulation, PFA = 0.1');
title('Simulated Probability of Detection for different PFA scenarios')
xlabel('SNR (dB)') 
ylabel('Probability of Detection, P_D')
grid on;

subplot(2,1,2)
plot(SNR_dB, P_D_1, 'r');%Plot of Theoretical Probability of Detection, PFA = 0.01
hold on;
plot(SNR_dB, P_D_2, 'c');%Plot of Theoretical Probability of Detection, PFA = 0.001
hold on;
plot(SNR_dB, P_D_3, 'y');%Plot of Theoretical Probability of Detection, PFA = 0.1
legend('Theory, PFA = 0.01','Theory, PFA = 0.001','Theory, PFA = 0.1');
title('Theoretical Probability of Detection for different PFA scenarios')
xlabel('SNR (dB)') 
ylabel('Probability of Detection, P_D')
grid on;
set(gcf,'position',[300,200,800,400])