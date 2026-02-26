clear; clc; close all;
s = settings;
s.matlab.appearance.figure.GraphicsTheme.TemporaryValue = "light";


% Cargar archivo

filename = 'datosRT_PbNaTe_Tarea2MQ.dat';
data = readmatrix(filename);

figure; hold on;
colors = hsv(size(data,2)/2);

for k = 1:(size(data,2)/2)
    
    H   = data(:,2*k-1);
    rho = data(:,2*k);
    
    valid = ~isnan(H) & ~isnan(rho);
    H = H(valid);
    rho = rho(valid);
    
    plot(1./H, rho, 'LineWidth', 1.2, 'Color', colors(k,:));
    
end

xlabel('$(\mu_0 H)^{-1}\,(T^{-1})$', Interpreter='latex');
xlim([0.06 0.13]);
legend(["2K","4K","6K","8K","10K","12K","14K","16K","18K","20K"]);
ylabel('$(\rho/\rho_0)-1$', Interpreter='latex');
grid on;
hold off;
%print(gcf, 'RawDataInvH', '-dpng', '-r300');


% Parámetros

minH = 1/0.13;      % Campo mínimo incluido en FFT
nTemps = size(data,2)/2;

Tlist = 2:2:20;

peakFrequency = zeros(nTemps,1);
peakAmplitude = zeros(nTemps,1);



% Cálculo FFT 

figure;
hold on;

for k = 1:nTemps
    
    colH   = 2*k - 1;
    colRho = 2*k;
    
    H   = data(:,colH);
    rho = data(:,colRho);
    
    valid = ~isnan(H) & ~isnan(rho);
    H   = H(valid);
    rho = rho(valid);
    
    arrsubs = [H rho];
    
    arrfourier = fourierdatav2(arrsubs, minH);
    
    f  = arrfourier(:,1);
    A  = arrfourier(:,2);
    
    % Buscar pico principal (ignorando frecuencia cero)

    [amp, idx] = max(A(2:end));
    idx = idx + 1;
    
    peakFrequency(k) = f(idx);
    peakAmplitude(k) = amp;
    
    plot(f, A, 'LineWidth', 1.2, 'Color', colors(k,:));
    
end

xlabel('Frequency (T)');
ylabel('FFT amplitude');
legend(["2K","4K","6K","8K","10K","12K","14K","16K","18K","20K"]);
xlim([0 500]);
grid on;
hold off;
%print(gcf, 'FFT', '-dpng', '-r300');



% Parte a): Frecuencia del componente oscilatorio

disp('Temperatura (K)   Frecuencia (T)   Amplitud');
disp([Tlist' peakFrequency peakAmplitude]);

Fvalid = peakFrequency;
Fmean  = mean(Fvalid);
Fstd   = std(Fvalid);              
Ferr   = Fstd/sqrt(length(Fvalid));

fprintf('F = %.1f ± %.1f T\n', Fmean, Ferr);


% Parte b): masa efectiva de ciclotrón


Tfit = Tlist(:);
Afit = peakAmplitude(:);

Heff = 11; % Campo promedio en la ventana usada en FFT

model = @(params,T) ...
    params(1) .* ...
    ( (14.7*params(2).*T./Heff) ./ ...
    sinh(14.7*params(2).*T./Heff) );

params0 = [max(Afit) 0.05];

[params,resnorm,~,~,~,~,J] = ...
    lsqcurvefit(model, params0, Tfit, Afit);

N = length(Afit);
p = length(params);

covariance = full((resnorm/(N-p)) * inv(J'*J));
errors = sqrt(diag(covariance));

m_ratio = params(2);
m_error = errors(2);

% Masa calculada

fprintf('m*/m0 = %.3f ± %.3f\n', m_ratio, m_error);

% Figura 

A_norm = peakAmplitude ./ peakAmplitude(1);

xdata = 14.7 * m_ratio * Tlist / (2*pi^2*Heff);

xfit = linspace(min(xdata)*0.5, max(xdata)*1.1, 200);
R_T = @(x) (2*pi^2*x)./sinh(2*pi^2*x);
yfit = R_T(xfit);
yfit = yfit./yfit(1);


figure;
plot(xfit,yfit, '--k'); hold on;
scatter(xdata,A_norm,80,'filled', 'diamond', 'o');

xlabel('$k_B T / \hbar \omega_c$', Interpreter='latex');
ylabel('Amplitude (arb. unit)');
legend("", "$m^* = (0.220 \pm 0.003)\,m_e$", Interpreter='latex');
grid on;
%print(gcf, 'MEff', '-dpng', '-r300');



% Parte c): Temperatura Dingle y tiempo libre medio


invH_max = 0.1;
nTemps = length(Tlist);

ThetaD_list = zeros(nTemps,1);
ThetaD_err  = zeros(nTemps,1);
tau_list    = zeros(nTemps,1);
tau_err    = zeros(nTemps,1);

hbar = 1.054e-34;
kB   = 1.381e-23;

for k = 1:nTemps
    
    colH   = 2*k - 1;
    colRho = 2*k;
    
    H   = data(:,colH);
    rho = data(:,colRho);
    
    valid = ~isnan(H) & ~isnan(rho);
    H = H(valid);
    rho = rho(valid);
    
    invH = 1./H;
    [invH_sorted, idx_sort] = sort(invH);
    rho_sorted = rho(idx_sort);
    
    mask = invH_sorted <= invH_max;
    invH_sorted = invH_sorted(mask);
    rho_sorted  = rho_sorted(mask);
    
    % Máximos de oscilación

    [peaks,locs] = findpeaks(rho_sorted,invH_sorted);
    
    Hi = 1./locs;
    T  = Tlist(k);
    
    modelD = @(params,H) ...
        params(1) .* ...
        ( (14.7*m_ratio*T./H) ./ ...
        sinh(14.7*m_ratio*T./H) ) .* ...
        exp( -14.7*m_ratio*params(2)./H );
    
    params0 = [max(abs(peaks)) 5];
    
    [params,resnorm,~,~,~,~,J] = ...
        lsqcurvefit(modelD, params0, Hi, abs(peaks));
    
    N = length(peaks);
    p = length(params);
    
    covariance = full((resnorm/(N-p)) * inv(J'*J));
    errors = sqrt(diag(covariance));
    
    ThetaD_list(k) = real(params(2));
    ThetaD_err(k)  = real(errors(2));
    
    tau_list(k) = hbar/(2*pi*kB*ThetaD_list(k))*1e12;
    tau_err(k) = tau_list(k)*ThetaD_err(k)/ThetaD_list(k);
    
end


% Cálculos 

disp('T (K)   Theta_D (K)   Error   tau (ps)   Error tau (ps)');
disp([Tlist' ThetaD_list ThetaD_err tau_list tau_err]);



% Parte d): Movilidad electrónica

e  = 1.602e-19;
m0 = 9.11e-31;

m_star = m_ratio * m0;

mu_list = e .* (tau_list/1e12) ./ m_star;   % en m^2/(V s)

% convertir a cm^2/(V s)
mu_list_cm = mu_list * 1e4;
mu_err = 1e4*e.*( (1e-12*tau_err./m_star).^2 + (1e-12*0.003*m0*tau_list./m_star^2).^2 ).^0.5;

% Cálculos

disp('T (K)   mu (cm^2/Vs) Error (cm^2/Vs)');
disp([Tlist' mu_list_cm mu_err]);

