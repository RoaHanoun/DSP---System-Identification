% DSP Project
% Leen Abu Omar 1190113
%Hana Kafri 119
%Roa Hanoun 119
% PART 1

close all;
clc
%%
% A
% DONE

N = 2000; % defining the Number Of Samples
n = 0:N-1; 
x = cos(0.03*pi*n); % defining the input signal

% Plotting the input signal
figure(1);
subplot(1,1,1);
plot(n,x);
xlabel('Sample index (n)');
ylabel('Amplitude');
title('Pt.1.A Input signal x[n]= cos(0.03\pi n)');



%%
% B
% DONE

[fir1, ww1] = freqz(h, 1, 'whole'); % returns a complex point frequency response corresponding to the second-order sections matrix 

figure
subplot(2,1,1)
plot(ww1/(2*pi), abs(fir1))
title('Pt.1.B mag1 (abs)')

subplot(2,1,2)
plot(ww1/(2*pi), angle(fir1))
title('phase1')


%% 
% C 
% DONE

N =2000; %NFFT-point DFT      
n = 0:N-1; 

x = cos(0.03*pi*n);

y=fft(x); % compute DFT using FFT        

figure
subplot(2, 1, 1);
plot(n,abs(y));      
title('Pt.1.C Double Sided FFT - without FFTShift');        
xlabel('Sample points (N-point DFT)')        
ylabel('DFT Values');

subplot(2, 1, 2);
plot(n,angle(y));      
title('Double Sided FFT - without FFTShift');        
xlabel('Sample points (N-point DFT)')        
ylabel('DFT Values');


%%
% D

N = 2000;
M = 3;
n = 0:N-1;

x = cos(0.03*pi*n);

h = [1, -2, 4];

mu = 0.01; %step size
w = [0 0 0];
wt = transpose(w);% returns the nonconjugate transpose of w

b = [1 -2 4];
a = 1;
d = filter(b,a,x);% filters the input data x using the numerator and denominator coefficients


for i = M:N
    c = x(i:-1:i-M+1);
    y(i) =   c * wt;
    e(i) = d(i) - y(i); %Computing the error signal
    w = w + 2*mu*e(i)*c; % Updating the weight vector
    wt = transpose(w); % take transpose of w
end


% part 2: j
j = e.^2; %calculating J

% part 3
lo = 10 * log10(j);

figure
subplot(3,1,1)
plot(abs(e)) ;
title('Pt.1.D Error curve') ;
ylim([0 2.5])

subplot(3,1,2)
plot(j) ;
title('j = e^2') ;

subplot(3,1,3)
plot(lo) ;
title('10 log10 (j)') ;

% E

[fir1, ww1] = freqz(h, 1, 'whole');% returns a complex point frequency response corresponding to the second-order sections matrix 
[fir2, ww2] = freqz(w, 1, 'whole');% returns a complex point frequency response corresponding to the second-order sections matrix 

figure
subplot(2,2,1)
plot(ww1/(2*pi), abs(fir1))
title('Pt.1.E mag1 (abs)')

subplot(2,2,3)
plot(ww1/(2*pi), angle(fir1))
title('phase1')

subplot(2,2,2)
plot(ww2/(2*pi), abs(fir2))
title('mag2 (abs)')

subplot(2,2,4)
plot(ww2/(2*pi), angle(fir2))
title('phase2')



%%
% F
% mu = 0.001

% repeat D

N = 2000;
M = 3;
n = 0:N-1;

x = cos(0.03*pi*n);

h = [1, -2, 4];

mu = 0.001; %step size
w = [0 0 0];
wt = transpose(w);% returns the nonconjugate transpose of w

b = [1 -2 4];
a = 1;
d = filter(b,a,x);

%for loop to apply the LMS Algorithm
for i = M:N
    c = x(i:-1:i-M+1);
    y(i) =   c * wt;
    e(i) = d(i) - y(i); %Computing the error signal
    w = w + 2*mu*e(i)*c; % Updating the weight vector
    wt = transpose(w); % take transpose of w
end


% part 2 : j
j = e.^2;%calculating J

% part 3
lo = 10 * log10(j);

figure
subplot(3,1,1)
plot(abs(e)) ;
title('Pt.1.F Error curve') ;
ylim([0 2.5])

subplot(3,1,2)
plot(j) ;
title('j = e^2') ;

subplot(3,1,3)
plot(lo) ;
title('10 log10 (j)') ;

% remove it 
% repeat E

[fir1, ww1] = freqz(h, 1, 'whole');% returns a complex point frequency response corresponding to the second-order sections matrix 
[fir2, ww2] = freqz(w, 1, 'whole');% returns a complex point frequency response corresponding to the second-order sections matrix 

figure
subplot(2,2,1)
plot(ww1/(2*pi), abs(fir1))
title('Pt.1.F mag1 (abs)')

subplot(2,2,3)
plot(ww1/(2*pi), angle(fir1))
title('phase1')

subplot(2,2,2)
plot(ww2/(2*pi), abs(fir2))
title('mag2 (abs)')

subplot(2,2,4)
plot(ww2/(2*pi), angle(fir2))
title('phase2')



%% 
% G = 40 db
% DONE

% repeat D

N = 2000;
M = 3;
n = 0:N-1;

x1 = cos(0.03*pi*n);

xG =awgn(x1,40,'measured'); % Add 40dB of zero-mean white Gaussian noise

h = [1, -2, 4];
%d = conv(x,h);

% plot it once with mu1, then change mu1 in the loop to mu2
mu1 = 0.01; %step size
mu2 = 0.001;

w = [0 0 0];
wt = transpose(w);% returns the nonconjugate transpose of w

b = [1 -2 4];
a = 1;
d = filter(b,a,xG);%to filter the input data xG using the numerator and denominator coefficients

%for loop to apply the LMS Algorithm using the 40dB of zeros mean white
%Gaussian noise 
for i = M:N
    c = xG(i:-1:i-M+1);
    y(i) =   c * wt;
    e(i) = d(i) - y(i); %Computing the error signal
    w = w + 2*mu1*e(i)*c; % Updating the weight vector
    wt = transpose(w); % take transpose of w
end

% 2 j
j = e.^2;

% 3
lo = 10 * log10(j);

figure
subplot(3,1,1)
plot(abs(e)) ;
title('Pt.1.G Error curve') ;
ylim([0 2.5])

subplot(3,1,2)
plot(j) ;
title('j = e^2') ;

subplot(3,1,3)
plot(lo) ;
title('10 log10 (j)') ;

% repeat E

[fir1, ww1] = freqz(h, 1, 'whole');
[fir2, ww2] = freqz(w, 1, 'whole');

figure
subplot(2,2,1)
plot(ww1/(2*pi), abs(fir1))
title('Pt.1.G mag1 (abs)')

subplot(2,2,3)
plot(ww1/(2*pi), angle(fir1))
title('phase1')

subplot(2,2,2)
plot(ww2/(2*pi), abs(fir2))
title('mag2 (abs)')

subplot(2,2,4)
plot(ww2/(2*pi), angle(fir2))
title('phase2')



%%
% H = 30 db
% DONE

% repeat D

N = 2000;
M = 3;
n = 0:N-1;

x1 = cos(0.03*pi*n);

xH =awgn(x1,30,'measured'); % Add 30dB of zero-mean white Gaussian noise

h = [1, -2, 4]; 

% plot it once with mu1, then change mu1 in the loop to mu2
mu1 = 0.01; %step size
mu2 = 0.001;

w = [0 0 0];
wt = transpose(w); % returns the nonconjugate transpose of w

%numerator and denominator coefficients
b = [1 -2 4];
a = 1;

d = filter(b,a,xH); % filters the input data xH using the numerator and denominator coefficients

%for loop to apply the LMS Algorithm using the 30dB of zeros mean white
%Gaussian noise 
for i = M:N
    c = xH(i:-1:i-M+1);
    y(i) =   c * wt;
    e(i) = d(i) - y(i); %Computing the error signal
    w = w + 2*mu1*e(i)*c; % Updating the weight vector
    wt = transpose(w); % take transpose of w
end

% 2 j
j = e.^2;

% 3
lo = 10 * log10(j);
f = 1:N-M;

figure
subplot(3,1,1)
plot(abs(e)) ;
title('Pt.1.H Error curve') ;
ylim([0 2.5])

subplot(3,1,2)
plot(j) ;
title('j = e^2') ;

subplot(3,1,3)
plot(lo) ;
title('10 log10 (j)') ;

% repeat E

[fir1, ww1] = freqz(h, 1, 'whole');
[fir2, ww2] = freqz(w, 1, 'whole');

figure
subplot(2,2,1)
plot(ww1/(2*pi), abs(fir1))
title('Pt.1.H mag1 (abs)')

subplot(2,2,3)
plot(ww1/(2*pi), angle(fir1))
title('phase1')

subplot(2,2,2)
plot(ww2/(2*pi), abs(fir2))
title('mag2 (abs)')

subplot(2,2,4)
plot(ww2/(2*pi), angle(fir2))
title('phase2')


%%
% I : J avg
% DONE

% repeat D

N = 2000;
M = 3;
n = 0:N-1;

x1 = cos(0.03*pi*n);
xG =awgn(x1,40,'measured'); % Add 40dB of zero-mean white Gaussian noise

h = [1, -2, 4];
%d = conv(x,h);
jSum = 0;
mu1 = 0.01; %step size
mu2 = 0.001;

I = 1000;
% loop for 1000
for v=1:I

    w = [0 0 0];
    wt = transpose(w);% returns the nonconjugate transpose of w

    b = [1 -2 4];
    a = 1;
    d = filter(b,a,xG);% filters the input data xG using the numerator and denominator coefficients
    % loop for 2000, to find e
    for i = M:N
        c = xG(i:-1:i-M+1);
        y(i) =   c * wt;
        e(i) = d(i) - y(i); %Computing the error signal
        w = w + 2*mu2*e(i)*c; % Updating the weight vector
        wt = transpose(w); % take transpose of w
    end

    % 2 j
    j = e.^2;

    jSum = jSum + j;% calculating the Summation of all J values

end

Javg = jSum / I;% calculating the average value of J

% 3
lo = 10 * log10(Javg);

figure
subplot(2,1,1)
plot(Javg) ;
title('Pt.1.I j = e^2') ;

subplot(2,1,2)
plot(lo) ;
title('10 log10 (j)') ;


