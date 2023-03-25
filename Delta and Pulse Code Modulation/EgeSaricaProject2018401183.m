%% TASK 1: PULSE CODE MODULATİON
t_end = 2; % time length
fs = 200; % sampling frequency

% Discrete time vector
t = linspace(0.005, t_end, t_end*fs);

% Create signal m(t) in discrete time points
m_dt = -2*cos(200*pi*t) + sin(50*pi*t);

% Divide range,2mp, into 128 intervals
q_lvl = linspace(-3, 3, 128);


indi = zeros(size(m_dt)); % preallocate indices array
%compare the value of m_dt with quantization levels. When you find
%m_dt<q_lvl, compare m_dt with the upper and lower quantization levels.
%Assign the value of m_dt to the closer quantization level.
for i = 1:length(m_dt)
    for j = 1:length(q_lvl)
        if m_dt(i) < q_lvl(j)
            % Check if sample is closer to upper or lower quantization lvl
            if abs(m_dt(i) - q_lvl(j-1)) < abs(m_dt(i) - q_lvl(j))
                indi(i) = j-1;
            else
                indi(i) = j;
            end
            break;
        end
    end
end

% Convert quantized levels to binary representation
bin_rep = dec2bin(indi);

% Display first 10 samples of binary representation
fprintf('PCM output: \n')
disp(bin_rep(1:10,:))

%% TASK 2: DELTA MODULATION
t_end = 2; % time length
fs = 4*200; % sampling frequency

% Discrete time vector
t = linspace(0.005, t_end, t_end*fs);

% Create signal m(t) in discrete time points
m_dt = -2*cos(200*pi*t) + sin(50*pi*t);

%finding the difference between samples m[k]-m[k-1]
for i = 2:length(m_dt)
    dq(i-1) = m_dt(i) - m_dt(i-1);
end

for i = 1:length(dq) %finding 1-bit DPCM
   if dq(i) > 0
      binart_out(i) = 1;
   else
      binart_out(i) = 0;
   end
end

fprintf('\nDelta Modulation output: \n')
disp(binart_out(1:20)) % print the binary representation of the first 20 samples
%%


