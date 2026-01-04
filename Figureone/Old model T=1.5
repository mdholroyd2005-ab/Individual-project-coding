%The below code is a 1D simulation of the Ising Model using Metropolis
%Hastings. It has been adapted from the code for question 1.
%sets temperature
rng(1)
T = 1.5;
N = 100;

%Setup model
n_samples = 5e6;
burn_in=1e4;
sample = zeros(n_samples,1);
accepted = zeros(n_samples,1);

% Track the middle spin
middle_spin_ = round(N/2);
spin_time_series = zeros(n_samples, 1);

%Generates a list of +1,-1s randomly
x0 = sign(randn(N,1));

%Generates the initial value and burns it in to remove effects of bad
%initial value
for i=1:burn_in
    [x0,~]=update_x(x0,T);
end

current_x = x0;

%generate the sample by iterating on the markov chain
sample(1) = mean(x0);
spin_time_series(1) = current_x(middle_spin_);

for i = 2:n_samples
    [current_x,accepted(i)] = update_x(current_x,T);
    sample(i) = mean(current_x);
    
    % Record the state of the middle spin at this step ---
    spin_time_series(i) = current_x(middle_spin_);
end

% Plot1: Histogram
figure(1); 
nbins = 50;
%plots histogram from samples
h=histogram(sample, nbins,'Normalization', 'pdf', 'DisplayName', 'Generated Sample');
ylim([0,max(h.Values)+0.5]);
xlim([-1,1]);
ylabel('Frequency');
xlabel('magnetism');
hold on;
t = title("1D Ising model simulation (T = " + string(T) + ")");
t.Interpreter= 'latex';  

% Plot2 - Single electron time series
figure(2); 
% Plot the first 10000 steps so transitions are visible
steps_to_view = 10000; 
stairs(1:steps_to_view, spin_time_series(1:steps_to_view), 'LineWidth',2 );
ylim([-1.2, 1.2]);
xlabel('Time Step');
ylabel('Spin State');
title(['Time Series of Spin at Index ' num2str(middle_spin_) ' (First ' num2str(steps_to_view) ' steps)']);
yticks([-1 1]);
yticklabels({'-1', '+1'});
grid on;

function [next_x, accepted] = update_x(x, T)
    %Choose a uniform number in {1...N} and flip sigma
    N = length(x);
    idx = randi(N); 
    proposed = x;
    proposed(idx) = -x(idx); % Flip spin
    %Energy Calculation
    % Energy E = -sum(sigma_i * sigma_{i+1})
    % Calculate energy for the current and proposed states
    % (Using line boundary conditions rather than ring)
    E_current = -sum(x(1:end-1) .* x(2:end));
    E_proposed = -sum(proposed(1:end-1) .* proposed(2:end));
    %Acceptance probabiility
    % Ratio = exp(-E'/T) / exp(-E/T) = exp(-(E' - E)/T) - constant cancels
    % here
    delta_E = E_proposed - E_current;
    accept_prob = min(1, exp(-delta_E / T));
    %Metropolis step
    if rand() <= accept_prob
        next_x = proposed;
        accepted = 1;
    else
        next_x = x;
        accepted = 0;
    end
end
