%The below code is a 1D simulation of the Ising Model using Metropolis
%Hastings. It has been adapted such that each spin interacts with all other
%spins, interaction strength varies with distance. This code now also
%tracks the 50th spin and creates a time series graph.
rng(2)
%sets temperature
T = 1.5;
N = 100;
a = 1.9;

%initialises data
n_samples = 5e6;
burn_in=1e4;
sample = zeros(n_samples,1);
accepted = zeros(n_samples,1);
spin50history = zeros(n_samples, 1);

%generates a list of +1,-1s randomly
x0 = sign(randn(N,1));

%generates the initial value and burns it in to remove effects of bad
%initial value
for i=1:burn_in
    [x0,~]=update_x(x0,T,a);
end
current_x = x0;
%generate the sample by iterating on the markov chain
sample(1) = mean(x0);
for i = 2:n_samples
    [current_x,accepted(i)] = update_x(current_x,T,a);
    sample(i) = mean(current_x);
    spin50history(i) = current_x(50);
end
nbins = 50;

%plots histogram from samples
h=histogram(sample, nbins,'Normalization', 'pdf', 'DisplayName', 'Generated Sample');
ylim([0,max(h.Values)+0.5]);
xlim([-1,1]);
ylabel('Frequency');
xlabel('Magnetism');
hold on;
t = title("1D Long-Range Ising (a = " + string(a) + ",T = " + string(T) + ")");
t.Interpreter= 'latex';  

%Plots time series of spin 50


figure;
plot(spin50history(1:2000), 'LineWidth', 1.5);
ylim([-1.5, 1.5]);
yticks([-1, 1]);
grid on;
ylabel('Spin Value');
xlabel('Time');
t2 = title("Time-Evolution of Spin 50 (a = " + string(a) + ",T = " + string(T) + ")");
t2.Interpreter = 'latex';

function [next_x, accepted] = update_x(x, T, a)
    %Choose a uniform number in {1...N} and flip sigma
    N = length(x);
    idx = randi(N); 

    %Energy Calculation changes here
    %spins interact with all others

    interaction_sum = 0;
    for j = 1:N
        if j ~= idx
            distance = abs(idx - j); % Calculate distance |i - j|
            % We weight interaction by distance to the power of a
            interaction_sum = interaction_sum + (x(j) / (distance^a));
        end
    end

    %Acceptance probability 
    % Ratio = exp(-E'/T) / exp(-E/T) = exp(-(E' - E)/T)
    % Adapted: Energy change for long-range single-flip
    delta_E = 2 * x(idx) * interaction_sum;

    accept_prob = min(1, exp(-delta_E / T));
    
    %Metropolis step
    if rand() <= accept_prob
        x(idx) = -x(idx);
        next_x = x;
        accepted = 1;
    else
        next_x = x;
        accepted = 0;
    end
end
