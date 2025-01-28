clear; clc;

% Simulation parameters
n = 100; % Number of particles
L = 25; % Length of the simulation box
mass = 1; % Mass of each particle
sigma = 1; % Lennard-Jones parameter
epsilon = 1; % Lennard-Jones parameter
dt = 0.002; % Time step (reduced for stability)
steps = 150; % Number of simulation steps
Volume = L^3;

% Initialize positions and velocities
positions = L * rand(n, 3); 
velocities = randn(n, 3); 
velocities = velocities - mean(velocities);

% Initialize arrays for storage
forces = zeros(n, 3);
Energy_storage = zeros(steps, 1);
KEnergy_storage = zeros(steps, 1);
Temperature_storage = zeros(steps, 1);
Pressure_storage = zeros(steps, 1);

% Verlet Initialization
prev_positions = positions - velocities * dt; 

for step = 1:steps
    % Reset forces and energies
    forces(:) = 0;
    PotentialEnergy = 0;
    KineaticEnergy = 0;

    % Compute forces and potential energy using Lennard-Jones potential
    for i = 1:n
        for j = i+1:n
            % Compute distance vector and distance
            r_vec = positions(i, :) - positions(j, :);
            % Apply minimum image convention for periodic boundary conditions
            r_vec = r_vec - L * round(r_vec / L);
            r = norm(r_vec);

            % Apply cutoff distance for efficiency
            if r < 2.5 * sigma
                r6 = (sigma / r)^6;
                r12 = r6^2;
                force_mag = 24 * epsilon * (2 * r12 - r6) / r^2;

                % Update forces
                forces(i, :) = forces(i, :) + force_mag * r_vec / r;
                forces(j, :) = forces(j, :) - force_mag * r_vec / r; % Newton's Third Law

                % Compute potential energy
                PotentialEnergy = PotentialEnergy + 4 * epsilon * (r12 - r6);
            end
        end

        % Compute kinetic energy
        v = norm(velocities(i, :));
        KineaticEnergy = KineaticEnergy + 0.5 * mass * v^2;
    end

    % Calculate total energy, temperature, and pressure
    Total_Energy = PotentialEnergy + KineaticEnergy;
    Temperature = (2 * KineaticEnergy) / (3 * n);
    Pressure = (n * Temperature / Volume) + (1 / (3 * Volume)) * sum(sum(forces .* positions));

    % Store values
    Energy_storage(step) = Total_Energy;
    KEnergy_storage(step) = KineaticEnergy;
    Temperature_storage(step) = Temperature;
    Pressure_storage(step) = Pressure;

    % Verlet integration
    new_positions = 2 * positions - prev_positions + (forces / mass) * dt^2;
    new_velocities = (new_positions - positions) / dt;

    % Apply periodic boundary conditions
    new_positions = mod(new_positions, L);

    % Update positions and velocities
    prev_positions = positions;
    positions = new_positions;
    velocities = new_velocities;
end

% Plot results
x = 1:steps;
figure;

subplot(3, 1, 1);
plot(x, Energy_storage, 'r', 'LineWidth', 1.5);
hold on;
plot(x, KEnergy_storage, 'b', 'LineWidth', 1.5);
legend('Total Energy', 'Kinetic Energy');
title('Energy vs. Time Steps');
xlabel('Time Step');
ylabel('Energy');
grid on;

subplot(3, 1, 2);
plot(x, Temperature_storage, 'g', 'LineWidth', 1.5);
title('Temperature vs. Time Steps');
xlabel('Time Step');
ylabel('Temperature');
grid on;

subplot(3, 1, 3);
plot(x, Pressure_storage, 'm', 'LineWidth', 1.5);
title('Pressure vs. Time Steps');
xlabel('Time Step');
ylabel('Pressure');
grid on;

disp('Simulation Complete!');
