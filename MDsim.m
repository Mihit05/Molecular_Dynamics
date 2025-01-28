
clear; clc;


n = 500; % Number of particles
L = 10; % Length of the simulation box
mass = 1; % Mass of each particle
sigma = 1; % Lennard-Jones parameter
epsilon = 1; % Lennard-Jones parameter
dt = 0.002; % Time step
steps = 5000; % Number of simulation steps
k = 1;
Volume = L*L*L;
PotentialEnergy = 0;
KineaticEnergy = 0;
% Initialize positions and velocities randomly within the box
positions = L * rand(n, 3); % Random positions
velocities = rand(n, 3); % Random velocities
velocities = velocities - mean(velocities); % Remove net momentum

forces = zeros(n, 3);
 Energy_storage = zeros(steps,1);
KEnergy_storage = zeros(steps,1);
Temperature_storage = zeros(steps,1);
 r_storage = zeros(steps,1);
R_val = zeros(steps,1);
% Set up the figure for visualization
figure;
hold on;
box on;
axis([0 L 0 L 0 L]); % Define the simulation box dimensions
xlabel('X'); ylabel('Y'); zlabel('Z');
grid on;

% Create scatter plot for particles
particle_plot = scatter3(positions(:, 1), positions(:, 2), positions(:, 3), ...
    100, 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r');
view(3);
% Verlet Initialization
prev_positions = positions - velocities * dt; % Backward step for Verlet

for step = 1:steps
    % Reset forces
    forces(:) = 0;
    
    r = 0;

    % Compute forces using Lennard-Jones potential
    for i = 1:n
        for j = i+1:n
            % Compute distance vector and distance
            r_vec = positions(i, :) - positions(j, :);
            r = norm(r_vec);
            v = norm(velocities);
            if r_vec>10
                continue;
            end
                            
            % Lennard-Jones force
            if r < 2.5 * sigma % Cutoff for computational efficiency
                r6 = (sigma / r)^6;
                r12 = r6^2;
                force_mag = 24 * epsilon * (2 * r12 - r6) / r^2;

                % Update forces
                forces(i, :) = forces(i, :) + force_mag * r_vec / r;
                forces(j, :) = forces(j, :) - force_mag * r_vec / r; % Newton's Third Law
    
            end
        end
        
        KineaticEnergy = KineaticEnergy + 0.5*mass*(v*v);
      
        PotentialEnergy = PotentialEnergy + 4*epsilon*((sigma/r)^12 - (sigma/r)^6);
          Total_Energy = PotentialEnergy + KineaticEnergy;
          Temperature = (2*KineaticEnergy)/(3*n);
          Pressure = (n*Temperature/(L*L*L));

    end
    %disp(Temperature);
    disp(Total_Energy);

    % disp(PotentialEnergy);
    % disp(KineaticEnergy);

    Energy_storage(step,1) = PotentialEnergy;
    KEnergy_storage(step,1) = KineaticEnergy;
    TotalEnergy_stg(step,1) = Total_Energy;
    Temperature_storage(step,1) = Temperature;
    %disp(TotalEnergy_stg);
    %fprintf("%d      %d     %d     %d\n", PotentialEnergy,KineaticEnergy, Total_Energy, Temperature );
    r_storage(step,1) = r;
    R_val(steps,1) = r;

    % Verlet integration
            new_positions = 2 * positions - prev_positions + (forces / mass) * dt^2;
            new_velocity = (new_positions - positions)/dt;

    % Apply rebounding boundary conditions
    for i = 1:n
        for dim = 1:3
            % Check if the particle crosses the lower boundary
            if new_positions(i, dim) < 0
                new_positions(i, dim) = -new_positions(i, dim); % Reflect position 
                velocities(i, dim) = -velocities(i, dim); % Invert velocity
            end

            % Check if the particle crosses the upper boundary
            if new_positions(i, dim) > L
                new_positions(i, dim) = 2 * L - new_positions(i, dim); % Reflect position
                velocities(i, dim) = -velocities(i, dim); % Invert velocity
            end
        end
    end


    % Update positions and previous positions
    prev_positions = positions;
    positions = new_positions;
    velocities = new_velocity;


    % Update the particle positions in the scatter plot
    set(particle_plot, 'XData', positions(:, 1), ...
                       'YData', positions(:, 2), ...
                       'ZData', positions(:, 3));

    % Update the title with the current step
    title(sprintf('Molecular Dynamics Simulation with Rebounding Boundaries\nStep: %d / %d', step, steps));

    % Pause for visualization
    pause(0.01);
end
disp('Simulation Complete!');