#include<iostream>
#include<cstdlib>
#include<vector>
#include<math.h>

int n =500;
int L = 10;
int steps = 500;
int dim =2;
double mass = 100;
double sigma = 1;
double epsilon = 0.001;
double dt = 0.002;
double Potential_Energy;
double KineaticEnergy;
//double r_vec_stg[n][n][dim];


struct Particles {
    double position[2];
    double next_position[2];
    double prev_position[2];
    double velocity[2];
    double force[2];
    double PEnergy;
    double KEnergy;

};

std::vector<Particles> particle(n);

void intital_positions(){
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < dim; j++) {
            particle[i].position[j] = L*(static_cast<double>(rand())/RAND_MAX);
            particle[i].velocity[j] = (static_cast<double>(rand())/RAND_MAX);
            particle[i].prev_position[j] = particle[i].position[j] - particle[i].velocity[j]*dt;
            particle[i].force[j] = 0;

        }
    }
}

/*void distance_between_particle() {
    for (int k =0;k<dim;k++) {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (i !=j) {
                    r_vec_stg[i][j][k] = particle[i].position[k] - particle[j].position[k];
                }

            }
        }
    }
}*/
void Force_onparticles() {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < dim; j++) {
            particle[i].force[j] = 0;
            particle[i].KEnergy = 0;
            particle[i].PEnergy = 0;
        }
    }
    KineaticEnergy = 0;
    Potential_Energy = 0;
    for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++) {
            double r_vec[dim];
            double rsq = 0;
            for (int k=0;k<dim;k++) {
                r_vec[k] = particle[i].position[k] - particle[j].position[k];
                rsq = rsq + r_vec[k]*r_vec[k];
            }
            double r = sqrt(rsq);
            if (r==0) {
                continue;
            }

            if (r<2.5*sigma) {
                double r6 = pow((sigma/r),6);
                double r12 = r6*r6;
                //std::cout<<r12<<" "<<r6<<std::endl;

                for (int k=0;k<dim;k++) {
                    particle[i].force[k] = particle[i].force[k] + (24 * epsilon*(2*r12-r6)/(r*r));
                }
                particle[i].PEnergy = particle[i].PEnergy + 4 * epsilon*(r12-r6);

            }
            else {
                continue;
            }
        }
        Potential_Energy += particle[i].PEnergy;
        double velocity = 0;
        for (int k=0;k<dim;k++) {
            velocity = velocity + particle[i].velocity[k]*particle[i].velocity[k];
        }
        particle[i].KEnergy = .5*(mass*velocity);
        KineaticEnergy = particle[i].KEnergy + KineaticEnergy;
    }
}

void update_ofparticle_position() {
    for (int i=0; i<n; i++) {
        for (int j=0; j<dim; j++) {
            particle[i].next_position[j] = 2*particle[i].position[j] - particle[i].prev_position[j] + (particle[i].force[j]/mass)*dt*dt;

            if (particle[i].next_position[j]<0) {
                particle[i].next_position[j] = fmod(particle[i].next_position[j], L) + L;
                particle[i].velocity[j] = -particle[i].velocity[j];
            }
            if (particle[i].next_position[j]>L) {
                particle[i].next_position[j] = fmod(particle[i].next_position[j], L) ;
                particle[i].velocity[j] = -particle[i].velocity[j];
            }
            particle[i].prev_position[j] = particle[i].position[j];
            particle[i].position[j] = particle[i].next_position[j];
        }
        //std::cout<<particle[1].next_position[i]<<std::endl;
    }
}

double compute_temperature() {
    double total_kinetic_energy = 0.0;
    for (int i = 0; i < n; i++) {
        double velocity_squared = 0.0;
        for (int j = 0; j < dim; j++) {
            velocity_squared += particle[i].velocity[j] * particle[i].velocity[j];
        }
        total_kinetic_energy += 0.5 * mass * velocity_squared;
    }

    // Boltzmann constant in reduced units is 1
    double temperature = (2.0 * total_kinetic_energy) / (dim * n);
    return temperature;
}

double compute_pressure(double temperature) {
    double virial_sum = 0.0;
    double volume = pow(L, dim);

    for (int i = 0; i < n; i++) {
        for (int j = 1; j < n; j++) {
            double r_vec[dim] = {0};
            double rsq = 0.0;

            for (int k = 0; k < dim; k++) {
                r_vec[k] = particle[i].position[k] - particle[j].position[k];
                rsq += r_vec[k] * r_vec[k];
            }


            double r = sqrt(rsq);
            if (r==0) {
                continue;
            }
            if (r < 2.5 * sigma) {
                double r6 = pow((sigma / r), 6);
                double r12 = r6 * r6;
                double force_mag = 24 * epsilon * (2 * r12 - r6) / (r * r);

                for (int k = 0; k < dim; k++) {
                    virial_sum += r_vec[k] * force_mag * r_vec[k] / r;
                }
            }
        }
    }

    double pressure = (n * temperature / volume) + (virial_sum / (dim * volume));
    return pressure;
}


int main() {
    intital_positions();

    for (int step = 0; step < steps; step++) {
        Force_onparticles();
        update_ofparticle_position();

        double temperature = compute_temperature();
        double pressure = compute_pressure(temperature);

        std::cout << "Step: " << step
                  << ", Temperature: " << temperature
                  << ", Pressure: " << pressure
                  << ", Potential Energy: " << Potential_Energy
                  << std::endl;
    }

    return 0;
}