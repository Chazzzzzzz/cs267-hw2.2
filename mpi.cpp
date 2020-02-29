#include "common.h"
#include <mpi.h>
#include <math>
#include <vector>

using namespace std;
// Put any static global variables here that you will use throughout the simulation.
typedef std::vector<particle_t*> bin_t;

double size;
int bin_row_count;
int bin_count;
double bin_size;
int bins_per_proc;
bin_t* bins;
std::vector<int> local_binID;
std::vector<particle_t> particles_to_send;
particle_t* particles;
int* migrate_size;
int* disp_size;
int proc_row_count;
int proc_count;

int inline get_bin_id(particle_t& particle) {
    int x, y;
    y = particle.y / bin_size;
    x = particle.x / bin_size;
    if (x == bin_row_count) {
        x--;
    }
    if (y == bin_row_count) {
        y--;
    }
    return y * bin_row_count + x;
}

int inline get_proc_id(int bin_id) {
    return bin_id / bins_per_proc;
}

void assign_bins(particle_t* parts, int num_parts, int rank) {
    for (int i = 0; i < num_parts; i++) {
        int bin_id = get_bin_id(parts[i]);
        int proc_id = get_proc_id(bin_id);
        if (proc_id == rank) {
            bins[bin_id].push_back(&parts[i]);
        }
    }
}

void rebin(particle_t* parts, int num_parts, int rank, int num_procs) {
    for (int binID: local_binID) {
        for (particle_t* p: bins[binID]) {
            int new_binID = get_bin_id(*p);
            int new_procID = get_proc_id(new_binID);
            if (new_procID != rank) {
                bins[binID].erase(remove(bins[binID].begin(), bins[binID].end(), p), bins[binID].end());
                particles_to_send.push_back(*p);
            } else if(new_binID != binID) {
                bins[binID].erase(remove(bins[binID].begin(), bins[binID].end(), p), bins[binID].end());
                bins[new_binID].push_back(p);
            }
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    int particles_to_send_size = particles_to_send.size();
    MPI_Allgather(&particles_to_send_size, 1, MPI::INT, migrate_size, 1, MPI::INT, MPI_COMM_WORLD);

    disp_size[0] = 0;
    for (int i = 1; i < n_proc; i++) {
        disp_size[i] = disp_size[i-1] + migrate_size[i-1];
    }

    MPI_Allgatherv(&particles_to_send[0], particles_to_send.size(), PARTICLE, particles, migrate_size, disp_size, PARTICLE, MPI_COMM_WORLD);
    
    int overall_migrate_num = 0;
    for (int i = 0; i < num_procs; i++) {
        overall_migrate_num += migrate_size[i];
    }

    for (int i = 0; i < overall_migrate_num; i++) {
        int bin_id = get_bin_id(particles[i]);
        int proc_id = get_proc_id(bin_id);
        if (proc_id == rank) {
            bins[bin_id].push_back(&particles[i]);
        }
    }
}
int max_partitions(int process_count) {
    int partitions = 1;
    for (int i = 1; i <= std::fmin(process_count, sqrt(bin_count)); i++) {
        if (bin_count % i == 0) {
            if (bin_count / i <= process_count) {
                return bin_count / i;
            } else {
                partitions = i;
            }
        }
    }
    return partitions;
}

void get_local_binID(int rank) {
    for (int i = 0; i < bins_per_proc; i++) {
        local_binID.push_back(rank * bins_per_proc + i);
    }
}

void init_simulation(particle_t* parts, int num_parts, double size_, int rank, int num_procs) {
	// You can use this space to initialize data objects that you may need
	// This function will be called once before the algorithm begins
	// Do not do any particle simulation here
    size = size_;
    bin_row_count = size / cutoff;
    bin_count = bin_row_count * bin_row_count;
    bin_size = size / bin_row_count;
    bins = new bin_t[bin_row_count * bin_row_count];
    // get bins_per_proc
    proc_count = max_partitions(num_procs);
    bins_per_proc = bin_count / proc_count;
    if (bin_row_count > bins_per_proc) {
        proc_row_count = bin_row_count / bins_per_proc;
    } else {
        proc_row_count = 1;
    }
    // assign bins
    get_local_binID(rank);
    assign_bins(parts, num_parts, rank);
    // init particles
    particles = (particle_t*) malloc(num_parts * sizeof(particle_t));
    migrate_size = (int*) malloc(num_procs * sizeof(int));
    disp_size = (int*) malloc(num_procs * sizeof(int));
}

void apply_force(particle_t& particle, particle_t& neighbor) {
    // Calculate Distance
    double dx = neighbor.x - particle.x;
    double dy = neighbor.y - particle.y;
    double r2 = dx * dx + dy * dy;

    // Check if the two particles should interact
    if (r2 > cutoff * cutoff)
        return;

    r2 = fmax(r2, min_r * min_r);
    double r = sqrt(r2);

    // Very simple short-range repulsive force
    double coef = (1 - cutoff / r) / r2 / mass;
    particle.ax += coef * dx;
    particle.ay += coef * dy;
}

void move(particle_t& p, double size) {
    // Slightly simplified Velocity Verlet integration
    // Conserves energy better than explicit Euler method
    p.vx += p.ax * dt;
    p.vy += p.ay * dt;
    p.x += p.vx * dt;
    p.y += p.vy * dt;

    // Bounce from walls
    while (p.x < 0 || p.x > size) {
        p.x = p.x < 0 ? -p.x : 2 * size - p.x;
        p.vx = -p.vx;
    }

    while (p.y < 0 || p.y > size) {
        p.y = p.y < 0 ? -p.y : 2 * size - p.y;
        p.vy = -p.vy;
    }
}


bool inline has_up_bin(int bin_id) {
    return bin_id - bin_row_count > -1;
}
bool inline has_down_bin(int bin_id) {
    return bin_id + bin_row_count < bin_count;
}
bool inline has_left_bin(int bin_id) {
    return bin_id % bin_row_count != 0;
}
bool inline has_right_bin(int bin_id) {
    return bin_id % bin_row_count != bin_row_count - 1;
}

bool inline has_up_proc(int proc_id) {
    return proc_id - proc_row_count > -1;
}
int inline get_up_proc(int proc_id) {
    return proc_id - proc_row_count;
}
bool inline has_down_proc(int proc_id) {
    return proc_id + proc_row_count < proc_count;
}
int inline get_down_proc(int proc_id) {
    return proc_id + proc_row_count;
}
bool inline has_left_proc(int proc_id) {
    return proc_id % proc_row_count != 0;
}
int inline get_left_proc(int proc_id) {
    return proc_id - 1;
}
bool inline has_right_proc(int proc_id) {
    return proc_id % proc_row_count != proc_row_count - 1;
}
int inline get_right_proc(int proc_id) {
    return proc_id + 1;
}

std::vector<int> get_up_border_bin_ids(int proc_id) {
    std::vector<int> res;
    if (bin_row_count > bins_per_proc) {
        for (int i = 0; i < bins_per_proc; i++) {
            res.push_back(proc_id * bins_per_proc + i);
        }
    } else {
        for (int i = 0; i < bin_row_count; i++) {
            res.push_back(proc_id * bins_per_proc + i);
        }
    }
    return res;
}
std::vector<int> get_down_border_bin_ids(int proc_id) {
    std::vector<int> res;
    if (bin_row_count > bins_per_proc) {
        for (int i = 0; i < bins_per_proc; i++) {
            res.push_back(proc_id * bins_per_proc + i);
        }
    } else {
        for (int i = 0; i < bin_row_count; i++) {
            res.push_back((proc_id + 1) * bins_per_proc - i - 1);
        }
    }
    return res;
}
int get_left_border_bin_id(int proc_id) {
    if (bin_row_count > bins_per_proc) {
        return proc_id * bins_per_proc;
    } else {
        return -1;
    }
}
int get_right_border_bin_id(int proc_id) {
    if (bin_row_count > bins_per_proc) {
        return (proc_id + 1) * bins_per_proc - 1;
    } else {
        return -1;
    }
}

void simulate_one_step(particle_t* parts, int num_parts, double size, int rank, int num_procs) {
    std::vector<MPI_Request> requests;
    if (has_up_proc(rank)) {
        std::vector<int> up_border = get_up_border_bin_ids(rank);
        std::vector<particle_t> to_send;
        for (auto bin_id : up_border) {
            for (auto part : bins[bin_id]) {
                to_send.push_back(*part);
            }
        }
        MPI_Request req;
        requests.push_back(req);
        MPI_Isend(&to_send[0], to_send.size(), PARTICLE, get_up_proc(rank), 0, MPI_COMM_WORLD, &req);
    }
    if (has_left_proc(rank)) {
        int left_border = get_left_border_bin_id(rank);
        std::vector<particle_t> to_send;
        for (auto part : bins[left_border]) {
            to_send.push_back(*part);
        }
        MPI_Request req;
        requests.push_back(req);
        MPI_Isend(&to_send[0], to_send.size(), PARTICLE, get_left_proc(rank), 0, MPI_COMM_WORLD, &req);
    }
    if (has_down_proc(rank)) {
        std::vector<int> down_border = get_down_border_bin_ids(rank);
        std::vector<particle_t> to_send;
        for (auto bin_id : down_border) {
            for (auto part : bins[bin_id]) {
                to_send.push_back(*part);
            }
        }
        MPI_Request req;
        requests.push_back(req);
        MPI_Isend(&to_send[0], to_send.size(), PARTICLE, get_down_proc(rank), 0, MPI_COMM_WORLD, &req);
    }
    if (has_right_proc(rank)) {
        int right_border = get_right_border_bin_id(rank);
        std::vector<particle_t> to_send;
        for (auto part : bins[right_border]) {
            to_send.push_back(*part);
        }
        MPI_Request req;
        requests.push_back(req);
        MPI_Isend(&to_send[0], to_send.size(), PARTICLE, get_right_proc(rank), 0, MPI_COMM_WORLD, &req);
    }
    if (has_up_proc(rank) && has_left_proc(rank)) {
        int left_border = get_left_border_bin_id(rank);
        std::vector<particle_t> to_send;
        for (auto part : bins[left_border]) {
            to_send.push_back(*part);
        }
        MPI_Request req;
        requests.push_back(req);
        MPI_Isend(&to_send[0], to_send.size(), PARTICLE, get_up_proc(rank) - 1, 0, MPI_COMM_WORLD, &req);
    }
    if (has_up_proc(rank) && has_right_proc(rank)) {
        int right_border = get_right_border_bin_id(rank);
        std::vector<particle_t> to_send;
        for (auto part : bins[right_border]) {
            to_send.push_back(*part);
        }
        MPI_Request req;
        requests.push_back(req);
        MPI_Isend(&to_send[0], to_send.size(), PARTICLE, get_up_proc(rank) + 1, 0, MPI_COMM_WORLD, &req);
    }
    if (has_down_proc(rank) && has_left_proc(rank)) {
        int left_border = get_left_border_bin_id(rank);
        std::vector<particle_t> to_send;
        for (auto part : bins[left_border]) {
            to_send.push_back(*part);
        }
        MPI_Request req;
        requests.push_back(req);
        MPI_Isend(&to_send[0], to_send.size(), PARTICLE, get_down_proc(rank) - 1, 0, MPI_COMM_WORLD, &req);
    }
    if (has_down_proc(rank) && has_right_proc(rank)) {
        int right_border = get_right_border_bin_id(rank);
        std::vector<particle_t> to_send;
        for (auto part : bins[right_border]) {
            to_send.push_back(*part);
        }
        MPI_Request req;
        requests.push_back(req);
        MPI_Isend(&to_send[0], to_send.size(), PARTICLE, get_down_proc(rank) + 1, 0, MPI_COMM_WORLD, &req);
    }
    /*
    for (auto req : requests) {
        MPI_Status status;
        MPI_Wait(&req, &status);
    }
    */
    if (has_up_proc(rank)) {
        MPI_Status status;
        particle_t recv_buff[10];
        MPI_Recv(&recv_buff, num_parts, PARTICLE, get_up_proc(rank), 0, MPI_COMM_WORLD, &status);
    }
}

void gather_for_save(particle_t* parts, int num_parts, double size, int rank, int num_procs) {
    // Write this function such that at the end of it, the master (rank == 0)
    // processor has an in-order view of all particles. That is, the array
    // parts is complete and sorted by particle id.


}
