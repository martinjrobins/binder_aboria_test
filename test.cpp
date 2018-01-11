#include "Aboria/src/Aboria.h"

const size_t N = 100;
const unsigned seed = 0;
const bool write_out = false;

/*
 * set parameters
 */
const int timesteps = 1e3;
const double r_cut =  std::sqrt(3.0/N);
const double c = 1e-3;

/*
 * Create a 2d particle container type with one 
 * additional variable "velocity", represented 
 * by a 2d double vector
 */
ABORIA_VARIABLE(velocity,double2,"velocity")
typedef Particles<std::tuple<velocity>,2> container_type;
typedef typename container_type::position position;

/*
 * create a particle set with size N
 */
container_type particles(N);

std::uniform_real_distribution<double> uni(0,1);
std::default_random_engine gen(seed);
for (int i = 0; i < N; ++i) {
    /*
     * set a random position, and initialise velocity
     */
    get<position>(particles)[i] = double2(uni(gen),uni(gen));
    get<velocity>(particles)[i] = double2(0,0);
}

/*
 * initiate neighbour search on a periodic 2d domain 
 * of side length 1
 */
particles.init_neighbour_search(double2(0,0),
                                double2(1,1),
                                r_cut,
                                bool2(true,true));


/*
 * create symbols and labels in order to use 
 * the Level 3 API
 */
Symbol<position> p;
Symbol<velocity> v;
Label<0,container_type> i(particles);
Label<1,container_type> j(particles);

/*
 * dx is a symbol representing the difference in 
 * positions of particles i and j.
 */
auto dx = create_dx(i,j);

/*
 * sum is a symbolic function that sums 
 * a sequence of 2d vectors
 */
Accumulate<std::plus<double2> > sum;

/*
 * perform timestepping
 */
for (int io = 0; io < timesteps; ++io) {

    /*
     * leap frog integrator
     */
    v[i] += c*sum(j, norm(dx)<r_cut && norm(dx)>0,
                -exp(-norm(dx))*dx/norm(dx)
                );
    p[i] += v[i];
}
