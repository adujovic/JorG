#include "solver.h"

int solver(char _basis[],char _supercell[],char _flippable[], size_t reference, size_t unique_flips, size_t ansatz){
#ifndef _SITESNUMBER
#define _SITESNUMBER 64
#endif
    constexpr size_t SITESNUMBER = _SITESNUMBER;
    int rank = -1;
    size_t unique_flips_per_iteration;
#ifdef _MPI
    char solution[SITESNUMBER];
    int nprocs;
    bool solution_found;
    MPI_Status status;

    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    unique_flips_per_iteration = std::ceil(unique_flips/static_cast<float>(nprocs));
#else
    unique_flips_per_iteration = unique_flips;
#endif //_MPI
// *****************************************************************************************
// **************************************   READ   *****************************************
// *****************************************************************************************
    aux::System<typename ising::IsingModel<SITESNUMBER>::VectorType> system;
    aux::read_basis(system.basis,_basis);
    aux::read_supercell(system.supercell,_supercell);
    aux::read_flippables(system.flippable,_flippable);
#ifndef _QUIET
    if (rank < 1) aux::greetings();
#endif
// *****************************************************************************************
// *******************************  INITIALIZE ENGINE    ***********************************
// *****************************************************************************************
    std::random_device randomDevice{};
    std::mt19937 generator{randomDevice()};
    std::normal_distribution<> gauss{0.0,1.0};
    ising::IsingModel<SITESNUMBER> model;
    model.set_basis(system.basis);
    model.set_supercell(system.flippable);
    model.set_reference(reference);
// *****************************************************************************************
// ********************************   MAP INTERACTIONS   ***********************************
// *****************************************************************************************
    double decayCoeff;
    aux::map_geometry(system.flippable,decayCoeff);
    std::bitset<SITESNUMBER> mask(std::string(SITESNUMBER,'0'));
    for (const auto& atom : system.flippable) mask.set(std::get<0>(atom));
    mask.reset(reference);
// *****************************************************************************************
// ********************************   CHECK MODEL   ****************************************
// *****************************************************************************************
#ifdef _VERBOSE
    if (rank < 1) { 
        std::cout<<"               ";
        for(int i=1; i<(SITESNUMBER-reference); ++i) std::cout<<" ";
        std::cout<<"| reference"<<std::endl;
        std::cout<<"Mask:          "<<mask<<std::endl;
    }
#endif
#ifdef _LOWERLIMIT
    size_t lowerlimit = _LOWERLIMIT;
#else
    size_t lowerlimit = 0.15*mask.count(); // doesn't allow less than 15% flipps
#endif
#ifdef _UPPERLIMIT
    size_t upperlimit = _UPPERLIMIT;
#else
    size_t upperlimit = 0.85*mask.count(); // doesn't allow more than 85% flipps
#endif
    upperlimit = upperlimit < 2U         ? 2U : upperlimit;
    lowerlimit = lowerlimit > upperlimit ? 1U : lowerlimit;
    lowerlimit = lowerlimit < 1U         ? 1U : lowerlimit;
    std::ofstream ostrm("best.flips", std::ios::out);
    std::unordered_set<std::bitset<SITESNUMBER>> solutions;
    size_t iteration = 0U;
    size_t nsolutions = 0U;
    for(auto n = 4.0; n<1024.0; n*=2){
        model.reset();
        auto d = aux::get_interactions(system,n*decayCoeff,ansatz);
        model.add_interaction(d);
#ifdef _VERBOSE
        if (rank < 1)
            std::cout<<"Limits are set to be: ["<<lowerlimit<<","<<upperlimit<<']'<<std::endl;
#endif
        for(size_t m = 0U; m<unique_flips_per_iteration; ++m){
            model.randomize_state();
            auto x = model.run(&mask);

            if(x.count() < lowerlimit || x.count() > upperlimit) {
#ifdef _MPI
                solution_found = false;
            } else solution_found = true;
            if(rank > 0) MPI_Send(&solution_found, 1, MPI_C_BOOL, 0, 0, MPI_COMM_WORLD);
#else
                continue;
            }
#endif // _MPI
            if(x.count() > 0.5*mask.count()) x.flip();
#ifdef _MPI
            if (rank > 0) {
                if (solution_found) {
                    MPI_Send(x.to_string().c_str(), SITESNUMBER, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
                }
            } else {
                if (solution_found) {
                    solutions.insert(x);
#ifndef _QUIET
                    std::cout<<"("<<iteration<<")\t"<<n*decayCoeff<<"\t";
                    aux::print_state(x,reference,mask);
#endif
                }
                for (int i = 1; i < nprocs; ++i) {
                    MPI_Recv(&solution_found, 1, MPI_C_BOOL, i, 0, MPI_COMM_WORLD, &status);
                    if (solution_found) {
                        MPI_Recv(solution, SITESNUMBER, MPI_CHAR, i, 0, MPI_COMM_WORLD, &status);
                        x = std::bitset<SITESNUMBER>(std::string(solution));
                        solutions.insert(x);
#ifndef _QUIET
                        std::cout<<"("<<iteration<<")\t"<<n*decayCoeff<<"\t";
                        aux::print_state(x,reference,mask);
#endif
                    }
                }
            }
#else
            solutions.insert(x);
#ifndef _QUIET
            std::cout<<"("<<iteration<<")\t"<<n*decayCoeff<<"\t";
            aux::print_state(x,reference,mask);
#endif
#endif // _MPI
            if (rank < 1) { nsolutions = solutions.size(); }
#ifdef _MPI
            MPI_Bcast(&nsolutions, 1, MPI_UINT64_T, 0, MPI_COMM_WORLD);
#endif // _MPI
            ++iteration;
            if(nsolutions >= unique_flips) break;
        }
        if(nsolutions >= unique_flips) break;
    }

    if (rank < 1) {
        for(const auto& x : solutions){
            for(unsigned b=0U; b<SITESNUMBER; ++b) ostrm<<x[b]<<" ";
            ostrm<<std::endl;
        }
    }

#ifdef _MPI
    MPI_Finalize();
#endif // _MPI
    return 0;
}
