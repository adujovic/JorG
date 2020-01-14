#ifndef _ISING_H
#define _ISING_H

#include <iostream>

#include <tuple>
#include <array>
#include <vector>
#include <bitset>

#include <algorithm>
#include <random>
#include <string>

#include <memory>
#include <limits>
#include <type_traits>

#ifdef _OPENMP
#include <valarray>
#include <omp.h>
#endif

#include "asa.h"
#include "arithmeticvector.h"

namespace ising{

template<size_t N> double isingEnergy (void* state);
template<size_t N> double isingMeasure(void* stateI,
                                       void* stateJ);
template<size_t N> void   isingStep   (const gsl_rng* random,
                                       void* state,
                                       double step);
template<size_t N> void   isingPrint  (void* state);


template<size_t N, class Model>
struct LatticeType{
    LatticeType() = delete;
    LatticeType(const LatticeType& rhs);
    LatticeType& operator=(const LatticeType& rhs);
    LatticeType(std::shared_ptr<std::mt19937>& _randomEngine,
            std::shared_ptr<std::uniform_int_distribution
                                  <unsigned long long int>>& _uniform,
            const Model* _model,
            std::bitset<N>* _nodes = nullptr);

    std::shared_ptr<std::mt19937> randomEngine;
    std::shared_ptr<std::uniform_int_distribution
                          <unsigned long long int>> uniform;
    const Model* model;
    std::bitset<N> nodes;
    size_t flipper;
};

class AbstractIsingModel{
public:
    virtual ~AbstractIsingModel() = 0;
};

#ifdef _OPENMP
//template<size_t N, size_t numOfThreads> // N-site lattie
template<size_t N> // N-site lattie
#else
template<size_t N> // N-site lattie
#endif
class IsingModel : public AbstractIsingModel{
public:
    IsingModel(bool QUICK_START=true):
        randomDevice(),
        randomEngine(std::make_shared<std::mt19937>(randomDevice())),
        uniform(std::make_shared<
                std::uniform_int_distribution<
                unsigned long long int>>(
                    std::numeric_limits<
                    unsigned long long int>::min(),
                    std::numeric_limits<
                    unsigned long long int>::max())),
        lattice(randomEngine,uniform,this){
            if(QUICK_START) hamiltonian.reserve(N*N/2); // maximum number of interactions (each site with all others)
            else        hamiltonian.reserve(N);     // minimum number of interactions (one per site)

            solver.set_energy (  isingEnergy<N> );
            solver.set_measure( isingMeasure<N> );
            solver.set_step   (    isingStep<N> );
#ifdef _VERBOSE
            solver.set_print  (   isingPrint<N> );
#endif
#ifdef _OPENMP
//            energy     = std::valarray<double>(0.0,numOfThreads);
//            probablity = std::valarray<double>(0.0,numOfThreads);
#endif
        }

    ~IsingModel() override {}

    typedef std::tuple<unsigned,unsigned,double>  TwoSiteInteraction;
    typedef std::vector<TwoSiteInteraction>       HamiltonianType;
    typedef celerium::ArithmeticVector            VectorType;

protected:
    std::random_device                              randomDevice;
    std::shared_ptr<std::mt19937>                   randomEngine;
    std::shared_ptr<std::uniform_int_distribution
                         <unsigned long long int>>    uniform;
    LatticeType<N,IsingModel>                         lattice;

    gsl::SimulatedAnnealing                           solver;
    HamiltonianType                                   hamiltonian;

    std::array<VectorType,3>                          basis;
    std::vector<std::tuple<size_t,VectorType,double>> supercell;

    size_t                                            referencePoint;

    const std::bitset<N>*                             mask;
#ifdef _OPENMP
//    std::valarray<double>                             energy;
//    std::valarray<double>                             probablity;
#endif

public:
    typename gsl::SimulatedAnnealing::Parameters&
    set_parameters(const typename gsl::SimulatedAnnealing::Parameters& _params){
            return solver.set_parameters(_params);
    }

    void set_basis(const std::array<VectorType,3>& _basis){
        basis = _basis;
    }

    void set_supercell(const std::vector<std::tuple<size_t,VectorType,double>>& _supercell){
        supercell = _supercell;
    }

    void set_reference(size_t _reference){
        referencePoint = _reference;
    }

    static unsigned randomize(std::bitset<N>& state, 
                              std::shared_ptr<std::mt19937>& randomEngine,
                              std::shared_ptr<
                              std::uniform_int_distribution
                                   <unsigned long long int>>& uniform,
                              const std::bitset<N>* mask=nullptr){ 
        unsigned ones  = (*uniform)(*randomEngine)%static_cast<size_t>(sqrt(N));
        unsigned zeros = N - ones;
        std::string face(ones,'1');
        face += std::string(zeros,'0');
        std::shuffle(face.begin(),face.end(),*randomEngine);
        state ^= std::bitset<N>(face);

        if(mask != nullptr) state &= *mask; // if mask exist, don't flip marked spins

        return 0;
    }

    static unsigned randomize(std::bitset<N>& state, 
                              std::mt19937& randomEngine,
                              size_t maxNumberOfFlips=
                                     static_cast<size_t>(1.0+2.0*log(N)),
                              const std::bitset<N>* mask=nullptr){ 
        unsigned ones  = std::min(maxNumberOfFlips,N);
        unsigned zeros = N - ones;
        std::string face(ones,'1');
        face += std::string(zeros,'0');
        std::shuffle(face.begin(),face.end(),randomEngine);
        state ^= std::bitset<N>(face);

        if(mask != nullptr) state &= *mask; // if mask exist, don't flip marked spins

        return 0;
    }

    const std::bitset<N>* get_mask() const{
        return mask;
    }

protected:
    // For SFINAE compiler-time evaulation
    template<class T>
    T tester(T t)const{
        if(std::is_integral<T>::value) return static_cast<unsigned>(t);
        return t;
    }

public:
    void add_interaction(...){
        std::cerr<<"Wrong input for auxiliary::IsingModel::add_interaction:"<<std::endl;
        std::cerr<<"       Either non-iterable or iterable of non <int,int,float> tuples."<<std::endl;
        std::cerr<<"       Nothing was added!"<<std::endl;
    }

    template<class intlike, class floatlike>
    auto add_interaction(intlike i, intlike j, floatlike J) -> decltype((unsigned)(tester<intlike>)(i),void()){
        hamiltonian.push_back(std::make_tuple(i,j,J));
    }

    template<class T>
    auto add_interaction(T interaction) -> decltype((TwoSiteInteraction&)(tester<T>)(interaction),void()){
        hamiltonian.push_back(interaction);
    }

    template<class Iterable>
    auto add_interaction(const Iterable& interactions) -> decltype((decltype(interactions.begin()))(std::begin)(interactions),void()){
        for(auto& interaction : interactions)
            hamiltonian.push_back(interaction);
    }

    std::bitset<N>* get_nodes_ptr(){
        return &(lattice.nodes);
    }

    void clear_hamiltonian(){
        hamiltonian.clear();
    }

    void reset(){
        this->clear_hamiltonian();
        this->randomize_state();
    }

    void randomize_state(){
        IsingModel<N>::generate_state(lattice.nodes,*(lattice.randomEngine),*(lattice.uniform));
    }

    const std::bitset<N>& get_nodes() const{
        return lattice.nodes;
    }

    const HamiltonianType& get_hamiltonian() const{
        return hamiltonian;
    }

    std::bitset<N>& get_nodes(){
        return lattice.nodes;
    }

    friend std::ostream& operator<<(std::ostream& stream, const IsingModel& model){
        return stream<<model.get_nodes();
    }

    static double energy(const LatticeType<N,IsingModel>* state){
        double E = 0.0;
#ifdef _OPENMP
        #pragma omp parallel
        {
        #pragma omp for reduction(+:E)
#endif
        for(size_t i = 0U; i < state->model->get_hamiltonian().size(); ++i){
            E += std::get<2>(state->model->get_hamiltonian()[i])
                *(state->nodes[std::get<0>(state->model->get_hamiltonian()[i])]-0.5)
                *(state->nodes[std::get<1>(state->model->get_hamiltonian()[i])]-0.5);
        }
#ifdef _OPENMP
        }
#endif
        return E;
    }

    static double measure(const std::bitset<N>& stateI, const std::bitset<N>& stateJ){
        std::bitset<N> output = ~stateI & stateJ;
        return output.count();
    }    

    std::bitset<N> run(){
        this->mask = nullptr;
#ifdef _VERBOSE
        std::cout<<"Starting from: "<<lattice.nodes<<std::endl;
#endif
        solver.run<LatticeType<N,IsingModel<N>>>(lattice,sizeof(lattice));
#ifdef _VERBOSE
        std::cout<<"Solution:      ";
        std::cout<<lattice.nodes<<std::endl;
#endif
        return lattice.nodes;
    }

#ifdef _VERBOSE
    std::bitset<N> run(std::bitset<N>* mask){
        this->mask = mask;
        lattice.nodes &= *mask;
        std::cout<<"Starting from: "<<lattice.nodes<<std::endl;
        solver.run<LatticeType<N,IsingModel<N>>>(lattice,sizeof(lattice));
        std::cout<<"Solution:      ";
        std::cout<<lattice.nodes<<std::endl;
        return lattice.nodes;
    }
#else
    std::bitset<N> run(std::bitset<N>* mask){
        this->mask = mask;
        lattice.nodes &= *mask;
        solver.run<LatticeType<N,IsingModel<N>>>(lattice,sizeof(lattice));
        return lattice.nodes;
    }
#endif

    static void generate_state(std::bitset<N>& state,
                               std::mt19937& engine,
                               std::uniform_int_distribution
                                    <unsigned long long int>& distribution){
        constexpr auto seedSize = 8*sizeof(unsigned long long int);

        state = std::bitset<N>(distribution(engine));
        auto currentSize = seedSize;

        while (currentSize < N){
            state <<= seedSize;
            state  |= std::bitset<N>(distribution(engine));
            currentSize += seedSize;
        }
    }
}; // end of class IsingModel

template<size_t N>
double isingEnergy (void* state){
    return IsingModel<N>::energy(static_cast<LatticeType<N,IsingModel<N>>*>(state));
}

template<size_t N>
double isingMeasure(void* stateI, void* stateJ){
    return IsingModel<N>::measure(
            static_cast<LatticeType<N,IsingModel<N>>*>(stateI)->nodes,
            static_cast<LatticeType<N,IsingModel<N>>*>(stateJ)->nodes
           );
}

template<size_t N>
void   isingStep   (const gsl_rng* random __attribute__((unused)), void* state, double step __attribute__((unused))){
    IsingModel<N>::randomize(
            static_cast<LatticeType<N,IsingModel<N>>*>(state)->nodes,
            static_cast<LatticeType<N,IsingModel<N>>*>(state)->randomEngine,
            static_cast<LatticeType<N,IsingModel<N>>*>(state)->uniform,
            static_cast<LatticeType<N,IsingModel<N>>*>(state)->model->get_mask()
           );
}

template<size_t N>
void   isingPrint  (void* state){
#ifndef _QUIET
    std::cout<<'\t'<<static_cast<LatticeType<N,IsingModel<N>>*>(state)->nodes;
#endif
}

template<size_t N, class Model>
LatticeType<N,Model>::LatticeType(const LatticeType<N,Model>& rhs){
    randomEngine = rhs.randomEngine;
    uniform      = rhs.uniform;
    nodes        = rhs.nodes;
    model        = rhs.model;
}

template<size_t N, class Model>
LatticeType<N,Model>& LatticeType<N,Model>::operator=(const LatticeType<N,Model>& rhs){
    randomEngine = rhs.randomEngine;
    uniform      = rhs.uniform;
    nodes        = rhs.nodes;
    model        = rhs.model;
}

template<size_t N, class Model>
LatticeType<N,Model>::LatticeType(
        std::shared_ptr<std::mt19937>& _randomEngine,
        std::shared_ptr<
          std::uniform_int_distribution
               <unsigned long long int>>& _uniform,
        const Model* _model,
        std::bitset<N>* _nodes):randomEngine(_randomEngine),
                             uniform(_uniform),
                             model(_model){
        if (_nodes == nullptr)                                 
          IsingModel<N>::generate_state(nodes,*_randomEngine,*_uniform);    
        else
          nodes = *_nodes;
    }

} //end of namespace ising

#endif
