#ifndef GA_PLUS_HPP
#define GA_PLUS_HPP

#include <vector>
#include <random>
#include <cstring>
#include <algorithm>

namespace ga_plus
{

enum struct ProblemType
{
    MAXIMIZE = 0,
    MINIMIZE
};

namespace individual
{
template<typename GeneType, typename FitnessType>
class Individual
{
    public:
    GeneType gene;
    FitnessType fitness;
};
}

namespace selection
{
template<typename GeneType, typename FitnessType, typename RandomEngine>
std::vector<individual::Individual<GeneType, FitnessType>> tournament(const std::vector<individual::Individual<GeneType, FitnessType>>& population, const std::size_t select_size, const std::size_t tournament_size, const ProblemType problem_type, RandomEngine&& engine)
{
    std::vector<individual::Individual<GeneType, FitnessType>> ret(select_size);

    std::vector<std::size_t> default_temp_indices(population.size());
    for(std::size_t i = 0; i < default_temp_indices.size(); i ++)
    {
        default_temp_indices[i] = i;
    }

    // Fisher-Yates
    std::vector<std::size_t> temp_indices(population.size());

    for(std::size_t i = 0; i < select_size; i ++)
    {
        std::memcpy(&(temp_indices[0]), &(default_temp_indices[0]), sizeof(std::size_t) * population.size());  // Initialize indices
        
        std::vector<individual::Individual<GeneType, FitnessType>> individuals(tournament_size);

        for(std::size_t j = 0; j < tournament_size; j ++)
        {
            std::uniform_int_distribution<std::size_t> distribution(j, population.size()-1);
            std::size_t index = distribution(engine);

            individuals[j].gene = population[temp_indices[index]].gene;
            individuals[j].fitness = population[temp_indices[index]].fitness;

            std::swap(temp_indices[j], temp_indices[index]);
        }

        if(problem_type == ProblemType::MAXIMIZE)
        {
            // Find the individual with the highest fitness
            auto best = std::max_element(individuals.begin(), individuals.end(),
                [](const individual::Individual<GeneType, FitnessType>& a, const individual::Individual<GeneType, FitnessType>& b)
                {
                    return a.fitness < b.fitness;
                });
            
            ret[i] = *best;
        }
        else
        {
            // Find the individual with the lowest fitness
            auto best = std::min_element(individuals.begin(), individuals.end(),
                [](const individual::Individual<GeneType, FitnessType>& a, const individual::Individual<GeneType, FitnessType>& b)
                {
                    return a.fitness < b.fitness;
                });
            
            ret[i] = *best;
        }
    }

    return ret;
}
}

}

#endif /* GA_PLUS_HPP */