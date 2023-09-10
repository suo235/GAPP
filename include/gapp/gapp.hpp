#ifndef GA_PLUS_HPP
#define GA_PLUS_HPP

#include <vector>
#include <random>
#include <cstring>
#include <algorithm>

namespace gapp
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
individual::Individual<GeneType, FitnessType> tournament(const std::vector<individual::Individual<GeneType, FitnessType>>& population, const std::size_t tournament_size, const ProblemType problem_type, RandomEngine&& engine)
{
    using Individual = individual::Individual<GeneType, FitnessType>;

    std::vector<Individual*> pointers(population.size());

    for(std::size_t i = 0; i < pointers.size(); i ++)
    {
        pointers[i] = const_cast<Individual*>(population.data()) + i;
    }

    // Fisher-Yates
    std::vector<std::size_t> temp_indices(population.size());

    for(std::size_t i = 0; i < tournament_size; i ++)
    {
        std::uniform_int_distribution<std::size_t> distribution(i, pointers.size()-1);
        std::size_t index = distribution(engine);
        
        std::swap(pointers[i], pointers[index]);
    }

    if(problem_type == ProblemType::MAXIMIZE)
    {
        // Find the individual with the highest fitness
        auto best = std::max_element(pointers.begin(), pointers.end(),
            [](const Individual* a, const Individual* b)
            {
                return a->fitness < b->fitness;
            });
        
        return *(*best);
    }
    else
    {
        // Find the individual with the lowest fitness
        auto best = std::min_element(pointers.begin(), pointers.end(),
            [](const Individual* a, const Individual* b)
            {
                return a->fitness < b->fitness;
            });
        
        return *(*best);
    }
}

template<typename GeneType, typename FitnessType, typename RandomEngine>
std::vector<individual::Individual<GeneType, FitnessType>> select_random_k(const std::vector<individual::Individual<GeneType, FitnessType>>& population, const std::size_t select_size, RandomEngine&& engine)
{
    std::vector<individual::Individual<GeneType, FitnessType>> ret(select_size);

    std::vector<std::size_t> temp_indices(population.size());
    for(std::size_t i = 0; i < temp_indices.size(); i ++)
    {
        temp_indices[i] = i;
    }

    // Fisher-Yates
    for(std::size_t i = 0; i < select_size; i ++)
    {
        std::uniform_int_distribution<std::size_t> distribution(i, population.size()-1);
        std::size_t index = distribution(engine);
        
        ret[i].gene = population[temp_indices[index]].gene;
        ret[i].fitness = population[temp_indices[index]].fitness;

        std::swap(temp_indices[i], temp_indices[index]);
    }

    return ret;
}

template<typename GeneType, typename FitnessType>
std::vector<individual::Individual<GeneType, FitnessType>> select_best_k(const std::vector<individual::Individual<GeneType, FitnessType>>& population, const std::size_t select_size, const ProblemType problem_type)
{
    using Individual = individual::Individual<GeneType, FitnessType>;

    std::vector<Individual> ret(select_size);

    std::vector<Individual*> pointers(population.size());

    // ポインタ配列を作成
    for(std::size_t i = 0; i < pointers.size(); i ++)
    {
        pointers[i] = const_cast<Individual*>(population.data()) + i;
    }

    if(problem_type == ProblemType::MAXIMIZE)
    {
        // fitnessが大きい順に上位select_size個のポインタを取得
        std::nth_element(pointers.begin(), pointers.begin()+select_size, pointers.end(),
            [](const Individual* a, const Individual* b)
            {
                return a->fitness > b->fitness;
            });
    }
    else
    {
        std::nth_element(pointers.begin(), pointers.begin()+select_size, pointers.end(),
            [](const Individual* a, const Individual* b)
            {
                return a->fitness < b->fitness;
            });
    }

    for(std::size_t i = 0; i < select_size; i ++)
    {
        ret[i] = *(pointers[i]);
    }

    return ret;
}
}

namespace crossover
{
template<typename GeneType, typename FitnessType, typename Engine>
std::vector<individual::Individual<GeneType, FitnessType>> k_points(const individual::Individual<GeneType, FitnessType>& a, const individual::Individual<GeneType, FitnessType>& b, const std::size_t k_points, Engine&& engine)
{
    std::vector<individual::Individual<GeneType, FitnessType>> ret{a, b};

    std::vector<std::size_t> indices(a.gene.size());
    for(std::size_t i = 0; i < indices.size(); i ++)
    {
        indices[i] = i;
    }

    // ランダムにk個のインデックスを選ぶ
    // Fisher-Yates
    for(std::size_t i = 0; i < k_points; i ++)
    {
        std::uniform_int_distribution<std::size_t> distribution(i, a.gene.size()-1);
        std::size_t index = distribution(engine);
        std::swap(indices[i], indices[index]);
    }

    if(k_points == 1)
    {
        std::swap_ranges(ret[0].gene.begin(), ret[0].gene.begin()+indices[0], ret[1].gene.begin());
        return ret;
    }

    // k個(k>=2)のインデックスを昇順ソート
    std::sort(indices.begin(), indices.begin()+k_points);

    // k_points=2: i=0
    // k_points=3: i=0
    // k_points=4: i=0, i=2
    for(std::size_t i = 0; i < k_points-1; i = i+2)
    {
        std::swap_ranges(ret[0].gene.begin()+indices[i], ret[0].gene.begin()+indices[i+1], ret[1].gene.begin()+indices[i]);
    }

    if((k_points%2) != 0)
    {
        std::swap_ranges(ret[0].gene.begin()+indices.back(), ret[0].gene.end(), ret[1].gene.begin()+indices.back());
    }

    return ret;
}

template<typename GeneType, typename FitnessType, typename Engine>
std::vector<individual::Individual<GeneType, FitnessType>> blx_alpha(const individual::Individual<GeneType, FitnessType>& a, const individual::Individual<GeneType, FitnessType>& b, const double alpha, Engine&& engine)
{
    std::vector<individual::Individual<GeneType, FitnessType>> ret{a, b};

    std::uniform_real_distribution<> generator(0.0, 1.0);

    for(std::size_t i = 0; i < a.gene.size(); i ++)
    {
        auto x_min = std::min(a.gene[i], b.gene[i]);
        auto x_max = std::max(a.gene[i], b.gene[i]);

        auto d = std::abs(a.gene[i] - b.gene[i]);

        auto r_min = x_min - alpha * d;
        auto r_max = x_max + alpha * d;

        ret[0].gene[i] = (r_max - r_min) * generator(engine) + r_min;
        ret[1].gene[i] = (r_max - r_min) * generator(engine) + r_min;
    }

    return ret;
}
}
}

#endif /* GA_PLUS_HPP */