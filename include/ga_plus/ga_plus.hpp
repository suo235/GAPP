#ifndef GA_PLUS_HPP
#define GA_PLUS_HPP

namespace ga_plus
{

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

}

#endif /* GA_PLUS_HPP */