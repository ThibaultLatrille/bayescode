#pragma once

#include <assert.h>
#include <set>
#include <string>
#include "lib/Array.hpp"
#include "lib/BidimArray.hpp"
#include "lib/BranchArray.hpp"
#include "lib/NodeArray.hpp"

class PoissonSuffStat;
class PoissonSuffStatBranchArray;
class IIDGamma;
class BranchIIDGamma;
class IIDDirichlet;
class MultiDirichlet;
class MultinomialAllocationVector;
class GammaWhiteNoise;
class GammaWhiteNoiseArray;
class Dirichlet;
class BidimIIDMultiGamma;
class IIDProfileMask;
class NodeAges;

/*
====================================================================================================
  Debug helpers
==================================================================================================*/
// TODO move into a logger class somewhere
#ifdef __GNUG__
#include <cxxabi.h>
static inline std::string demangle(const char* name) {
    int status{0};
    std::unique_ptr<char, void (*)(void*)> res{
        abi::__cxa_demangle(name, NULL, NULL, &status), std::free};
    return (status == 0) ? res.get() : name;
}
#else
static std::string demangle(const char* name) { return name; }
#endif

/*
====================================================================================================
  RegistrarBase class
==================================================================================================*/
using filter_t = std::set<std::string>;

#define CONVERT_REF_TO(FROM, TO)                                         \
    template <class... Args>                                             \
    void register_element(std::string s, FROM& target, Args&&... args) { \
        convert_ref<TO>(s, target, std::forward<Args>(args)...);         \
    }

template <class T>
class RegistrarBase {
    mutable filter_t _filter;

    template <class TargetType, class OriginalType, class... Args>
    void convert_ref(std::string s, OriginalType& target, Args&&... args) {
        try {
            auto& converted_ref = dynamic_cast<TargetType&>(target);
            static_cast<T*>(this)->register_element(s, converted_ref, std::forward<Args>(args)...);
        } catch (const std::bad_cast& e) {
            std::cerr << "\e[31mError\e[0m| Failed to convert from type \e[32m"
                      << demangle(typeid(OriginalType).name()) << "\e[0m to type \e[32m"
                      << demangle(typeid(TargetType).name()) << "\e[0m\n";
            exit(1);
        }
    }

  public:
    CONVERT_REF_TO(IIDGamma, SimpleArray<double>);
    CONVERT_REF_TO(Dirichlet, SimpleArray<double>);
    CONVERT_REF_TO(BranchIIDGamma, SimpleBranchArray<double>);
    CONVERT_REF_TO(IIDDirichlet, SimpleArray<std::vector<double>>);
    CONVERT_REF_TO(MultiDirichlet, SimpleArray<std::vector<double>>);
    CONVERT_REF_TO(MultinomialAllocationVector, SimpleArray<int>);
    CONVERT_REF_TO(GammaWhiteNoise, SimpleBranchArray<double>);
    CONVERT_REF_TO(NodeAges, SimpleNodeArray<double>);
    CONVERT_REF_TO(GammaWhiteNoiseArray, Array<GammaWhiteNoise>);
    CONVERT_REF_TO(PoissonSuffStatBranchArray, SimpleBranchArray<PoissonSuffStat>);
    CONVERT_REF_TO(BidimIIDMultiGamma, SimpleBidimArray<std::vector<double>>);
    CONVERT_REF_TO(IIDProfileMask, SimpleArray<std::vector<int>>);

    template <class Elem, class... Args>
    void register_element(std::string s, std::vector<std::vector<Elem>>& vv, Args&&... args) {
        for (size_t i = 0; i < vv.size(); i++) {
            static_cast<T*>(this)->register_element(
                s + "_" + std::to_string(i), vv.at(i), std::forward<Args>(args)...);
        }
    }

    template <class Target, class... Args>
    void register_element(std::string s, SimpleBranchArray<Target>& target, Args&&... args) {
        static_cast<T*>(this)->register_element(s, target.GetArray(), std::forward<Args>(args)...);
    }

    template <class Target, class... Args>
    void register_element(std::string s, SimpleNodeArray<Target>& target, Args&&... args) {
        static_cast<T*>(this)->register_element(s, target.GetArray(), std::forward<Args>(args)...);
    }

    template <class Target, class... Args>
    void register_element(std::string s, SimpleArray<Target>& target, Args&&... args) {
        static_cast<T*>(this)->register_element(s, target.GetArray(), std::forward<Args>(args)...);
    }

    template <class Target, class... Args>
    void register_element(std::string s, SimpleBidimArray<Target>& target, Args&&... args) {
        static_cast<T*>(this)->register_element(s, target.GetArray(), std::forward<Args>(args)...);
    }

    template <class Target, class... Args>
    void register_element(std::string s, Target&, Args... args) {
        std::stringstream ss;
        ss << "\e[1m\e[31mError\e[0m| \e[33mregister_element\e[0m overload for element \e[32m" << s
           << "\e[0m not implemented for class \e[32m" << demangle(typeid(T).name())
           << "\e[0m\n     | reference type is \e[32m" << demangle(typeid(Target).name())
           << "\e[0m\n     \\ other parameters are: \e[32m";
        std::vector<std::string> types{demangle(typeid(Args).name())...};
        for (auto t : types) { ss << t << " "; }
        ss << "\e[0m\n";
        std::cerr << ss.str();
        exit(1);
    }

    RegistrarBase() = default;

    template <class... Args>
    void add(std::string name, Args&&... args) {
        if ((_filter.size() == 0) or (_filter.find(name) != _filter.end())) {
            static_cast<T*>(this)->register_element(name, std::forward<Args>(args)...);
        }
    }

    template <class M>
    void register_from_method(
        M& ref, void (M::*f)(RegistrarBase<T>&), filter_t filter = filter_t{}) {
        /* -- */
        assert(filter.size() > 0);
        _filter = filter;  // local copy
        (ref.*f)(*this);
        _filter.clear();
    }
};
