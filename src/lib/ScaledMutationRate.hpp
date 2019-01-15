#pragma once

/**
 * \brief An interface for the scaled mutation rate (theta=4*Ne*u) at a given taxon of tree, where
 * Ne is the population size and u is the mutation rate.
 *
 * This abstract class is meant as a general read-only interface returning theta=4*Ne*u.
 *
 * A good example is the class PolyProcess. All calculations performed by PolyProcess do not depend
 * on the specific implementation (and calculation) of theta=4*Ne*u.
 * As a consequence, the only thing required by PolyProcess is a generic interface  specifying which
 * theta should be used for a given taxon.
 *
 * In practice, this interface will hide two distinct implementations: theta=4*Ne*u can be either a
 * free parameter or a compound parameter.
 * In the case theta is a free parameter (homogeneous per
 * taxon or heterogeneous), Ne and u are not distinguishable parameters.
 * In the case theta is a compound parameter (homogeneous per
 * taxon or heterogeneous), Ne and u are distinguishable parameters and
 * theta=4*Ne*u is computed as a product of both (see CompoundScaledMutationRate in
 * NodeMultivariateProcess.hpp).
 */
class ScaledMutationRate {
  public:
    virtual ~ScaledMutationRate() = default;

    //! return theta for a given taxon
    virtual double GetTheta(int taxon) const = 0;
};

/**
 * \brief A ScaledMutationRate that returns the scaled mutation rate (theta=4*Ne*u) for any taxon.
 *
 * Useful for implementing models assuming the same Ne (population size) and u (mutation rate
 * per generation) all taxa. With these assumptions theta=4*Ne*u is a free parameter and not a
 * compound parameter.
 */
class HomogeneousScaledMutationRate : public ScaledMutationRate {
  public:
    //! \brief Constructor, taking as its arguments the value to be
    //! returned for any taxon
    explicit HomogeneousScaledMutationRate(double const &intheta_scale)
        : theta_scale{intheta_scale} {};

    ~HomogeneousScaledMutationRate() override = default;

    double GetTheta() const { return theta_scale; }
    double GetTheta(int taxon) const override { return theta_scale; }

  private:
    double const &theta_scale;
};
