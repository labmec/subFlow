//
// Created by Giovane Avancini on 01/09/2025
//

#pragma once

#include "Material/DarcyFlow/TPZMixedDarcyFlow.h"

class TSFMixedDarcy : public TPZMixedDarcyFlow {

    using TBase = TPZMixedDarcyFlow;

public:
    /**
     * @brief Default constructor
     */
    TSFMixedDarcy();

    /**
	 * @brief Class constructor
	 * @param [in] id material id
	 * @param [in] dim problem dimension
	 */
    [[maybe_unused]] TSFMixedDarcy(int id, int dim);

    /**
             copy constructor
     */
    TSFMixedDarcy(const TSFMixedDarcy &copy);
    /**
             copy constructor
     */
    TSFMixedDarcy &operator=(const TSFMixedDarcy &copy);
    /**
	 * @brief Returns a 'std::string' with the name of the material
	 */
    [[nodiscard]] std::string Name() const override { return "TSFMixedDarcy"; }

    /**
     * @brief It computes a contribution to the stiffness matrix and load vector at one integration point
     * @param[in] datavec stores all input data
     * @param[in] weight is the weight of the integration rule
     * @param[out] ek is the element matrix
     * @param[out] ef is the rhs vector
     */
    void Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek,
                    TPZFMatrix<STATE> &ef) override;

    /**
     * @brief It computes a contribution to the stiffness matrix and load vector at one BC integration point
     * @param[in] datavec stores all input data
     * @param[in] weight is the weight of the integration rule
     * @param[out] ek is the element matrix
     * @param[out] ef is the rhs vector
     * @param[in] bc is the boundary condition material
     */
    void ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek,
                      TPZFMatrix<STATE> &ef, TPZBndCondT<STATE> &bc) override;

    /*
     * @brief Fill requirements for volumetric contribute
     */
    void FillDataRequirements(TPZVec<TPZMaterialDataT<STATE> > &datavec) const override;

    /*
     * @brief Fill requirements for boundary contribute
     */
    void FillBoundaryConditionDataRequirements(int type, TPZVec<TPZMaterialDataT<STATE> > &datavec) const override;

};
