/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | DLBFoam: Dynamic Load Balancing 
   \\    /   O peration     | for fast reactive simulations
    \\  /    A nd           | 
     \\/     M anipulation  | 2020, Aalto University, Finland
-------------------------------------------------------------------------------
License
    This file is part of DLBFoam library, derived from OpenFOAM.
    https://github.com/blttkgl/DLBFoam
    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.
    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.
    
\*---------------------------------------------------------------------------*/

#include "noChemistrySolver.H"
#include "EulerImplicit.H"
#include "ode.H"
#include "LoadBalancedChemistryModel.H"

#include "psiReactionThermo.H"
#include "rhoReactionThermo.H"

#include "forCommonGases.H"
#include "forCommonLiquids.H"
#include "forPolynomials.H"
#include "makeChemistrySolver.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define defineChemistrySolversA(ReactionThermo, ThermoPhysics)                  \
    defineChemistrySolver                                                      \
    (                                                                          \
        LoadBalancedChemistryModel,                                            \
        ReactionThermo,                                                        \
        ThermoPhysics                                                          \
    );                                                                         \

#define makeChemistrySolversA(Solver, ReactionThermo, ThermoPhysics)            \
    makeChemistrySolver                                                        \
    (                                                                          \
        Solver,                                                                \
        LoadBalancedChemistryModel,                                            \
        ReactionThermo,                                                        \
        ThermoPhysics                                                          \
    );                                                                         \

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    //using standard ode

    //gasses
    forCommonGases(defineChemistrySolversA, psiReactionThermo);
    forCommonGases(defineChemistrySolversA, rhoReactionThermo);

    //forCommonGases(makeChemistrySolversA, noChemistrySolver, psiReactionThermo);
    //forCommonGases(makeChemistrySolversA, noChemistrySolver, rhoReactionThermo);
    //forCommonGases(makeChemistrySolversA, EulerImplicit, psiReactionThermo);
    //forCommonGases(makeChemistrySolversA, EulerImplicit, rhoReactionThermo);
    forCommonGases(makeChemistrySolversA, ode, psiReactionThermo);
    forCommonGases(makeChemistrySolversA, ode, rhoReactionThermo);

    /*
    //liquids
    forCommonLiquids(defineChemistrySolversA, rhoReactionThermo);

    forCommonLiquids(makeChemistrySolversA, noChemistrySolver, rhoReactionThermo);
    forCommonLiquids(makeChemistrySolversA, EulerImplicit, rhoReactionThermo);
    forCommonLiquids(makeChemistrySolversA, ode, rhoReactionThermo);

    //polynomials
    forPolynomials(defineChemistrySolversA, rhoReactionThermo);

    forPolynomials(makeChemistrySolversA, noChemistrySolver, rhoReactionThermo);
    forPolynomials(makeChemistrySolversA, EulerImplicit, rhoReactionThermo);
    forPolynomials(makeChemistrySolversA, ode, rhoReactionThermo);
    */
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //