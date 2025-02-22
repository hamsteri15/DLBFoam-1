/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | DLBFoam: Dynamic Load Balancing 
   \\    /   O peration     | for fast reactive simulations
    \\  /    A nd           | 
     \\/     M anipulation  | 2020, Aalto University, Finland
-------------------------------------------------------------------------------
License
    This file is part of DLBFoam library, derived from OpenFOAM.

    https://github.com/Aalto-CFD/DLBFoam

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

Class
    Foam::LoadBalancedChemistryModel

Description
    Extends standardChemistryModel model by adding a load balancing algorithm.

SourceFiles
    LoadBalancedChemistryModel.C

\*---------------------------------------------------------------------------*/

#ifndef LoadBalancedChemistryModel_H
#define LoadBalancedChemistryModel_H

#include "ChemistryProblem.H"
#include "ChemistrySolution.H"
#include "LoadBalancer.H"
#include "OFstream.H"
#include "IOmanip.H"
#include "StandardChemistryModel.H"
#include "clockTime.H"

#include <execution> //execution policies

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class LoadBalancedChemistryModel Declaration
\*---------------------------------------------------------------------------*/

template<class ReactionThermo, class ThermoType>
class LoadBalancedChemistryModel
: 
    public StandardChemistryModel<ReactionThermo, ThermoType>
{

private:
    // Private member data

        // Load balancing object
        LoadBalancer balancer_;


        // Field containing chemistry CPU time information    
        volScalarField cpuTimes_;

        // A file to output the balancing stats
        autoPtr<OFstream>        cpuSolveFile_;


    // Private Member Functions


        //- Create a load balancer object
        LoadBalancer createBalancer();

        //- Get the list of problems to be solved
        template<class DeltaTType>
        DynamicList<ChemistryProblem> getProblems(const DeltaTType& deltaT);

        //- Solve a list of chemistry problems and return a list of solutions
        DynamicList<ChemistrySolution> 
        solveList(const UList<ChemistryProblem>& problems) const;

        //- Solve the problem buffer coming from the balancer
        RecvBuffer<ChemistrySolution> 
        solveBuffer(RecvBuffer<ChemistryProblem>& problems) const;

        //- Update the reaction rate based on the solution. Assumes that the solution.cellid is on the correct process
        void updateResult(const ChemistrySolution& solution);

        //- Update the reaction rates from a list of solutions
        scalar updateReactionRates(const RecvBuffer<ChemistrySolution>& solutions);

        //- Solve the reaction system for the given time step
        //  of given type and return the characteristic time
        template<class DeltaTType>
        scalar solve(const DeltaTType& deltaT);




public:

    //- Runtime type information
    TypeName("loadBalanced");


    // Constructors

        //- Construct from thermo
        LoadBalancedChemistryModel(const ReactionThermo& thermo);

        //- Disallow default bitwise copy construction
        LoadBalancedChemistryModel(const LoadBalancedChemistryModel&);


    //- Destructor
    virtual ~LoadBalancedChemistryModel();


    // Member Functions

        //- Solve a single ChemistryProblem and put the solution to ChemistrySolution 
        ChemistrySolution solveSingle(const ChemistryProblem& problem) const;
        
        
        // Chemistry model functions (overriding functions in
        // StandardChemistryModel to use the private solve function)

            //- Solve the reaction system for the given time step
            //  and return the characteristic time
            virtual scalar solve(const scalar deltaT) override;

            //- Solve the reaction system for the given time step
            //  and return the characteristic time
            virtual scalar solve(const scalarField& deltaT) override;

        // ODE functions (overriding abstract functions in ODE.H)
        virtual void solve
        (
            scalar&      p,
            scalar&      T,
            scalarField& c,
            const label  li,
            scalar&      deltaT,
            scalar&      subDeltaT
        ) const = 0;

        //- Create and return a log file of the given name
        inline autoPtr<OFstream> logFile(const word& name) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const LoadBalancedChemistryModel&) = delete;


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam
#include "LoadBalancedChemistryModelI.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "LoadBalancedChemistryModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
