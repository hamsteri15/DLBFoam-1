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

\*---------------------------------------------------------------------------*/

#include "LoadBalancedChemistryModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template <class ReactionThermo, class ThermoType>
Foam::LoadBalancedChemistryModel<ReactionThermo, ThermoType>::
    LoadBalancedChemistryModel(const ReactionThermo& thermo)
    : 
        StandardChemistryModel<ReactionThermo, ThermoType>(thermo),
        balancer_(createBalancer()), 
        cpuTimes_
        (
            IOobject
            (
                thermo.phasePropertyName("cellCpuTimes"),
                this->time().timeName(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            this->mesh(),
            scalar(0.0)
        )
    {
        if(balancer_.log())
        {
            cpuSolveFile_ = logFile("cpu_solve.out");
            cpuSolveFile_() << "                  time" << tab
                            << "           getProblems" << tab  
                            << "           updateState" << tab
                            << "               balance" << tab
                            << "           solveBuffer" << tab
                            << "             unbalance" << tab
                            << "               rank ID" << endl;
        }

    }

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template <class ReactionThermo, class ThermoType>
Foam::LoadBalancedChemistryModel<ReactionThermo, ThermoType>::
    ~LoadBalancedChemistryModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //



template <class ReactionThermo, class ThermoType>
Foam::LoadBalancer
Foam::LoadBalancedChemistryModel<ReactionThermo, ThermoType>::createBalancer()
{
    const IOdictionary chemistryDict_tmp
        (
            IOobject
            (
                this->thermo().phasePropertyName("chemistryProperties"),
                this->thermo().db().time().constant(),
                this->thermo().db(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        );

    return LoadBalancer(chemistryDict_tmp);
}


template <class ReactionThermo, class ThermoType>
template <class DeltaTType>
Foam::scalar Foam::LoadBalancedChemistryModel<ReactionThermo, ThermoType>::solve
(
    const DeltaTType& deltaT
)
{
    // CPU time analysis
    clockTime timer;
    scalar t_getProblems(0);
    scalar t_updateState(0);
    scalar t_balance(0);
    scalar t_solveBuffer(0);
    scalar t_unbalance(0);

    BasicChemistryModel<ReactionThermo>::correct();

    if(!this->chemistry_)
    {
        return great;
    }

    timer.timeIncrement();
    DynamicList<ChemistryProblem> allProblems = getProblems(deltaT);
    t_getProblems = timer.timeIncrement();

    RecvBuffer<ChemistrySolution> incomingSolutions;

    if(balancer_.active())
    {
        timer.timeIncrement();
        balancer_.updateState(allProblems);
        t_updateState = timer.timeIncrement();

        timer.timeIncrement();
        auto guestProblems = balancer_.balance(allProblems);
        auto ownProblems = balancer_.getRemaining(allProblems);
        t_balance = timer.timeIncrement();

        timer.timeIncrement();
        auto ownSolutions = solveList(ownProblems);
        auto guestSolutions = solveBuffer(guestProblems);
        t_solveBuffer = timer.timeIncrement();

        timer.timeIncrement();      
        incomingSolutions = balancer_.unbalance(guestSolutions);
        incomingSolutions.append(ownSolutions);
        t_unbalance = timer.timeIncrement();
    }
    else
    {
        timer.timeIncrement();
        incomingSolutions.append(solveList(allProblems));
        t_solveBuffer = timer.timeIncrement();
    }
        
    if(balancer_.log())
    {
        if(balancer_.active())
        {
            balancer_.printState();
        }
        cpuSolveFile_() << setw(22)
                        << this->time().timeOutputValue()<<tab
                        << setw(22) << t_getProblems<<tab
                        << setw(22) << t_updateState<<tab
                        << setw(22) << t_balance<<tab
                        << setw(22) << t_solveBuffer<<tab
                        << setw(22) << t_unbalance<<tab
                        << setw(22) << Pstream::myProcNo()
                        << endl;
    }

    return updateReactionRates(incomingSolutions);
}



template <class ReactionThermo, class ThermoType>
Foam::ChemistrySolution Foam::LoadBalancedChemistryModel<ReactionThermo, ThermoType>::solveSingle
(
    const ChemistryProblem& problem
) const
{
    ChemistrySolution solution(this->nSpecie());

    scalar timeLeft = problem.deltaT;

    // Timer begins
    clockTime time;
    time.timeIncrement();

    // Define a const label to pass as the cell index placeholder
    const label arbitrary = 0;

    scalar p = problem.pi;
    scalar T = problem.Ti;
    scalarField c = problem.c;
    scalarField c0 = problem.c;
    scalar deltaTChem = problem.deltaTChem;

    // Calculate the chemical source terms
    while(timeLeft > small)
    {
        scalar dt = timeLeft;

        this->solve(
            p,
            T,
            c,
            arbitrary,
            dt,
            deltaTChem);


        timeLeft -= dt;
    }

    solution.c_increment = (c - c0) / problem.deltaT;
    solution.deltaTChem = min(deltaTChem, this->deltaTChemMax_);

    // Timer ends
    solution.cpuTime = time.timeIncrement();

    solution.cellid = problem.cellid;
    solution.rhoi = problem.rhoi;



    return solution;
}

template <class ReactionThermo, class ThermoType>
Foam::scalar
Foam::LoadBalancedChemistryModel<ReactionThermo, ThermoType>::updateReactionRates
(
    const RecvBuffer<ChemistrySolution>& solutions
)
{
    scalar deltaTMin = great;

    auto policy = std::execution::par_unseq;   

    for(const auto& array : solutions)
    {
        std::for_each(policy, array.begin(), array.end(), [&](const ChemistrySolution& s){this->updateResult(s);});
        
        /*
        const auto& solution = *std::min_element(
            array.begin(),
            array.end(),
            [](const ChemistrySolution& lhs, const ChemistrySolution& rhs) {return lhs.deltaTChem < rhs.deltaTChem;}
            );

        deltaTMin = min(solution.deltaTChem, deltaTMin);
        */
        //TODO: replace this with std::min_element using a parallel policy
        for(const auto& solution : array)
        {
            deltaTMin = min(solution.deltaTChem, deltaTMin);
        }
    }

    return deltaTMin;
}


template <class ReactionThermo, class ThermoType>
Foam::scalar Foam::LoadBalancedChemistryModel<ReactionThermo, ThermoType>::solve
(
    const scalarField& deltaT
)
{
    return this->solve<scalarField>(deltaT);
}


template <class ReactionThermo, class ThermoType>
Foam::scalar Foam::LoadBalancedChemistryModel<ReactionThermo, ThermoType>::solve
(
    const scalar deltaT
)
{
    // Don't allow the time-step to change more than a factor of 2
    return min(
        this->solve<UniformField<scalar>>(UniformField<scalar>(deltaT)),
        2 * deltaT);
}


template <class ReactionThermo, class ThermoType>
Foam::RecvBuffer<Foam::ChemistrySolution>
Foam::LoadBalancedChemistryModel<ReactionThermo, ThermoType>::solveBuffer
(
    RecvBuffer<ChemistryProblem>& problems
) const
{
    // allocate the solutions buffer
    RecvBuffer<ChemistrySolution> solutions;
    
    for(auto& p : problems)
    {

        solutions.append(solveList(p));
    }
    return solutions;
}


template <class ReactionThermo, class ThermoType>
Foam::DynamicList<Foam::ChemistrySolution>
Foam::LoadBalancedChemistryModel<ReactionThermo, ThermoType>::solveList
(
    const UList<ChemistryProblem>& problems
) const
{
    DynamicList<ChemistrySolution> solutions(
        problems.size(), ChemistrySolution(this->nSpecie_));

    auto policy = std::execution::par_unseq;


    std::transform(policy, problems.begin(), problems.end(), solutions.begin(), 
        [&](const ChemistryProblem& p) {return this->solveSingle(p);}
    );



    return solutions;
}



template <class ReactionThermo, class ThermoType>
template<class DeltaTType>
Foam::DynamicList<Foam::ChemistryProblem>
Foam::LoadBalancedChemistryModel<ReactionThermo, ThermoType>::getProblems
(
    const DeltaTType& deltaT
)
{
    const scalarField& T = this->thermo().T();
    const scalarField& p = this->thermo().p();
    tmp<volScalarField> trho(this->thermo().rho());
    const scalarField& rho = trho();
    


    DynamicList<ChemistryProblem> solved_problems;

    solved_problems.resize(p.size(), ChemistryProblem(this->nSpecie_));

    scalarField concentration(this->nSpecie_);

    label counter = 0;
    forAll(T, celli)
    {

        if(T[celli] > this->Treact())
        {
            for(label i = 0; i < this->nSpecie_; i++)
            {
                concentration[i] = rho[celli] * this->Y_[i][celli] / this->specieThermos_[i].W();
            }
            
            ChemistryProblem problem;
            problem.c = concentration;     
            problem.Ti = T[celli];
            problem.pi = p[celli];
            problem.rhoi = rho[celli];
            problem.deltaTChem = this->deltaTChem_[celli];
            problem.deltaT = deltaT[celli];
            problem.cpuTime = cpuTimes_[celli];
            problem.cellid = celli;


            solved_problems[counter] = problem;
            counter++;

        }
        else
        {
            for(label i = 0; i < this->nSpecie(); i++)
            {
                this->RR_[i][celli] = 0;
            }
        }

    }

    //the real size is set here
    solved_problems.setSize(counter);

    

    return solved_problems;
}



template <class ReactionThermo, class ThermoType>
void Foam::LoadBalancedChemistryModel<ReactionThermo, ThermoType>::updateResult
(
    const ChemistrySolution& solution
)
{
    const label cellid = solution.cellid;

    for(label j = 0; j < this->nSpecie_; j++)
    {
        this->RR_[j][cellid] = solution.c_increment[j] * this->specieThermos_[j].W();
    }
    this->deltaTChem_[cellid] = min(solution.deltaTChem, this->deltaTChemMax_);

    cpuTimes_[solution.cellid] = solution.cpuTime;




}



