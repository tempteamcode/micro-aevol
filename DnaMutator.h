// ***************************************************************************************************************
//
//          Mini-Aevol is a reduced version of Aevol -- An in silico experimental evolution platform
//
// ***************************************************************************************************************
//
// Copyright: See the AUTHORS file provided with the package or <https://gitlab.inria.fr/rouzaudc/mini-aevol>
// Web: https://gitlab.inria.fr/rouzaudc/mini-aevol
// E-mail: See <jonathan.rouzaud-cornabas@inria.fr>
// Original Authors : Jonathan Rouzaud-Cornabas
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// ***************************************************************************************************************


#ifndef RAEVOL_CUDA_DNAMUTATOR_H
#define RAEVOL_CUDA_DNAMUTATOR_H

#include <vector>

#include "Threefry.h"
class Organism;

/**
 * Class that generates the mutation events for a given Organism
 */
class DnaMutator {
public:

/**
 * Constructor for DnaMutator class
 *
 * Generate mutations of the DNA of an Organism
 *
 * @param mut_prng : PRNG to simulate the mutation
 * @param length  : Size of the DNA at the initialization
 * @param mutation_rate : Mutation rate of the organisms
 */
DnaMutator(Threefry::Gen&& mut_prng, int length, double mutation_rate, int indiv_id)
: mut_prng_(std::move(mut_prng))
, length_(length)
, mutation_rate_(mutation_rate)
{
}

    void generate_mutations();

    void apply_mutations(Organism& organism);

    inline bool hasMutate() const { return (!mutation_list_.empty()); }

private:
    void generate_next_mutation(int length);

    std::vector<int> mutation_list_;

    Threefry::Gen mut_prng_;
    int length_;

    double mutation_rate_;

    //--------------------------- Mutation counters
    int nb_swi_;
    int nb_mut_;
};


#endif //RAEVOL_CUDA_DNAMUTATOR_H
