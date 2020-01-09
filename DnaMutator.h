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
#include "MutationEvent.h"


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
 * @param indiv_id : Unique identification number for the Organism
 */
DnaMutator(Threefry::Gen&& mut_prng, int length, double mutation_rate, int indiv_id)
: mut_prng_(std::move(mut_prng))
, length_(length)
, mutation_rate_(mutation_rate)
, id_(indiv_id)
{
}

    void generate_mutations();

    void generate_next_mutation(int length);

    bool mutation_available() const { return (cpt_mut_ > 0); }

    std::vector<MutationEvent> mutation_list_;

    inline bool hasMutate() const { return hasMutate_; }

    void setMutate(bool mutate) { hasMutate_ = mutate; }

    static int mod(int a, int b) {
        assert(b > 0);

        while (a < 0) a += b;
        while (a >= b) a -= b;

        return a;
    }

// private:
    int id_;
    Threefry::Gen mut_prng_;
    int length_;

    double mutation_rate_;


    //--------------------------- Mutation counters
    int nb_swi_;
    int nb_mut_;

    int cpt_mut_;

    int min_genome_length_ = 10;
    int max_genome_length_ = 10000000;

    bool hasMutate_ = false;
};


#endif //RAEVOL_CUDA_DNAMUTATOR_H
