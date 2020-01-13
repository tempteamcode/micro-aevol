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

#include "DnaMutator.h"
#include "Organism.h"

/**
 * Generate both type of the mutations (see below)
 */
void DnaMutator::generate_mutations() {
    mutation_list_.clear();
    nb_swi_ = mut_prng_.binomial_random(length_, mutation_rate_);
    nb_mut_ = nb_swi_;

    for (int cpt_mut_ = nb_mut_; cpt_mut_ > 0; cpt_mut_--) {
        generate_next_mutation(cpt_mut_);
    }
}

/**
 * Generate the next mutation event for an organism.
 */
void DnaMutator::generate_next_mutation(int cpt_mut_) {
    int random_value = mut_prng_.random(cpt_mut_);

    if (random_value < nb_swi_) {
        nb_swi_--;

        int pos = mut_prng_.random(length_);
        mutation_list_.push_back(pos);
    }
}

void DnaMutator::apply_mutations(Organism& organism) {
    for (int mutation_pos: mutation_list_) {
        organism.apply_mutation(mutation_pos);
    }
}

