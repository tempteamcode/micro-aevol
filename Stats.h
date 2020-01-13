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


#ifndef PDC_MINI_AEVOL_STATS_H
#define PDC_MINI_AEVOL_STATS_H

#include <cstdint>
#include <fstream>
#include <limits>

class ExpManager;
class Organism;

/**
 * Class to manage and generate the Stats (and the related file) of a simulation
 */
class Stats {
public:
    inline Stats() : generation_(0) { }

    Stats(int generation, bool best_or_not);

    ~Stats() {
        if (is_indiv_) {
            statfile_best_.flush();
            statfile_best_.close();
        } else {
            statfile_mean_.flush();
            statfile_mean_.close();
        }
    }

    void compute_best(const Organism& best_indiv);
    void compute_average(ExpManager& exp_m);

    void write_best();
    void write_average();

    void reinit(int generation, bool is_indiv);
    void prepare();

    bool is_indiv() const { return is_indiv_; }

private:
    int generation_;

    bool is_indiv_;

    int pop_size_;

    double fitness_;
    double mean_fitness_;
    double metabolic_error_;
    double mean_metabolic_error_;

    int amount_of_dna_;
    float mean_amount_of_dna_;
    int nb_coding_rnas_;
    float mean_nb_coding_rnas_;
    int nb_non_coding_rnas_;
    float mean_nb_non_coding_rnas_;

    int nb_functional_genes_;
    float mean_nb_functional_genes_;
    int nb_non_functional_genes_;
    float mean_nb_non_functional_genes_;

    int nb_mut_;
    float mean_nb_mut_;
    int nb_switch_;
    float mean_nb_switch_;

    bool is_computed_ = false;

    // Stats

    std::ofstream statfile_best_;
    std::ofstream statfile_mean_;
};


#endif //PDC_MINI_AEVOL_STATS_H
