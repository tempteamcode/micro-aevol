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


#ifndef PDC_MINI_AEVOL_ORGANISM_H
#define PDC_MINI_AEVOL_ORGANISM_H

#include <map>
#include <memory>
#include <set>
#include <zlib.h>
#include <list>
#include <unordered_map>

#include "Promoter.h"
#include "RNA.h"
#include "Protein.h"
#include "Dna.h"

//class ExpManager;

/**
 * Class that implements an organism and its related DNAs, RNAs, Protein and Phenotype
 */
class Organism {

public:

/**
 * Constructor to generate a random organism (i.e. an organism with a random DNA)
 *
 * @param rng : random generator
 * @param length : Length of the generated random DNA
 * @param indiv_id : Unique Identification Number
 */
Organism(Threefry::Gen&& rng, int length, int indiv_id)
//: rna_count_(0)
: dna_(length, std::move(rng))
, parent_length_(length)
, indiv_id_(indiv_id)
{
}

/**
 * Create an organism with a given genome
 *
 * @param length : Length of the generated random DNA
 * @param genome : Genome to assign to the organism
 * @param indiv_id : Unique Identification Number
 */
/*Organism(int length, char *genome, int indiv_id)
//: rna_count_(0)
: parent_length_(length)
, dna_(length, genome)
, indiv_id_(indiv_id)
{
}*/

    Organism(const Organism& other);

    Organism(gzFile backup_file);

    ~Organism() = default;

    void save(gzFile backup_file) const;
    void load(gzFile backup_file);

    inline int length() const { return dna_.length(); };

    //void apply_mutation(int pos);
    void apply_mutation(std::vector<int> mutation_list);

    void start_stop_RNA();
    void compute_RNA();
    void opt_prom_compute_RNA();
    void start_protein();
    void compute_protein();
    void translate_protein(double w_max);

    inline void compute_phenotype() { } // void compute_phenotype();
    void compute_fitness(double selection_pressure, double* target);

    void reset_mutation_stats();
    void compute_protein_stats();

    void opt_compute_RNA();

//private:
    // Map position (int) to Promoter
    std::map<int, Promoter> promoters_;

    std::set<int> terminators;
    std::vector<RNA> rnas;

    std::vector<Protein> proteins;

    std::vector<int> mutation_list;

public:
    double fitness;
    double metaerror;

    Dna dna_;
    int parent_length_;

    int indiv_id_;
    int parent_id_;

    int global_id = -1;

    int usage_count_ = 1;

    // Stats
    int nb_genes_activ = 0;
    int nb_genes_inhib = 0;
    int nb_func_genes = 0;
    int nb_non_func_genes = 0;
    int nb_coding_RNAs = 0;
    int nb_non_coding_RNAs = 0;

    int nb_swi_ = 0;
    int nb_mut_ = 0;

private:
    void remove_all_promoters();

    bool remove_promoters_around(int32_t pos_1, int32_t pos_2);

    bool remove_promoters_starting_between(int32_t pos_1, int32_t pos_2);

    bool remove_promoters_starting_after(int32_t pos);

    bool remove_promoters_starting_before(int32_t pos);

    void look_for_new_promoters_around(int32_t pos_1, int32_t pos_2);

    void look_for_new_promoters_starting_between(int32_t pos_1, int32_t pos_2);

    void look_for_new_promoters_starting_after(int32_t pos);

    void look_for_new_promoters_starting_before(int32_t pos);

    void add_new_promoter(int32_t position, int8_t error);

    void remove_all_terminators();

    bool remove_terminators_around(int32_t pos_1, int32_t pos_2);

    bool remove_terminators_starting_between(int32_t pos_1, int32_t pos_2);

    bool remove_terminators_starting_after(int32_t pos);

    bool remove_terminators_starting_before(int32_t pos);

    void look_for_new_terminators_around(int32_t pos_1, int32_t pos_2);

    void look_for_new_terminators_starting_between(int32_t pos_1, int32_t pos_2);

    void look_for_new_terminators_starting_after(int32_t pos);

    void look_for_new_terminators_starting_before(int32_t pos);

    inline int32_t mod(int32_t a, int32_t b) {
        assert(b > 0);

        while (a < 0) a += b;
        while (a >= b) a -= b;

        return a;
        //return m >= 0 ? m % n : ( n - abs ( m%n ) ) % n;
    }

    inline int64_t mod(int64_t a, int64_t b) {
        assert(b > 0);

        while (a < 0) a += b;
        while (a >= b) a -= b;

        return a;
        //return m >= 0 ? m % n : ( n - abs ( m%n ) ) % n;
    }
};


#endif //PDC_MINI_AEVOL_ORGANISM_H
