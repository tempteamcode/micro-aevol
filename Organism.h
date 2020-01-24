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
#include <atomic>

#include "RNA.h"
#include "Protein.h"
#include "Dna.h"

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
// * @param indiv_id : Unique Identification Number
 */
Organism(Threefry::Gen&& rng, int length)
: dna_(length, std::move(rng))
// , parent_length_(length)
//, indiv_id_(indiv_id)
{
  //ids_[0] = indiv_id;
  ids_[3] = length;
}

/**
 * Create an organism with a given genome
 *
 * @param length : Length of the generated random DNA
 * @param genome : Genome to assign to the organism
// * @param indiv_id : Unique Identification Number
 */
/*Organism(int length, char *genome)
: dna_(length, genome)
// , parent_length_(length)
//, indiv_id_(indiv_id)
{
  //ids_[0] = indiv_id;
  ids_[3] = length;
}*/


    Organism(const Organism& other) = delete;
    Organism(Organism&& other) = delete;

    inline Organism(const Organism& other, int)
    : dna_(other.dna_)
    , promoters_(other.promoters_)
    // , parent_length_(other.length())
    {
      ids_[3] = other.ids_[3];
    }

    inline Organism(Organism&& other, int)
    : dna_(std::move(other.dna_))
    , promoters_(std::move(other.promoters_))
    // , parent_length_(other.length())
    {
      ids_[3] = other.ids_[3];
    }


    Organism(gzFile backup_file);

    ~Organism() = default;

    void save(gzFile backup_file, int nb_indivs_) const;
    void load(gzFile backup_file);

    inline int length() const { return dna_.length(); };

    void apply_mutation(int pos);

    void start_stop_RNA();
    void compute_RNA();
    void opt_prom_compute_RNA();
    //void start_protein();
    //void compute_protein();
    void compute_proteins();
private:
    void compute_protein(RNA& rna, int protein_start);
public:
    void translate_protein(double w_max);

    // void compute_phenotype();
    // void compute_fitness(double selection_pressure, double* target);
    void compute_phenotype_fitness(double selection_pressure, double* target);

/**
 * Reset the stats variable (used when an organism is a perfect clone of its parent, it means no mutation)
 */
inline void reset_mutation_stats() { nb_swi_ = 0; nb_mut_ = 0; }
    void compute_protein_stats();

private:
    // Map of promoters (pos, err)
    std::map<int, int> promoters_;

    std::set<int> terminators;
    std::vector<RNA> rnas;
    std::vector<Protein> proteins;

public:
    double fitness;
    double metaerror;

    Dna dna_;

    int ids_[4];
    //int& indiv_id_ = ids_[0];
    //int& parent_id_ = ids_[1];
    //int& global_id_ = ids_[2];
    //int& parent_length_ = ids_[3];

    std::atomic<int> usage_count_{1};

    // Stats
    // int nb_genes_activ;
    // int nb_genes_inhib;
    int nb_func_genes;
    int nb_non_func_genes;
    int nb_coding_RNAs;
    int nb_non_coding_RNAs;
    int nb_swi_;
    int nb_mut_;

private:
    void remove_all_promoters();

    // void remove_promoters_around(int32_t pos);

    void remove_promoters_around(int32_t pos_1, int32_t pos_2);

    void remove_promoters_starting_between(int32_t pos_1, int32_t pos_2);

    void remove_promoters_starting_after(int32_t pos);

    void remove_promoters_starting_before(int32_t pos);

    // void locate_promoters();

    // void look_for_new_promoters_around(int32_t pos);

    void look_for_new_promoters_around(int32_t pos_1, int32_t pos_2);

    void look_for_new_promoters_starting_between(int32_t pos_1, int32_t pos_2);

    void look_for_new_promoters_starting_after(int32_t pos);

    void look_for_new_promoters_starting_before(int32_t pos);

    void add_new_promoter(int32_t position, int8_t error);

};


#endif //PDC_MINI_AEVOL_ORGANISM_H
