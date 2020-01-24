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



#include <cmath>
#include <cstring>
#include "Organism.h"
#include "ExpManager.h"

#include <iostream>

using namespace std;

/**
 * Constructor to create a clone of a given Organism
 *
 * @param clone : The organism to clone
 */
Organism::Organism(const Organism& other)
//: rna_count_(0)
: parent_length_(other.length())
, dna_(other.dna_)
, promoters_(other.promoters_)
//, terminators(other.terminators)
{
}

/**
 * Create an Organism from a backup/checkpointing file
 *
 * @param backup_file : gzFile to read from
 */
Organism::Organism(gzFile backup_file)
//: rna_count_(0)
{
    load(backup_file);
}

/**
 * Save the organism to backup/checkpointing file
 *
 * @param backup_file : where to the save the organism
 */
void Organism::save(gzFile backup_file) const {
    gzwrite(backup_file, &indiv_id_, sizeof(indiv_id_));
    gzwrite(backup_file, &parent_id_, sizeof(parent_id_));
    gzwrite(backup_file, &global_id, sizeof(global_id));

    gzwrite(backup_file, &parent_length_, sizeof(parent_length_));

    dna_.save(backup_file);
}

/**
 * Load the organism from backup/checkpointing file
 *
 * @param backup_file : from where restore the organism
 */
void Organism::load(gzFile backup_file) {
    gzread(backup_file, &indiv_id_, sizeof(indiv_id_));
    gzread(backup_file, &parent_id_, sizeof(parent_id_));
    gzread(backup_file, &global_id, sizeof(global_id));

    gzread(backup_file, &parent_length_, sizeof(parent_length_));

    dna_ = Dna();
    dna_.load(backup_file);
}

/**
 * Reset the stats variable (used when an organism is a perfect clone of its parent, it means no mutation)
 */
void Organism::reset_mutation_stats() {
    nb_swi_ = 0;
    nb_mut_ = 0;
}

void Organism::compute_protein_stats() {
    nb_genes_activ = 0;
    nb_genes_inhib = 0;
    nb_func_genes = 0;
    nb_non_func_genes = 0;
    nb_coding_RNAs = 0;
    nb_non_coding_RNAs = 0;

    for (const RNA& rna : rnas) {
    //for(const RNA* rna : starts_rnas)
        if (rna.is_coding_)
            nb_coding_RNAs++;
        else
            nb_non_coding_RNAs++;
    }

    for (const Protein& protein : proteins) {
        if (protein.is_functional) {
            nb_func_genes++;
        } else {
            nb_non_func_genes++;
        }
        if (protein.h > 0) {
            nb_genes_activ++;
        } else {
            nb_genes_inhib++;
        }
    }
}

/**
 * Switch the DNA base-pair at a given position
 *
 * @param pos : the position where to switch the base-pair
 * @return
 */

//void Organism::apply_mutation(int pos) {
void Organism::apply_mutation(vector<int> mutation_list) {
    this->mutation_list=mutation_list;

    //#pragma omp parallel for
    //for(int pos : mutation_list){
    for(int i=0;i<mutation_list.size();i++){
        int pos = mutation_list[i];
        dna_.do_switch(pos);

        // Remove promoters containing the switched base
        bool rem_prom = remove_promoters_around(pos, mod(pos + 1, length()));
        //bool rem_term = remove_terminators_around(pos, mod(pos + 1, length()));

        // Look for potential new promoters containing the switched base
        if (length() >= PROM_SIZE)
            look_for_new_promoters_around(pos, mod(pos + 1, length()));

        // Look for potential new terminators containing the switched base
        /*if(length() >= TERM_SIZE)
            look_for_new_terminators_around(pos, mod(pos + 1, length()));*/

        //#pragma omp atomic
        nb_swi_++;

        //#pragma omp atomic
        nb_mut_++;
    }

    //cout << "TEST: " << promoters_.size() << " : " << terminators.size() << endl;
}


/**
Optimize promoters search
 **/
bool Organism::remove_promoters_around(int32_t pos_1, int32_t pos_2) {
    if (mod(pos_1 - pos_2, dna_.length()) >= PROM_SIZE) {
        return remove_promoters_starting_between(mod(pos_1 - PROM_SIZE + 1,
                                              dna_.length()),
                                          pos_2);
    } else {
        remove_all_promoters();
        return true;
    }
}

void Organism::look_for_new_promoters_around(int32_t pos_1, int32_t pos_2) {
    if (dna_.length() >= PROM_SIZE) {
        look_for_new_promoters_starting_between(
                mod(pos_1 - PROM_SIZE + 1,
                    dna_.length()), pos_2);
    }
}

void Organism::remove_all_promoters() {
    promoters_.clear();
}

/** LEADING promoters **/
/** REMOVE **/
bool Organism::remove_promoters_starting_between(int32_t pos_1, int32_t pos_2) {
    if (pos_1 > pos_2) {
        return remove_promoters_starting_after(pos_1) || remove_promoters_starting_before(pos_2);
    } else {
        // suppression is in [pos1, pos_2[, pos_2 is excluded
        int init_size = promoters_.size();
        promoters_.erase(promoters_.lower_bound(pos_1), promoters_.upper_bound(pos_2-1));
        return init_size != promoters_.size();
    }
}

bool Organism::remove_promoters_starting_after(int32_t pos) {
    int init_size = promoters_.size();
    promoters_.erase(promoters_.lower_bound(pos), promoters_.end());
    return init_size != promoters_.size();
}

bool Organism::remove_promoters_starting_before(int32_t pos) {
    // suppression is in [0, pos[, pos is excluded
    int init_size = promoters_.size();
    promoters_.erase(promoters_.begin(), promoters_.upper_bound(pos-1));
    return init_size != promoters_.size();
}

void Organism::add_new_promoter(int32_t position, int8_t error) {
    // TODO: Insertion should not always occur, especially if promoter become better or worse ?
    // Promoters are deleted anyway if victim of mutation. the IF stays unnecessary
    if(promoters_.find(position) == promoters_.end())
        promoters_[position] = Promoter(position, error);
}

void Organism::look_for_new_promoters_starting_between(int32_t pos_1, int32_t pos_2) {
    // When pos_1 > pos_2, we will perform the search in 2 steps.
    // As positions  0 and dna_.length() are equivalent, it's preferable to
    // keep 0 for pos_1 and dna_.length() for pos_2.

    if (pos_1 >= pos_2) {
        look_for_new_promoters_starting_after(pos_1);
        look_for_new_promoters_starting_before(pos_2);
        return;
    }
    // Hamming distance of the sequence from the promoter consensus

    for (int32_t i = pos_1; i < pos_2; i++) {
        int8_t dist = dna_.promoter_at(i);

        if (dist <= 4) { // dist takes the hamming distance of the sequence from the consensus
            add_new_promoter(i, dist);
        }
    }
}

void Organism::look_for_new_promoters_starting_after(int32_t pos) {
    for (int32_t i = pos; i < dna_.length(); i++) {
        int dist = dna_.promoter_at(i);

        if (dist <= 4) { // dist takes the hamming distance of the sequence from the consensus
            add_new_promoter(i, dist);
        }
    }
}

void Organism::look_for_new_promoters_starting_before(int32_t pos) {
    // Hamming distance of the sequence from the promoter consensus

    for (int32_t i = 0; i < pos; i++) {
        int dist = dna_.promoter_at(i);
        if (dist <= 4) { // dist takes the hamming distance of the sequence from the consensus
            add_new_promoter(i, dist);
        }
    }
}


/**
 * Dealing with terminators
 */
void Organism::remove_all_terminators(){
    terminators.clear();
}

bool Organism::remove_terminators_around(int32_t pos_1, int32_t pos_2){
    if (mod(pos_1 - pos_2, dna_.length()) >= TERM_SIZE) {
        return remove_terminators_starting_between(mod(pos_1 - TERM_SIZE + 1,
                                              dna_.length()),
                                          pos_2);
    } else {
        remove_all_terminators();
        return true;
    }
}

bool Organism::remove_terminators_starting_between(int32_t pos_1, int32_t pos_2){
    if (pos_1 > pos_2) {
        return remove_terminators_starting_after(pos_1) || remove_terminators_starting_before(pos_2);
    } else {
        // suppression is in [pos1, pos_2[, pos_2 is excluded
        int init_size = terminators.size();
        terminators.erase(terminators.lower_bound(pos_1), terminators.upper_bound(pos_2-1));
        return init_size != terminators.size();
    }
}

bool Organism::remove_terminators_starting_after(int32_t pos){
    int init_size = terminators.size();
    terminators.erase(terminators.lower_bound(pos), terminators.end());
    return init_size != terminators.size();
}

bool Organism::remove_terminators_starting_before(int32_t pos){
    int init_size = terminators.size();
    terminators.erase(terminators.begin(), terminators.upper_bound(pos-1));
    return init_size != terminators.size();
}

void Organism::look_for_new_terminators_around(int32_t pos_1, int32_t pos_2){
    if (dna_.length() >= PROM_SIZE) {
        look_for_new_terminators_starting_between(
                mod(pos_1 - TERM_SIZE + 1,
                    dna_.length()), pos_2);
    }
}

void Organism::look_for_new_terminators_starting_between(int32_t pos_1, int32_t pos_2){
    if (pos_1 >= pos_2) {
        look_for_new_terminators_starting_after(pos_1);
        look_for_new_terminators_starting_before(pos_2);
        return;
    }
    // Hamming distance of the sequence from the promoter consensus

    for (int32_t i = pos_1; i < pos_2; i++) {
        int8_t dist = dna_.terminator_at(i);

        if (dist == 4) { // dist takes the hamming distance of the sequence from the consensus
            terminators.insert(i);
        }
    }
}

void Organism::look_for_new_terminators_starting_after(int32_t pos){
    for (int32_t i = pos; i < dna_.length(); i++) {
        int dist = dna_.terminator_at(i);

        if (dist == 4) { // dist takes the hamming distance of the sequence from the consensus
            terminators.insert(i);
        }
    }
}

void Organism::look_for_new_terminators_starting_before(int32_t pos){
    for (int32_t i = 0; i < pos; i++) {
        int dist = dna_.terminator_at(i);

        if (dist == 4) { // dist takes the hamming distance of the sequence from the consensus
            terminators.insert(i);
        }
    }
}


/**
 * Search for Promoters and Terminators (i.e. beginning and ending of a RNA) within the whole DNA of an Organism
 */
void Organism::start_stop_RNA() {
    if (dna_.length() < PROM_SIZE) return;

    //#pragma omp parallel for
    for (int dna_pos = 0; dna_pos < dna_.length(); dna_pos++) {
        int dist_lead = dna_.promoter_at(dna_pos);
        if (dist_lead <= 4) {
            //#pragma omp critical
            {
                add_new_promoter(dna_pos, dist_lead);
            }
        }

        // Computing if a terminator exists at that position
        int dist_term_lead = dna_.terminator_at(dna_pos);
        if (dist_term_lead == 4) {
            //#pragma omp critical
            {
                terminators.insert(dna_pos);
            }
	    }
    }
}

/**
 * Create the list of RNAs based on the found promoters and terminators on the DNA of an Organism
 */
void Organism::compute_RNA() {
    rnas.clear();
    rnas.reserve(promoters_.size());

    for (const auto &prom_pair: promoters_) {
        int prom_pos = prom_pair.first;

        int k = prom_pos + 22;
        k = k >= dna_.length()
            ? k - dna_.length()
            : k;

        auto it_rna_end = terminators.lower_bound(k);
        if (it_rna_end == terminators.end()) {
            it_rna_end = terminators.begin();
        }

        int rna_end = *it_rna_end + 10 >= dna_.length()
                      ? *it_rna_end + 10 - dna_.length()
                      : *it_rna_end + 10;

        int rna_length;
        if (prom_pos > rna_end)
            rna_length = dna_.length() - prom_pos + rna_end;
        else
            rna_length = rna_end - prom_pos;
        rna_length -= 21;

        if (rna_length >= 0) {
            rnas.emplace_back(
                prom_pos,
                rna_end,
                1.0 - std::fabs(((float) prom_pair.second.error)) / 5.0,
                rna_length);
            //rna_count_++;
        }
    }
}

/**
 * Optimize version that do not need to search the whole Dna for promoters
 */
void Organism::opt_prom_compute_RNA() {
    proteins.clear();
    rnas.clear();
    terminators.clear();

    rnas.reserve(promoters_.size());

    for (const auto &prom_pair: promoters_) {
        int prom_pos = prom_pair.first;

        /* Search for terminators */
        int cur_pos = prom_pos + 22;
        cur_pos = cur_pos >= dna_.length()
                  ? cur_pos - dna_.length()
                  : cur_pos;

        int start_pos = cur_pos;

        bool terminator_found = false;

        while (!terminator_found) {
            int term_dist_leading = dna_.terminator_at(cur_pos);

            if (term_dist_leading == 4)
                terminator_found = true;
            else {
                cur_pos = cur_pos + 1 >= dna_.length()
                          ? cur_pos + 1 - dna_.length()
                          : cur_pos + 1;

                if (cur_pos == start_pos) {
                    break;
                }
            }
        }

        if (terminator_found) {
            int32_t rna_end = cur_pos + 10 >= dna_.length()
                              ? cur_pos + 10 - dna_.length()
                              : cur_pos + 10;

            int32_t rna_length = 0;

            if (prom_pos > rna_end)
                rna_length = dna_.length() - prom_pos + rna_end;
            else
                rna_length = rna_end - prom_pos;

            rna_length -= 21;

            if (rna_length > 0) {
                rnas.emplace_back(
                    prom_pos,
                    rna_end,
                    1.0 - std::fabs(((float) prom_pair.second.error)) / 5.0,
                    rna_length);
                //rna_count_++;
            }
        }
    }
}

void Organism::opt_compute_RNA(){
    proteins.clear();
    rnas.clear();

    rnas.reserve(promoters_.size());

    set<int>::iterator term = terminators.begin();

    for(const auto &prom_pair: promoters_) {
        int prom_pos = prom_pair.first;

        int cur_pos = prom_pos + 22;
        cur_pos = cur_pos >= dna_.length()
                  ? cur_pos - dna_.length()
                  : cur_pos;

        while(cur_pos > *term && term != terminators.end())
            ++term;

        if(term == terminators.end())
            term = terminators.begin();

        int start_pos = *term;

        int32_t rna_end = start_pos + 10 >= dna_.length()
                          ? start_pos + 10 - dna_.length()
                          : start_pos + 10;

        int32_t rna_length = 0;

        if (prom_pos > rna_end)
            rna_length = dna_.length() - prom_pos + rna_end;
        else
            rna_length = rna_end - prom_pos;

        rna_length -= 21;

        if (rna_length > 0) {
            rnas.emplace_back(
                    prom_pos,
                    rna_end,
                    1.0 - std::fabs(((float) prom_pair.second.error)) / 5.0,
                    rna_length);
        }
    }
}

/**
 * Search for Shine Dal sequence and Start sequence deliminating the start of genes within one of the RNA of an Organism
 */
void Organism::start_protein() {
    for (RNA& rna : rnas) {
        if (rna.length < 22) continue;

        int c_pos = rna.begin + 22;
        c_pos = c_pos >= dna_.length()
                ? c_pos - dna_.length()
                : c_pos;

        while (c_pos != rna.end) {
            if (dna_.shine_dal_start(c_pos)) {
                rna.start_prot.push_back(c_pos);
            }

            c_pos++;
            c_pos = c_pos >= dna_.length()
                    ? c_pos - dna_.length()
                    : c_pos;
        }
    }
}

/**
 * Compute the list of genes/proteins of an Organism
 */
void Organism::compute_protein() {
    int resize_to = 0;
    for (const RNA& rna : rnas) {
        resize_to += rna.start_prot.size();
    }

    proteins.clear();
    proteins.reserve(resize_to);

    for (RNA& rna : rnas) {
        for (int protein_idx = 0; protein_idx < rna.start_prot.size(); protein_idx++) {
            int protein_start = rna.start_prot[protein_idx];
            int current_position = protein_start + 13;

            current_position = current_position >= dna_.length()
                               ? current_position - dna_.length()
                               : current_position;

            int transcribed_start = rna.begin + 22;
            transcribed_start = transcribed_start >= dna_.length()
                                ? transcribed_start - dna_.length()
                                : transcribed_start;

            int transcription_length;
            if (transcribed_start <= protein_start) {
                transcription_length = protein_start - transcribed_start;
            } else {
                transcription_length = dna_.length() - transcribed_start + protein_start;
            }
            transcription_length += 13;

            while (rna.length - transcription_length >= 3) {
                if (dna_.protein_stop(current_position)) {
                    int prot_length;

                    int protein_end = current_position + 2 >= dna_.length() ?
                                      current_position - dna_.length() + 2 :
                                      current_position + 2;

                    if (protein_start + 13 < protein_end) {
                        prot_length = protein_end - (protein_start + 13);
                    } else {
                        prot_length = dna_.length() - (protein_start + 13) + protein_end;
                    }

                    if (prot_length >= 3) {
                        proteins.emplace_back
                                           (protein_start,
                                            protein_end,
                                            prot_length,
                                            rna.e);
                        //protein_count_++;

                        rna.is_coding_ = true;
                    }
                    break;
                }

                current_position += 3;
                current_position = current_position >= dna_.length()
                                   ? current_position - dna_.length()
                                   : current_position;
                transcription_length += 3;
            }
        }
    }
}

/**
 * Compute the pseudo-chimical model (i.e. the width, height and location in the phenotypic space) of a genes/protein
 *
 * @param w_max : Maximum width of the triangle generated by a Protein
 */
void Organism::translate_protein(double w_max) {
    for (Protein& protein : proteins) {
        if (!protein.is_init_) continue;

        int c_pos = protein.protein_start;
        c_pos += 13;
        c_pos = c_pos >= dna_.length()
                ? c_pos - dna_.length()
                : c_pos;

        int codon_list[64] = {};
        int codon_idx = 0;
        int count_loop = 0;

        //printf("Codon list : ");
        while (count_loop < protein.protein_length / 3 &&
               codon_idx < 64) {
            codon_list[codon_idx] = dna_.codon_at(c_pos);
            //printf("%d ",codon_list[codon_idx]);
            codon_idx++;

            count_loop++;
            c_pos += 3;
            c_pos = c_pos >= dna_.length()
                    ? c_pos - dna_.length()
                    : c_pos;
        }
        //printf("\n");

        double M = 0.0;
        double W = 0.0;
        double H = 0.0;

        int nb_m = 0;
        int nb_w = 0;
        int nb_h = 0;

        bool bin_m = false; // Initializing to false will yield a conservation of the high weight bit
        bool bin_w = false; // when applying the XOR operator for the Gray to standard conversion
        bool bin_h = false;


        for (int i = 0; i < codon_idx; i++) {
            switch (codon_list[i]) {
                case CODON_M0 : {
                    // M codon found
                    nb_m++;

                    // Convert Gray code to "standard" binary code
                    bin_m ^= false; // as bin_m was initialized to false, the XOR will have no effect on the high weight bit

                    // A lower-than-the-previous-lowest weight bit was found, make a left bitwise shift
                    //~ M <<= 1;
                    M *= 2;

                    // Add this nucleotide's contribution to M
                    if (bin_m) M += 1;

                    break;
                }
                case CODON_M1 : {
                    // M codon found
                    nb_m++;

                    // Convert Gray code to "standard" binary code
                    bin_m ^= true; // as bin_m was initialized to false, the XOR will have no effect on the high weight bit

                    // A lower-than-the-previous-lowest bit was found, make a left bitwise shift
                    //~ M <<= 1;
                    M *= 2;

                    // Add this nucleotide's contribution to M
                    if (bin_m) M += 1;

                    break;
                }
                case CODON_W0 : {
                    // W codon found
                    nb_w++;

                    // Convert Gray code to "standard" binary code
                    bin_w ^= false; // as bin_m was initialized to false, the XOR will have no effect on the high weight bit

                    // A lower-than-the-previous-lowest weight bit was found, make a left bitwise shift
                    //~ W <<= 1;
                    W *= 2;

                    // Add this nucleotide's contribution to W
                    if (bin_w) W += 1;

                    break;
                }
                case CODON_W1 : {
                    // W codon found
                    nb_w++;

                    // Convert Gray code to "standard" binary code
                    bin_w ^= true; // as bin_m was initialized to false, the XOR will have no effect on the high weight bit

                    // A lower-than-the-previous-lowest weight bit was found, make a left bitwise shift
                    //~ W <<= 1;
                    W *= 2;

                    // Add this nucleotide's contribution to W
                    if (bin_w) W += 1;

                    break;
                }
                case CODON_H0 :
                case CODON_START : // Start codon codes for the same amino-acid as H0 codon
                {
                    // H codon found
                    nb_h++;

                    // Convert Gray code to "standard" binary code
                    bin_h ^= false; // as bin_m was initialized to false, the XOR will have no effect on the high weight bit

                    // A lower-than-the-previous-lowest weight bit was found, make a left bitwise shift
                    //~ H <<= 1;
                    H *= 2;

                    // Add this nucleotide's contribution to H
                    if (bin_h) H += 1;

                    break;
                }
                case CODON_H1 : {
                    // H codon found
                    nb_h++;

                    // Convert Gray code to "standard" binary code
                    bin_h ^= true; // as bin_m was initialized to false, the XOR will have no effect on the high weight bit

                    // A lower-than-the-previous-lowest weight bit was found, make a left bitwise shift
                    //~ H <<= 1;
                    H *= 2;

                    // Add this nucleotide's contribution to H
                    if (bin_h) H += 1;

                    break;
                }
            }
        }

        protein.protein_length = codon_idx;


        //  ----------------------------------------------------------------------------------
        //  2) Normalize M, W and H values in [0;1] according to number of codons of each kind
        //  ----------------------------------------------------------------------------------
        protein.m =
                nb_m != 0 ? M / (pow(2, nb_m) - 1) : 0.5;
        protein.w =
                nb_w != 0 ? W / (pow(2, nb_w) - 1) : 0.0;
        protein.h =
                nb_h != 0 ? H / (pow(2, nb_h) - 1) : 0.5;

        //  ------------------------------------------------------------------------------------
        //  3) Normalize M, W and H values according to the allowed ranges (defined in macros.h)
        //  ------------------------------------------------------------------------------------
        // x_min <= M <= x_max
        // w_min <= W <= w_max
        // h_min <= H <= h_max
        protein.m =
                (X_MAX - X_MIN) *
                protein.m +
                X_MIN;
        protein.w =
                (w_max - W_MIN) *
                protein.w +
                W_MIN;
        protein.h =
                (H_MAX - H_MIN) *
                protein.h +
                H_MIN;

        protein.is_functional = !(nb_m == 0 || nb_w == 0 || nb_h == 0 || protein.w == 0.0 || protein.h == 0.0);
    }


    std::map<int, Protein *> lookup;

    for (Protein& protein : proteins) {
        if (protein.is_init_) {
            if (lookup.find(protein.protein_start) == lookup.end()) {
                lookup[protein.protein_start] = &protein;
            } else {
                lookup[protein.protein_start]->e += protein.e;
                protein.is_init_ = false;
            }
        }
    }
}

/**
 * From the list of proteins, build the phenotype of an organism
 */
/*
void Organism::compute_phenotype() {
    double activ_phenotype[300]{};
    double inhib_phenotype[300]{};

    for (const Protein& protein : proteins) {
        if (protein.is_init_ &&
            fabs(protein.w) >= 1e-15 &&
            fabs(protein.h) >= 1e-15 &&
            protein.is_functional) {
            // Compute triangle points' coordinates
            double x0 =
                    protein.m -
                    protein.w;

            double x1 = protein.m;
            double x2 =
                    protein.m +
                    protein.w;

            int ix0 = (int) (x0 * 300);
            int ix1 = (int) (x1 * 300);
            int ix2 = (int) (x2 * 300);

            if (ix0 < 0) ix0 = 0; else if (ix0 > (299)) ix0 = 299;
            if (ix1 < 0) ix1 = 0; else if (ix1 > (299)) ix1 = 299;
            if (ix2 < 0) ix2 = 0; else if (ix2 > (299)) ix2 = 299;

            // Compute the first equation of the triangle
            double incY =
                    (protein.h *
                     protein.e) /
                    (ix1 - ix0);
            int count = 1;

            // Updating value between x0 and x1
            for (int i = ix0 + 1; i < ix1; i++) {
                if (protein.h > 0)
                    activ_phenotype[i] += (incY * (count++));
                else
                    inhib_phenotype[i] += (incY * (count++));
            }


            if (protein.h > 0)
                activ_phenotype[ix1] += (protein.h *
                                         protein.e);
            else
                inhib_phenotype[ix1] += (protein.h *
                                         protein.e);


            // Compute the second equation of the triangle
            incY = (protein.h *
                    protein.e) /
                   (ix2 - ix1);
            count = 1;

            // Updating value between x1 and x2
            for (int i = ix1 + 1; i < ix2; i++) {
                if (protein.h > 0)
                    activ_phenotype[i] += ((protein.h *
                                            protein.e) -
                                           (incY * (count++)));
                else
                    inhib_phenotype[i] += ((protein.h *
                                            protein.e) -
                                           (incY * (count++)));
            }
        }
    }

    for (int fuzzy_idx = 0; fuzzy_idx < 300; fuzzy_idx++) {
        if (activ_phenotype[fuzzy_idx] > 1)
            activ_phenotype[fuzzy_idx] = 1;
        if (inhib_phenotype[fuzzy_idx] < -1)
            inhib_phenotype[fuzzy_idx] = -1;
    }

    for (int fuzzy_idx = 0; fuzzy_idx < 300; fuzzy_idx++) {
        phenotype[fuzzy_idx] = activ_phenotype[fuzzy_idx] + inhib_phenotype[fuzzy_idx];
        if (phenotype[fuzzy_idx] < 0)
            phenotype[fuzzy_idx] = 0;
    }
}
*/

/**
 * From the phenotype of an organism, compute its metabolic error and fitness
 *
 * @param selection_pressure : Selection pressure used during the selection process
 */
/*
void Organism::compute_fitness(double selection_pressure, double* target) {
    for (int fuzzy_idx = 0; fuzzy_idx < 300; fuzzy_idx++) {
        if (phenotype[fuzzy_idx] > 1)
            phenotype[fuzzy_idx] = 1;
        if (phenotype[fuzzy_idx] < 0)
            phenotype[fuzzy_idx] = 0;

        delta[fuzzy_idx] = phenotype[fuzzy_idx] - target[fuzzy_idx];
    }

    metaerror = 0;
    for (int fuzzy_idx = 0; fuzzy_idx < 299; fuzzy_idx++) {
        metaerror +=
                ((std::fabs(delta[fuzzy_idx]) +
                  std::fabs(delta[fuzzy_idx + 1])) /
                 (600.0));
    }

    fitness = exp(-selection_pressure * ((double) metaerror));
}
*/


/**
 * From the phenotype of an organism, compute its metabolic error and fitness
 *
 * @param selection_pressure : Selection pressure used during the selection process
 */
void Organism::compute_fitness(double selection_pressure, double* target) {
    double activ_phenotype[300]{};
    double inhib_phenotype[300]{};

    for (const Protein& protein : proteins) {
        if (protein.is_init_ &&
            fabs(protein.w) >= 1e-15 &&
            fabs(protein.h) >= 1e-15 &&
            protein.is_functional) {
            // Compute triangle points' coordinates
            double x0 = protein.m - protein.w;
            double x1 = protein.m;
            double x2 = protein.m + protein.w;

            double h_times_e = protein.h * protein.e;

            int ix0 = (int) (x0 * 300);
            int ix1 = (int) (x1 * 300);
            int ix2 = (int) (x2 * 300);

            if (ix0 < 0) ix0 = 0; else if (ix0 > (299)) ix0 = 299;
            if (ix1 < 0) ix1 = 0; else if (ix1 > (299)) ix1 = 299;
            if (ix2 < 0) ix2 = 0; else if (ix2 > (299)) ix2 = 299;

            // Compute the first equation of the triangle
            double incY = h_times_e / (ix1 - ix0);
            int count = 1;

            // Updating value between x0 and x1
            for (int i = ix0 + 1; i < ix1; i++) {
                if (protein.h > 0)
                    activ_phenotype[i] += (incY * (count++));
                else
                    inhib_phenotype[i] += (incY * (count++));
            }

            if (protein.h > 0)
                activ_phenotype[ix1] += h_times_e;
            else
                inhib_phenotype[ix1] += h_times_e;

            // Compute the second equation of the triangle
            incY = h_times_e / (ix2 - ix1);
            count = 1;

            // Updating value between x1 and x2
            for (int i = ix1 + 1; i < ix2; i++) {
                if (protein.h > 0)
                    activ_phenotype[i] += (h_times_e - (incY * (count++)));
                else
                    inhib_phenotype[i] += (h_times_e - (incY * (count++)));
            }
        }
    }

    for (int fuzzy_idx = 0; fuzzy_idx < 300; fuzzy_idx++) {
        if (activ_phenotype[fuzzy_idx] > 1)
            activ_phenotype[fuzzy_idx] = 1;
        if (inhib_phenotype[fuzzy_idx] < -1)
            inhib_phenotype[fuzzy_idx] = -1;
    }

    metaerror = 0;

    double abs_delta_prev;
    double abs_delta_curr;
    for (int fuzzy_idx = 0; fuzzy_idx < 300; fuzzy_idx++) {
        abs_delta_prev = abs_delta_curr;
        
        abs_delta_curr = activ_phenotype[fuzzy_idx] + inhib_phenotype[fuzzy_idx];
        if (abs_delta_curr < 0) abs_delta_curr = 0;
        if (abs_delta_curr > 1) abs_delta_curr = 1;

        abs_delta_curr -= target[fuzzy_idx];
        abs_delta_curr = std::fabs(abs_delta_curr);

        if (fuzzy_idx != 0) {
            metaerror += (abs_delta_prev + abs_delta_curr) / (600.0);
        }
    }

    fitness = exp(-selection_pressure * ((double) metaerror));
}


