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
#include <cstdint>
#include <cstring>
#include <zlib.h>
#include "Organism.h"
#include "ExpManager.h"

inline int32_t mod_once(int32_t a, int32_t b) {
    return a < 0 ? a + b : (a >= b ? a - b : a);
}
inline int64_t mod_once(int64_t a, int64_t b) {
    return a < 0 ? a + b : (a >= b ? a - b : a);
}

void Organism::compute_protein_stats() {
    // nb_genes_activ = 0;
    // nb_genes_inhib = 0;
    nb_func_genes = 0;
    nb_non_func_genes = 0;
    nb_coding_RNAs = 0;
    nb_non_coding_RNAs = 0;

    for (const RNA& rna : rnas) {
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
    }
}

/**
 * Switch the DNA base-pair at a given position
 *
 * @param pos : the position where to switch the base-pair
 * @return
 */

void Organism::apply_mutation(std::vector<int> mut_list) {
    const int dna_length = dna_.seq_.size();

    for(int i=0;i<mut_list.size();i++){
        int pos = mut_list[i];
        dna_.do_switch(pos);

        // Remove promoters containing the switched base
        remove_promoters_around(pos, mod_once(pos + 1, dna_length));

        // Look for potential new promoters containing the switched base
        if (dna_length >= PROM_SIZE)
            look_for_new_promoters_around(pos, mod_once(pos + 1, dna_length));

        nb_swi_++;
        nb_mut_++;
    }
}


/**
Optimize promoters search
 **/

void Organism::remove_promoters_around(int32_t pos_1, int32_t pos_2) {
    const int dna_length = dna_.seq_.size();

    if (mod_once(pos_1 - pos_2, dna_length) >= PROM_SIZE) {
        remove_promoters_starting_between(mod_once(pos_1 - PROM_SIZE + 1, dna_length), pos_2);
    } else {
        remove_all_promoters();
    }
}

void Organism::look_for_new_promoters_around(int32_t pos_1, int32_t pos_2) {
    const int dna_length = dna_.seq_.size();

    if (dna_length >= PROM_SIZE) {
        look_for_new_promoters_starting_between(mod_once(pos_1 - PROM_SIZE + 1, dna_length), pos_2);
    }
}

void Organism::remove_all_promoters() {
    promoters_.clear();
}

/** LEADING promoters **/
/** REMOVE **/
void Organism::remove_promoters_starting_between(int32_t pos_1, int32_t pos_2) {
    if (pos_1 > pos_2) {
        remove_promoters_starting_after(pos_1);
        remove_promoters_starting_before(pos_2);
    } else {
        // suppression is in [pos1, pos_2[, pos_2 is excluded
        promoters_.erase(promoters_.lower_bound(pos_1), promoters_.upper_bound(pos_2-1));
    }
}

void Organism::remove_promoters_starting_after(int32_t pos) {
    promoters_.erase(promoters_.lower_bound(pos), promoters_.end());
}

void Organism::remove_promoters_starting_before(int32_t pos) {
    // suppression is in [0, pos[, pos is excluded
    promoters_.erase(promoters_.begin(), promoters_.upper_bound(pos-1));
}

void Organism::look_for_new_promoters_starting_between(int32_t pos_1, int32_t pos_2) {
    // When pos_1 > pos_2, we will perform the search in 2 steps.
    // As positions  0 and dna_length are equivalent, it's preferable to
    // keep 0 for pos_1 and dna_length for pos_2.

    if (pos_1 >= pos_2) {
        look_for_new_promoters_starting_after(pos_1);
        look_for_new_promoters_starting_before(pos_2);
        return;
    }
    // Hamming distance of the sequence from the promoter consensus

    dna_.seq_.forSequences(pos_1, pos_2-pos_1, PROM_SEQ_LEN, [&] (size_t pos_plus_len, int_t sequence) {
       int dist = hamming_distance(sequence, PROM_SEQ);
       if (dist <= 4) { // dist takes the hamming distance of the sequence from the consensus
         promoters_[pos_plus_len < PROM_SEQ_LEN ? dna_.seq_.size() + pos_plus_len - PROM_SEQ_LEN : pos_plus_len - PROM_SEQ_LEN] = dist;
       }
    });
}

void Organism::look_for_new_promoters_starting_after(int32_t pos) {
    dna_.seq_.forSequences(pos, dna_.seq_.size() - pos, PROM_SEQ_LEN, [&] (size_t pos_plus_len, int_t sequence) {
       int dist = hamming_distance(sequence, PROM_SEQ);
       if (dist <= 4) { // dist takes the hamming distance of the sequence from the consensus
         promoters_[pos_plus_len < PROM_SEQ_LEN ? dna_.seq_.size() + pos_plus_len - PROM_SEQ_LEN : pos_plus_len - PROM_SEQ_LEN] = dist;
       }
    });
}

void Organism::look_for_new_promoters_starting_before(int32_t pos) {
    // Hamming distance of the sequence from the promoter consensus

    dna_.seq_.forSequences(0, pos, PROM_SEQ_LEN, [&] (size_t pos_plus_len, int_t sequence) {
       int dist = hamming_distance(sequence, PROM_SEQ);
       if (dist <= 4) { // dist takes the hamming distance of the sequence from the consensus
         promoters_[pos_plus_len < PROM_SEQ_LEN ? dna_.seq_.size() + pos_plus_len - PROM_SEQ_LEN : pos_plus_len - PROM_SEQ_LEN] = dist;
       }
    });
}


/**
 * Search for Promoters and Terminators (i.e. beginning and ending of a RNA) within the whole DNA of an Organism
 */
void Organism::start_stop_RNA() {
    if (dna_.seq_.size() < PROM_SIZE) return;

    dna_.seq_.forSequences(0, dna_.seq_.size(), PROM_SEQ_LEN, [&] (size_t pos_plus_len, int_t sequence) {
        int dist = hamming_distance(sequence, PROM_SEQ);
        if (dist <= 4) { // dist takes the hamming distance of the sequence from the consensus
            int dna_pos = pos_plus_len < PROM_SEQ_LEN ? dna_.seq_.size() + pos_plus_len - PROM_SEQ_LEN : pos_plus_len - PROM_SEQ_LEN;
            promoters_[dna_pos] = dist;
        }

        // Computing if a terminator exists at that position
        int_t subseq_rev = lookuptable[(sequence >> 7) & 0xf0];
        int_t getseq_rev = (sequence >> 22-4);
        int dist_term_lead = hamming_distance(getseq_rev, subseq_rev);
        if (dist_term_lead == 4) {
            int dna_pos = pos_plus_len < PROM_SEQ_LEN ? dna_.seq_.size() + pos_plus_len - PROM_SEQ_LEN : pos_plus_len - PROM_SEQ_LEN;
            terminators.insert(dna_pos);
        }
    });
}


/**
 * Create the list of RNAs based on the found promoters and terminators on the DNA of an Organism
 */
void Organism::compute_RNA() {
    const int dna_length = dna_.seq_.size();

    rnas.clear();
    rnas.reserve(promoters_.size());

    for (const auto &prom_pair: promoters_) {
        int prom_pos = prom_pair.first;

        int k = prom_pos + 22;
        k = k >= dna_length
            ? k - dna_length
            : k;

        auto it_rna_end = terminators.lower_bound(k);
        if (it_rna_end == terminators.end()) {
            it_rna_end = terminators.begin();
        }

        int rna_end = *it_rna_end + 10 >= dna_length
                      ? *it_rna_end + 10 - dna_length
                      : *it_rna_end + 10;

        int rna_length;
        if (prom_pos > rna_end)
            rna_length = dna_length - prom_pos + rna_end;
        else
            rna_length = rna_end - prom_pos;
        rna_length -= 21;

        if (rna_length >= 0) {
            rnas.emplace_back(
                prom_pos,
                rna_end,
                1.0 - std::fabs(((float) prom_pair.second)) / 5.0,
                rna_length);
            //rna_count_++;
        }
    }
}

/**
 * Optimize version that do not need to search the whole Dna for promoters
 */
void Organism::opt_prom_compute_RNA() {
    const int dna_length = dna_.seq_.size();

    proteins.clear();
    rnas.clear();
    terminators.clear();

    rnas.reserve(promoters_.size());

    for (const auto &prom_pair: promoters_) {
        int prom_pos = prom_pair.first;

        /* Search for terminators */
        int cur_pos = prom_pos + 22;
        cur_pos = cur_pos >= dna_length
                  ? cur_pos - dna_length
                  : cur_pos;

        int start_pos = cur_pos;

        bool terminator_found = false;

        while (!terminator_found) {
            int term_dist_leading = dna_.terminator_at(cur_pos);

            if (term_dist_leading == 4)
                terminator_found = true;
            else {
                cur_pos = cur_pos + 1 >= dna_length
                          ? cur_pos + 1 - dna_length
                          : cur_pos + 1;

                if (cur_pos == start_pos) {
                    break;
                }
            }
        }

        if (terminator_found) {
            int32_t rna_end = cur_pos + 10 >= dna_length
                              ? cur_pos + 10 - dna_length
                              : cur_pos + 10;

            int32_t rna_length = 0;

            if (prom_pos > rna_end)
                rna_length = dna_length - prom_pos + rna_end;
            else
                rna_length = rna_end - prom_pos;

            rna_length -= 21;

            if (rna_length > 0) {
                rnas.emplace_back(
                    prom_pos,
                    rna_end,
                    1.0 - std::fabs(((float) prom_pair.second)) / 5.0,
                    rna_length);
                //rna_count_++;
            }
        }
    }
}


/**
 * Search for Shine Dal sequence and Start sequence deliminating the start of genes within one of the RNA of an Organism
 * Compute the list of genes/proteins of an Organism
 */
void Organism::compute_proteins() {
    const int dna_length = dna_.seq_.size();

    proteins.clear();

    for (RNA& rna : rnas) {
        int c_pos = rna.begin;
        int c_len = rna.length;
        if (c_len > PROM_SEQ_LEN) {
            c_pos += PROM_SEQ_LEN;
            if (c_pos >= dna_length) c_pos -= dna_length;
            
            dna_.seq_.forSequences(c_pos, c_len - PROM_SEQ_LEN, SHINE_DAL_SEQ_LEN, [&] (size_t pos_plus_len, int_t sequence) {
                if ((sequence & SHINE_DAL_SEQ_MASK) == SHINE_DAL_SEQ_BITS) {
                    int dna_pos = pos_plus_len < SHINE_DAL_SEQ_LEN ? dna_length + pos_plus_len - SHINE_DAL_SEQ_LEN : pos_plus_len - SHINE_DAL_SEQ_LEN;
                    compute_protein(rna, dna_pos);
                }
            });
        }
    }
}

void Organism::compute_protein(RNA& rna, int protein_start) {
    const int dna_length = dna_.seq_.size();

    int current_position = protein_start + 13;

    current_position = current_position >= dna_length
                       ? current_position - dna_length
                       : current_position;

    int transcribed_start = rna.begin + 22;
    transcribed_start = transcribed_start >= dna_length
                        ? transcribed_start - dna_length
                        : transcribed_start;

    int transcription_length;
    if (transcribed_start <= protein_start) {
        transcription_length = protein_start - transcribed_start;
    } else {
        transcription_length = dna_length - transcribed_start + protein_start;
    }
    transcription_length += 13;

    while (rna.length - transcription_length >= 3) {
        if (dna_.protein_stop(current_position)) {
            int prot_length;

            int protein_end = current_position + 2 >= dna_length ?
                              current_position - dna_length + 2 :
                              current_position + 2;

            if (protein_start + 13 < protein_end) {
                prot_length = protein_end - (protein_start + 13);
            } else {
                prot_length = dna_length - (protein_start + 13) + protein_end;
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
        current_position = current_position >= dna_length
                           ? current_position - dna_length
                           : current_position;
        transcription_length += 3;
    }
}


/**
 * Compute the pseudo-chimical model (i.e. the width, height and location in the phenotypic space) of a genes/protein
 *
 * @param w_max : Maximum width of the triangle generated by a Protein
 */
void Organism::translate_protein(double w_max) {
    const int dna_length = dna_.seq_.size();

    for (Protein& protein : proteins) {
        if (!protein.is_init_) continue;

        int c_pos = protein.protein_start;
        c_pos += 13;
        c_pos = c_pos >= dna_length
                ? c_pos - dna_length
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
            c_pos = c_pos >= dna_length
                    ? c_pos - dna_length
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

void Organism::compute_phenotype_fitness(double selection_pressure, double* target) {
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


