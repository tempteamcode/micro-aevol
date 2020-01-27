//
// Created by arrouan on 01/10/18.
//

#include "Dna.h"
#include "ExpManager.h"
#include "bits.h"

// Generate a random genome
Dna::Dna(int length, Threefry::Gen& rng) : seq_(length, [&rng] () { return (rng.random(NB_BASE) != 0); }) { }

Dna::Dna(int length, char* genome) : seq_(genome, length) { }

/*
Dna::Dna(int length) : seq_(length) {
  seq_.set_all(false);
}
*/

void Dna::do_switch(int pos) {
  seq_.flip(pos);
}

int Dna::promoter_at(int pos) {
    return hamming_distance(seq_.getSequenceRev(pos, PROM_SEQ_LEN), PROM_SEQ_REV);
}

int Dna::terminator_at(int pos) {
    int_t subseq_rev = seq_.getSequence(pos+10-4+1, 4);
    return hamming_distance(seq_.getSequenceRev(pos, 4), subseq_rev);
}

bool Dna::shine_dal_start(int pos) {
    return (seq_.getSequenceRev(pos, SHINE_DAL_SEQ_LEN) & SHINE_DAL_SEQ_MASK_REV) == SHINE_DAL_SEQ_BITS_REV;
}

bool Dna::protein_stop(int pos) {
    return seq_.getSequenceRev(pos, PROTEIN_END_LEN) == PROTEIN_END_REV;
}

int Dna::codon_at(int pos) {
    return seq_.getSequence(pos, CODON_LEN);
}

