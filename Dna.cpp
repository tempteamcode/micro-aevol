//
// Created by arrouan on 01/10/18.
//

#include "Dna.h"
#include "ExpManager.h"

#include "bits.h"

// Generate a random genome
Dna::Dna(int length, Threefry::Gen& rng) : seq_(length, [&rng] () { return (rng.random(NB_BASE) != 0); }) { }

Dna::Dna(int length, char* genome) : seq_(genome, length) { }

Dna::Dna(int length) : seq_(length) {
  seq_.set_all(false);
}

int Dna::length() const {
  return seq_.size();
}

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
    const int_t first_part_rev = (SHINE_DAL_SEQ_REV&0b000111111);
    const int_t second_part_rev = (SHINE_DAL_SEQ_REV&0b111000000)>>6;

    return seq_.getSequenceRev(pos, 6) == first_part_rev &&
           seq_.getSequenceRev(pos+10, 3) == second_part_rev;
}

bool Dna::protein_stop(int pos) {
    return seq_.getSequenceRev(pos, PROTEIN_END_LEN) == PROTEIN_END_REV;
}

int Dna::codon_at(int pos) {
    return seq_.getSequence(pos, CODON_SIZE);
}

