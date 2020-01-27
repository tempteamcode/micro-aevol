//
// Created by arrouan on 01/10/18.
//

#pragma once

#include <cstring>
#include <vector>
#include <zlib.h>

#include "Threefry.h"
#include "bitset.h"
#include "bits.h"

typedef uint32_t int_t;

constexpr const int NB_BASE = 2;

constexpr const int_t PROM_SEQ     = 0b0101011001110010010110;
constexpr const int_t PROM_SEQ_REV = 0b0110100100111001101010;
constexpr const int PROM_SEQ_LEN = 22;

constexpr const int_t SHINE_DAL_SEQ_BITS = 0b0110110000000;
constexpr const int_t SHINE_DAL_SEQ_MASK = 0b1111110000111;
constexpr const int_t SHINE_DAL_SEQ_BITS_REV = 0b0000000110110;
constexpr const int_t SHINE_DAL_SEQ_MASK_REV = 0b1110000111111;
constexpr const int SHINE_DAL_SEQ_LEN = 13;

constexpr const int_t PROTEIN_END     = 0b001; // CODON_STOP
constexpr const int_t PROTEIN_END_REV = 0b100;
constexpr const int PROTEIN_END_LEN = 3;

constexpr const int CODON_LEN = 3;


class Dna {
public:
  // inline Dna() : seq_(0) { } // used in Organism.cpp

  inline Dna(int length, Threefry::Gen& rng) : seq_(length, [&rng] () { return (rng.random(NB_BASE) != 0); }) { }

  inline Dna(int length, char* genome) : seq_(genome, length) { }

  inline Dna(int length) : seq_(length) { seq_.set_all(false); }

  void save(gzFile backup_file) const;
  friend Dna&& Dna_load(gzFile backup_file);

  inline void do_switch(int pos) { seq_.flip(pos); }

  inline int promoter_at(int pos) { return hamming_distance(seq_.getSequenceRev(pos, PROM_SEQ_LEN), PROM_SEQ_REV); }

  inline int terminator_at(int pos) { return hamming_distance(seq_.getSequenceRev(pos, 4), seq_.getSequence(pos+10-4+1, 4)); }

  inline bool shine_dal_start(int pos) { return (seq_.getSequenceRev(pos, SHINE_DAL_SEQ_LEN) & SHINE_DAL_SEQ_MASK_REV) == SHINE_DAL_SEQ_BITS_REV; }

  inline bool protein_stop(int pos) { return seq_.getSequenceRev(pos, PROTEIN_END_LEN) == PROTEIN_END_REV; }

  inline int codon_at(int pos) { return seq_.getSequence(pos, CODON_LEN); }

public:
  own_bitset seq_; // used in Algorithms.cu, Organism.cpp
};

inline void Dna::save(gzFile backup_file) const {
    int dna_length = seq_.size();
    std::string data = seq_.export_string();
    gzwrite(backup_file, &dna_length, sizeof(dna_length));
    gzwrite(backup_file, data.c_str(), dna_length * sizeof(char));
}

inline Dna&& Dna_load(gzFile backup_file) {
    int dna_length;
    gzread(backup_file, &dna_length, sizeof(dna_length));

    char tmp_seq[dna_length];
    gzread(backup_file, tmp_seq, dna_length * sizeof(tmp_seq[0]));

    return std::move(Dna(dna_length, tmp_seq));
}

