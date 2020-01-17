//
// Created by arrouan on 01/10/18.
//

#pragma once

#include <cstdio>
#include <cstring>
#include <cassert>
#include <cstdint>
#include <vector>
#include <zlib.h>

#include "Threefry.h"
#include "bitset.h"


typedef uint32_t int_t;

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

constexpr int8_t CODON_LEN = 3;


class ExpManager;

class Dna {

 public:
  Dna(int length, Threefry::Gen& rng);

  Dna(int length, char* genome);

  Dna(int length);

  ~Dna() = default;

  int length() const;

  void save(gzFile backup_file);
  friend Dna Dna_load(gzFile backup_file);

  void do_switch(int pos);

  int promoter_at(int pos);

  int terminator_at(int pos);

  bool shine_dal_start(int pos);

  bool protein_stop(int pos);

  int codon_at(int pos);

  own_bitset seq_;
};

inline void Dna::save(gzFile backup_file) {
    int dna_length = length();
    std::string data = seq_.export_string();
    gzwrite(backup_file, &dna_length, sizeof(dna_length));
    gzwrite(backup_file, data.c_str(), dna_length * sizeof(char));
}

inline Dna Dna_load(gzFile backup_file) {
    int dna_length;
    gzread(backup_file, &dna_length, sizeof(dna_length));

    char tmp_seq[dna_length];
    gzread(backup_file, tmp_seq, dna_length * sizeof(tmp_seq[0]));

    return Dna(dna_length, tmp_seq);
}

