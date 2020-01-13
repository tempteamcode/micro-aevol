//
// Created by arrouan on 01/10/18.
//

#pragma once

#include <cstring>
#include <vector>
#include <zlib.h>

#include "Threefry.h"

constexpr int CODON_SIZE = 3;

constexpr const char* PROM_SEQ = "0101011001110010010110";
constexpr const char* SHINE_DAL_SEQ = "011011000";
constexpr const char* PROTEIN_END = "001"; // CODON_STOP

class Dna {
public:
  Dna() = default;

  Dna(const Dna& other) = default;

  inline Dna(int length) : seq_(length) { }

  // inline Dna(int length, char* genome) : seq_(length) { strcpy(seq_.data(), genome); }

  Dna(int length, Threefry::Gen&& rng);

  ~Dna() = default;

  inline int length() const { return seq_.size(); };

  void save(gzFile backup_file);
  void load(gzFile backup_file);

  inline void do_switch(int pos) {
    if (seq_[pos] == '0')
      seq_[pos] = '1';
    else
      seq_[pos] = '0';
  }

  int promoter_at(int pos);

  int terminator_at(int pos);

  bool shine_dal_start(int pos);

  bool protein_stop(int pos);

  int codon_at(int pos);

public:
  std::vector<char> seq_; // accessed in Algorithms.cu
};

