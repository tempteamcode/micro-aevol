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
  inline Dna() : len_(0), seq_(nullptr) { }

  inline Dna(const Dna& other) : len_(other.len_), seq_(other.len_ > 0 ? new char[other.len_] : nullptr) { if (len_ > 0) std::memcpy(seq_, other.seq_, len_); }

  inline Dna(Dna&& other) : len_(other.len_), seq_(other.seq_) { other.len_ = 0; other.seq_ = nullptr; }

  Dna& operator=(const Dna& other)
  {
    if (&other == this) return *this;
    if (len_ > 0) delete[] seq_;
    if ((len_ = other.len_) > 0)
    {
      seq_ = new char[len_];
      std::memcpy(seq_, other.seq_, len_);
    }
    else
    {
      seq_ = nullptr;
    }
    return *this;
  }

  Dna& operator=(Dna&& other)
  {
    if (&other == this) return *this;
    if (len_ > 0) delete[] seq_;
    len_ = other.len_; other.len_ = 0;
    seq_ = other.seq_; other.seq_ = nullptr;
    return *this;
  }

  inline Dna(int length) : len_(length), seq_(new char[length]) { }

  // inline Dna(int length, char* sequence) : len_(length), seq_(new char[length]) { std::memcpy(seq_, sequence, length); }

  inline Dna(int length, char* sequence) : len_(length), seq_(sequence) { }

  Dna(int length, Threefry::Gen&& rng);

  inline ~Dna() { if (len_ > 0) delete[] seq_; }

  inline int length() const { return len_; }
  inline const char* data() const { return seq_; }

  void save(gzFile backup_file) const;

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

private:
  int len_;
  char* seq_;
};

inline Dna Dna_load(gzFile backup_file) {
  int length;
  gzread(backup_file, &length, sizeof(length));

  char* seq = new char[length];
  gzread(backup_file, seq, length * sizeof(seq[0]));

  return Dna(length, std::move(seq));
}

