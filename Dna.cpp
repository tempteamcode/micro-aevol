//
// Created by arrouan on 01/10/18.
//

#include "Dna.h"
#include "ExpManager.h"

Dna::Dna(int length, Threefry::Gen&& rng)
: len_(length)
, seq_(new char[length])
{
  // Generate a random genome
  for (int i = 0; i < length; i++) {
    seq_[i] = '0' + rng.random(NB_BASE);
  }
}

void Dna::save(gzFile backup_file) const {
    int len = len_;
    gzwrite(backup_file, &len, sizeof(len));
    gzwrite(backup_file, seq_, len * sizeof(seq_[0]));
}

int Dna::promoter_at(int pos) {
  int dist_lead = 0;

  for (int motif_id = 0; motif_id < 22; motif_id++) {
    // Searching for the promoter
    if (PROM_SEQ[motif_id] !=
        seq_[
            pos + motif_id >= len_
            ? pos + motif_id - len_
            : pos + motif_id
        ]
       )
      dist_lead++;
  }

  return dist_lead;
}

int Dna::terminator_at(int pos) {
  int dist_term_lead = 0;

  for (int motif_id = 0; motif_id < 4; motif_id++) {
    // Search for the terminators
    if(
        seq_[
            pos + motif_id >= len_
            ? pos + motif_id - len_
            : pos + motif_id
            ]
        !=
        seq_[
            pos - motif_id + 10 >= len_
            ? pos - motif_id + 10 - len_
            : pos - motif_id + 10
            ]
    )
      dist_term_lead++;
  }

  return dist_term_lead;
}

bool Dna::shine_dal_start(int pos) {
  int t_pos, k_t;

  for (int k = 0; k < 9; k++) {
    k_t = k >= 6 ? k + 4 : k;
    t_pos = pos + k_t >= len_
          ? pos + k_t - len_
          : pos + k_t;

    if (seq_[t_pos] != SHINE_DAL_SEQ[k]) {
      return false;
    }
  }

  return true;
}

bool Dna::protein_stop(int pos) {
  int t_k;

  for (int k = 0; k < 3; k++) {
    t_k = pos + k >= len_
          ? pos - len_ + k
          : pos + k;

    if (seq_[t_k] != PROTEIN_END[k]) {
      return false;
    }
  }

  return true;
}

int Dna::codon_at(int pos) {
  int value = 0;

  for (int i = 0; i < 3; i++) {
    int t_pos =
        pos + i >= len_
        ? pos + i - len_
        : pos + i;
    if (seq_[t_pos] == '1')
      value += 1 << (CODON_SIZE - i - 1);
  }

  return value;
}

