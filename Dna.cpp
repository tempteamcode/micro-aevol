//
// Created by arrouan on 01/10/18.
//

#include "Dna.h"
#include "ExpManager.h"

Dna::Dna(int length, Threefry::Gen&& rng) : seq_(length) {
  // Generate a random genome
  for (int32_t i = 0; i < length; i++) {
    seq_[i] = '0' + rng.random(NB_BASE);
  }
}

void Dna::save(gzFile backup_file) const {
    int dna_length = length();
    gzwrite(backup_file, &dna_length, sizeof(dna_length));
    gzwrite(backup_file, seq_.data(), dna_length * sizeof(seq_[0]));
}

void Dna::load(gzFile backup_file) {
    int dna_length;
    gzread(backup_file, &dna_length, sizeof(dna_length));

    char tmp_seq[dna_length];
    gzread(backup_file, tmp_seq, dna_length * sizeof(tmp_seq[0]));

    seq_ = std::vector<char>(tmp_seq, tmp_seq+dna_length);
}

int Dna::promoter_at(int pos) {
  int dist_lead = 0;

  for (int motif_id = 0; motif_id < 22; motif_id++) {
    // Searching for the promoter
    if (PROM_SEQ[motif_id] !=
        seq_[
            pos + motif_id >= seq_.size()
            ? pos + motif_id - seq_.size()
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
            pos + motif_id >= seq_.size()
            ? pos + motif_id - seq_.size()
            : pos + motif_id
            ]
        !=
        seq_[
            pos - motif_id + 10 >= seq_.size()
            ? pos - motif_id + 10 - seq_.size()
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
    t_pos = pos + k_t >= seq_.size()
          ? pos + k_t - seq_.size()
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
    t_k = pos + k >= seq_.size()
          ? pos - seq_.size() + k
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
        pos + i >= seq_.size()
        ? pos + i - seq_.size()
        : pos + i;
    if (seq_[t_pos] == '1')
      value += 1 << (CODON_SIZE - i - 1);
  }

  return value;
}

