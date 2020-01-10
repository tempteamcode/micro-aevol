//
// Created by arrouan on 01/10/18.
//

#include "Dna.h"
#include "ExpManager.h"

#define PROM_SEQ_BIT(n) ((PROM_SEQ>>(PROM_SEQ_LEN-n-1)) & 0b1)
#define SHINE_DAL_SEQ_BIT(n) ((SHINE_DAL_SEQ>>(SHINE_DAL_SEQ_LEN-n-1)) & 0b1)
#define PROTEIN_END_BIT(n) ((PROTEIN_END>>(PROTEIN_END_LEN-n-1)) & 0b1)

Dna::Dna(const Dna& clone) : seq_(clone.seq_) {
}

Dna::Dna(int length, Threefry::Gen& rng) : seq_(length) {
  // Generate a random genome
  seq_.generate(length, [&rng] () { return (rng.random(NB_BASE) != 0); });
}

Dna::Dna(char* genome, int length) : seq_(length) {
  seq_.import_string(genome, length);
}

Dna::Dna(int length) : seq_(length) {
  seq_.set_all(false);
}

int Dna::length() const {
  return seq_.size();
}

void Dna::save(gzFile backup_file) {
    int dna_length = length();
    std::string data = seq_.export_string();
    gzwrite(backup_file, &dna_length, sizeof(dna_length));
    gzwrite(backup_file, data.c_str(), dna_length * sizeof(char));
}

void Dna::load(gzFile backup_file) {
    int dna_length;
    gzread(backup_file, &dna_length, sizeof(dna_length));

    char tmp_seq[dna_length];
    gzread(backup_file, tmp_seq, dna_length * sizeof(tmp_seq[0]));

    seq_.import_string(tmp_seq, dna_length);
}

/**
 * Remove the DNA inbetween pos_1 and pos_2
 *
 * @param pos_1
 * @param pos_2
 */
void Dna::remove(int pos_1, int pos_2) {
  assert(pos_1 >= 0 && pos_2 >= pos_1 && pos_2 <= seq_.size());
  seq_.range_erase(pos_1, pos_2);
}

/**
 * Insert a sequence of a given length at a given position into the DNA of the Organism
 *
 * @param pos : where to insert the sequence
 * @param seq : the sequence itself
 * @param seq_length : the size of the sequence
 */
void Dna::insert(int pos, const own_dynamic_bitset& seq) {
// Insert sequence 'seq' at position 'pos'
  assert(pos >= 0 && pos < seq_.size());

  seq_.range_insert(pos, seq);
}

/**
 * Insert a sequence of a given length at a given position into the DNA of the Organism
 *
 * @param pos : where to insert the sequence
 * @param seq : the sequence itself
 * @param seq_length : the size of the sequence
 */
void Dna::insert(int pos, Dna* seq) {
  insert(pos, seq->seq_);
}

void Dna::do_switch(int pos) {
  seq_.flip(pos);
}

void Dna::do_duplication(int pos_1, int pos_2, int pos_3) {
  // Duplicate segment [pos_1; pos_2[ and insert the duplicate before pos_3
  char* duplicate_segment = NULL;

  int32_t seg_length;

  if (pos_1 < pos_2) {
    //
    //       pos_1         pos_2                   -> 0-
    //         |             |                   -       -
    // 0--------------------------------->      -         -
    //         ===============                  -         - pos_1
    //           tmp (copy)                      -       -
    //                                             -----      |
    //                                             pos_2    <-'
    //
    own_dynamic_bitset seq_dupl(seq_, pos_1, pos_2);

    insert(pos_3, seq_dupl);
  } else { // if (pos_1 >= pos_2)
    // The segment to duplicate includes the origin of replication.
    // The copying process will be done in two steps.
    //
    //                                            ,->
    //    pos_2                 pos_1            |      -> 0-
    //      |                     |                   -       - pos_2
    // 0--------------------------------->     pos_1 -         -
    // ======                     =======            -         -
    //  tmp2                        tmp1              -       -
    //                                                  -----
    //
    //
    own_dynamic_bitset seq_dupl_p1(seq_, pos_1, seq_.size());
    own_dynamic_bitset seq_dupl_p2(seq_, 0, pos_2);

    insert(pos_3, seq_dupl_p2);
    insert(pos_3, seq_dupl_p1);
  }
}

int Dna::promoter_at(int pos) {
  int prom_dist[22];

  for (int motif_id = 0; motif_id < 22; motif_id++) {
    // Searching for the promoter
    prom_dist[motif_id] =
        PROM_SEQ_BIT(motif_id) ==
        seq_[
            pos + motif_id >= seq_.size() ? pos +
                                            motif_id -
                                            seq_.size()
                                          : pos +
                                            motif_id]
        ? 0
        : 1;

  }


  // Computing if a promoter exists at that position
  int dist_lead = prom_dist[0] +
                  prom_dist[1] +
                  prom_dist[2] +
                  prom_dist[3] +
                  prom_dist[4] +
                  prom_dist[5] +
                  prom_dist[6] +
                  prom_dist[7] +
                  prom_dist[8] +
                  prom_dist[9] +
                  prom_dist[10] +
                  prom_dist[11] +
                  prom_dist[12] +
                  prom_dist[13] +
                  prom_dist[14] +
                  prom_dist[15] +
                  prom_dist[16] +
                  prom_dist[17] +
                  prom_dist[18] +
                  prom_dist[19] +
                  prom_dist[20] +
                  prom_dist[21];

  return dist_lead;
}

int Dna::terminator_at(int pos) {
  int term_dist[4];
  for (int motif_id = 0; motif_id < 4; motif_id++) {

    // Search for the terminators
    term_dist[motif_id] =
        seq_[
            pos + motif_id >= seq_.size() ? pos +
                                            motif_id -
                                            seq_.size() :
            pos + motif_id] !=
        seq_[
            pos - motif_id + 10 >= seq_.size() ?
            pos - motif_id + 10 - seq_.size() :
            pos -
            motif_id +
            10] ? 1
                : 0;
  }
  int dist_term_lead = term_dist[0] +
                       term_dist[1] +
                       term_dist[2] +
                       term_dist[3];

  return dist_term_lead;
}

bool Dna::shine_dal_start(int pos) {
  bool start = false;
  int t_pos, k_t;

  for (int k = 0; k < 9; k++) {
    k_t = k >= 6 ? k + 4 : k;
    t_pos = pos + k_t >= seq_.size() ? pos + k_t -
                                       seq_.size()
                                     : pos + k_t;

    if (seq_[t_pos] ==
        SHINE_DAL_SEQ_BIT(k)) {
      start = true;
    } else {
      start = false;
      break;
    }
  }

  return start;
}

bool Dna::protein_stop(int pos) {
  bool is_protein;
  int t_k;

  for (int k = 0; k < 3; k++) {
    t_k = pos + k >= seq_.size() ?
          pos - seq_.size() + k :
          pos + k;

    if (seq_[t_k] ==
        PROTEIN_END_BIT(k)) {
      is_protein = true;
    } else {
      is_protein = false;
      break;
    }
  }

  return is_protein;
}

int Dna::codon_at(int pos) {
  int value = 0;

  int t_pos;

  for (int i = 0; i < 3; i++) {
    t_pos =
        pos + i >= seq_.size() ? pos + i -
                                 seq_.size()
                               : pos + i;
    if (seq_[t_pos] ==
        '1')
      value += 1 << (CODON_SIZE - i - 1);
  }

  return value;
}
