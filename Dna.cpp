//
// Created by arrouan on 01/10/18.
//

#include "Dna.h"
#include "ExpManager.h"

Dna::Dna(const Dna& clone) : seq_(clone.seq_) {
}

// Generate a random genome
Dna::Dna(int length, Threefry::Gen& rng) : seq_(length, [&rng] () { return (rng.random(NB_BASE) != 0); })
{
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

void Dna::do_switch(int pos) {
  seq_.flip(pos);
}

int Dna::promoter_at(int pos) {
  return seq_.getHammingDistance(pos, PROM_SEQ_LEN, PROM_SEQ);
}

int Dna::terminator_at(int pos) {
    unsigned int subseq = seq_.getSequence(pos+10-4+1, 4);
    unsigned int temp = (lookuptable[ subseq & 0xff ]<<24 |
                         lookuptable[ (subseq >> 8) & 0xff ]<<16 |
                         lookuptable[ (subseq >> 16 )& 0xff ]<< 8 |
                         lookuptable[ (subseq >>24 ) & 0xff ])
            >> (sizeof(int)*8-4);

    return seq_.getHammingDistance(pos, 4, temp);
}

bool Dna::shine_dal_start(int pos) {
  unsigned int first_part = (SHINE_DAL_SEQ&0b111111000)>>3;
  unsigned int second_part = (SHINE_DAL_SEQ&0b000000111);

  bool start = seq_.search(pos, 6, first_part);

  if(start)
      start = seq_.search(pos+10, 3, second_part);

  return start;
}

bool Dna::protein_stop(int pos) {
  return seq_.search(pos, PROTEIN_END_LEN, PROTEIN_END);
}

int Dna::codon_at(int pos) {
  return seq_.getSequence(pos, CODON_SIZE);
}
