

#ifndef FASTA_H
#define FASTA_H

#include "pan_type.h"
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <sstream>
#include <string>
#include <exception>
#include "string_split.h"
#include <cstring>
#include <algorithm>
#include "exceptions.h"

namespace pan{

using std::vector;
using std::unordered_map;
using std::string;

class Fasta
{
public:
    enum STRAND{ POSITIVE, NEGATIVE };

    /* 
    Constract a Fasta sequence
    fasta_file_name     -- Input fasta file name
    strict              -- If true, raise Unexpected_Error when duplicate sequences are found
    */
    Fasta(const string &fasta_file_name, bool strict = false);
    
    /*
    Return the raw chrID list
    */
    StringArray keys() const;

    /*
    Return a map of raw chrID => chrLen
    */
    MapStringuLONG lens() const;
    
    /*
    Return the length of a given chrID
    key                 -- The chrID
    */
    uLONG len(const string &key) const;

    /*
    Return the number of sequence
    */
    uLONG size() const;

    /*
    Return the sequence of given chrID
    key                 -- The chrID
    */
    const string &seq(const string &key) const;

    /*
    Return the annotation of given chrID
    key                 -- The chrID
    */
    const string &annotation(const string &key) const;

    /*
    Return a fragment of the sequence of given chrID
    key                 -- The chrID
    start               -- The start position, 0-based
    len                 -- The sub-sequence length
    strand              -- Strand
    */
    string seq_frag(const string &key, uLONG start, uLONG len=string::npos, STRAND strand=POSITIVE) const;

    /*
    Test if a chrID is in the Fasta
    key                 -- The chrID
    */
    bool has(const string &key) const { return sequence.find(key) != sequence.cend(); }
    
    // Iterators
    MapStringString::const_iterator cbegin() const { return sequence.cbegin(); }
    MapStringString::const_iterator cend() const { return sequence.cend(); }
    MapStringString::iterator begin(){ return sequence.begin(); }
    MapStringString::iterator end(){ return sequence.end(); }

private:

    // The chrID list, keep the raw order
    StringArray chr_ids;
    // A map of chrID => sequence 
    MapStringString sequence;
    // A map of chrID => annotation
    MapStringString sequence_annotation;
};


struct Amino_Acid
{
    string animo_acid;
    string three_letter;
    string one_letter;
    
    string side_chain_name;
    string side_chain_polarity;
    string side_chain_charge;
    double hydropathy_index;
};

ostream& operator<<(ostream& OUT, const Amino_Acid &aa);

void load_amino_acid(MapStringT<Amino_Acid> &aa_vec, 
        const string &file_name="conf_data/amino_acid.txt");

void load_codon_map(MapStringString &codon_map, 
        const string &file_name="conf_data/codon_table.txt");

string reverse_comp(const string &raw_seq);
StringArray translation(const string &raw_seq, const MapStringString &codon_map);

string cut_seq(const string &raw_seq, size_t seg_length=60);

int global_align(const string &a,
        const string &b,
        string &a_aligned, 
        string &b_aligned,
        const int gap_penalty=5,
        const int mismatch_penalty=1);









}
#endif