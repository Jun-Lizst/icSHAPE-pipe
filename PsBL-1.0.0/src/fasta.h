

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

    Fasta(const string &fasta_file_name, bool strict = false);
    
    StringArray keys() const;
    MapStringuLONG lens() const;
    uLONG len(const string &key) const;
    uLONG size() const;
    const string &seq(const string &key) const;
    const string &annotation(const string &key) const;
    string seq_frag(const string &key, uLONG start, uLONG len=string::npos, STRAND strand=POSITIVE) const;
    bool has(const string &key) const { return sequence.find(key) != sequence.cend(); }
    //bool has_annotation(const string &key) const { return sequence.find(key) != sequence.cend(); }

    MapStringString::const_iterator cbegin() const { return sequence.cbegin(); }
    MapStringString::const_iterator cend() const { return sequence.cend(); }
    MapStringString::iterator begin(){ return sequence.begin(); }
    MapStringString::iterator end(){ return sequence.end(); }

    /*
    static void test_this_class(const string &file_name)
    {
        using namespace std;

        Fasta fasta(file_name);
        MapStringuLONG seq_lens = fasta.lens();
        StringArray seq_keys = fasta.keys();
        cout << "Total " << fasta.size() << endl;
        for(auto this_key: seq_keys)
        {
            cout << this_key << "\t" << seq_lens.at(this_key) << "\t" << fasta.annotation(this_key) << "\t";
            cout << fasta.seq_frag(this_key, 0, 10) << "\n";
        }
    }
    */

private:
    StringArray chr_ids; // keep the raw order of chr_ids
    MapStringString sequence;
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

