
#include "fasta.h"

using namespace std;

namespace pan{

Fasta::Fasta(const string &fasta_file_name, bool strict)
{
    ifstream SEQ(fasta_file_name, ifstream::in);
    if(not SEQ){
        throw runtime_error( "Bad_Input_File: "+fasta_file_name );
    }

    string this_line; 
    string cur_seq;
    while(getline(SEQ, this_line))
    {
        if(this_line.empty()){ continue; }
        if(this_line[0] == '>')
        {
            istringstream head_line(this_line);
            
            head_line >> cur_seq;
            cur_seq = cur_seq.substr(1);
            auto pos = sequence.find(cur_seq);
            if( pos != sequence.end() )
            {
                if(strict)
                    throw Unexpected_Error("Error: duplicate Fatsa Sequence: "+cur_seq);
                else
                {
                    cerr << "Warning: " << "duplicate fatsa sequence: " << cur_seq << "; only preserve the last one" << endl;
                    sequence[cur_seq] = "";
                    sequence_annotation[cur_seq] = "";
                    continue;
                }
            }
            sequence[cur_seq];

            string annotation;
            while( head_line.good() )
            {
                string cur_annotation;
                head_line >> cur_annotation;
                annotation.append( (annotation.empty() ? "" : " ")+cur_annotation );
            }
            sequence_annotation[cur_seq] = annotation;
            chr_ids.push_back(cur_seq);

        }else{
            sequence[cur_seq].append(this_line);
        }
    }

    SEQ.close();
}

StringArray Fasta::keys() const
{
    return chr_ids;
}

uLONG Fasta::len(const string &key) const
{
    return sequence.at(key).size();
}

MapStringuLONG Fasta::lens() const
{
    MapStringuLONG len;
    for(auto iter=sequence.cbegin(); iter!=sequence.cend(); iter++)
        len[iter->first] = iter->second.size();
    return len;
}

uLONG Fasta::size() const
{
    return sequence.size();
}

const string &Fasta::seq(const string &key) const
{
    return sequence.at(key);
}

const string &Fasta::annotation(const string &key) const
{
    return sequence_annotation.at(key);
}

string Fasta::seq_frag(const string &key, uLONG start, uLONG len, STRAND strand) const
{
    if(strand == Fasta::POSITIVE)
        return sequence.at(key).substr(start, len);
    else{
        return reverse_comp(sequence.at(key).substr(start, len));
    }
}

void load_amino_acid(MapStringT<Amino_Acid> &aa_vec, const string &file_name)
{
    ifstream IN(file_name, ifstream::in);
    if(not IN)
    {
        throw runtime_error( "Bad_Input_File: "+file_name );
    }

    aa_vec.clear();
    string cur_line;
    while(getline(IN, cur_line))
    {
        if(cur_line.empty())
            continue;
        istringstream line_in(cur_line);
        Amino_Acid aa;
        line_in >>  aa.animo_acid >> aa.three_letter >> aa.one_letter >> 
                    aa.side_chain_name >> aa.side_chain_polarity >> aa.side_chain_charge >>
                    aa.hydropathy_index;
        aa_vec[ aa.three_letter ] = aa;
    }
    IN.close();
}

void load_codon_map(MapStringString &codon_map, const string &file_name)
{
    ifstream IN(file_name, ifstream::in);
    if(not IN)
    {
        throw runtime_error( "Bad_Input_File: "+file_name );
    }

    codon_map.clear();
    string cur_line;
    while(getline(IN, cur_line))
    {
        if(cur_line.empty())
            continue;

        auto aa_codon = split(cur_line, '=');
        if(aa_codon.size() != 2)
        {
            throw runtime_error("Bad Codon_table " + file_name);
        }
        auto codons = split(aa_codon.at(1), '|');
        for(string codon: codons)
        {
            codon_map[ codon ] = aa_codon.at(0);
            std::replace(codon.begin(), codon.end(), 'U', 'T');
            codon_map[ codon ] = aa_codon.at(0);
        }
        
    }
    IN.close();
}

ostream& operator<<(ostream& OUT, const Amino_Acid &aa)
{
    OUT << aa.animo_acid << "\t" << aa.three_letter << "\t" << aa.one_letter
        << "\t" << aa.side_chain_name << "\t" << aa.side_chain_polarity 
        << "\t" << aa.side_chain_charge << "\t" << aa.hydropathy_index;

    return OUT;
}

string reverse_comp(const string &raw_seq)
{
    string rev_comp_seq;
    for(auto iter=raw_seq.crbegin(); iter!=raw_seq.crend(); iter++)
    {
        switch(*iter)
        {
            case 'A': case 'a':
                rev_comp_seq += 'T';
                break;
            case 'T': case 'U': case 't': case 'u':
                rev_comp_seq += 'A';
                break;
            case 'C': case 'c':
                rev_comp_seq += 'G';
                break;
            case 'G': case 'g':
                rev_comp_seq += 'C';
                break;
            default:
                rev_comp_seq += 'N';
        }

        /*
        if(*iter=='A' or *iter=='a')
            rev_comp_seq += 'T';
        else if(*iter=='T' or *iter=='U' or )
            rev_comp_seq += 'A';
        else if(*iter=='C')
            rev_comp_seq += 'G';
        else if(*iter=='G')
            rev_comp_seq += 'C';
        else
            rev_comp_seq += 'N';
        */
    }
    return rev_comp_seq;
}

StringArray translation(const string &raw_seq, const MapStringString &codon_map)
{
    StringArray protein;
    for(size_t idx=0; idx<raw_seq.size()-2; idx+=3)
    {
        protein.push_back( codon_map.at(raw_seq.substr(idx, 3)) );
    }
    return protein;
}

string cut_seq(const string &raw_seq, size_t seg_length)
{
    if(seg_length >= raw_seq.size())
        return raw_seq+"\n";

    string cutted_seq;
    size_t idx=0;
    for(; idx<=raw_seq.size()-seg_length; idx+=seg_length)
    {
        cutted_seq += raw_seq.substr(idx, seg_length);
        cutted_seq += "\n";
    }
    if(idx!=raw_seq.size())
    {
        cutted_seq += raw_seq.substr(idx);
        cutted_seq += "\n";
    }
    return cutted_seq;
}

int global_align(const string &a,
        const string &b,
        string &a_aligned, 
        string &b_aligned,
        const int gap_penalty,
        const int mismatch_penalty)
{
    a_aligned.clear();
    b_aligned.clear();

    struct DP_node
    {
        enum DIRECT{ left, top, left_top };
        int penalty;
        DIRECT direct;
        
        DP_node(int p=0, DIRECT d=left_top): penalty(p), direct(d){}
    };

    size_t n = a.size();
    size_t m = b.size();

    vector<vector<DP_node>> A(n + 1, vector<DP_node>(m + 1));

    for (size_t i = 0; i <= m; ++i)
        A[0][i] = DP_node(gap_penalty * i, DP_node::top);
    for (size_t i = 0; i <= n; ++i)
        A[i][0] = DP_node(gap_penalty * i, DP_node::left);

    for (size_t i = 1; i <= n; ++i)
    {
        for (size_t j = 1; j <= m; ++j)
        {
            char x_i = a[i-1];
            char y_j = b[j-1];

            int penalty_lt = A[i-1][j-1].penalty + ( (x_i == y_j) ? 0 : mismatch_penalty );
            int penalty_l = A[i-1][j].penalty + gap_penalty;
            int penalty_t = A[i][j-1].penalty + gap_penalty;

            if(penalty_lt <= penalty_l and penalty_lt <= penalty_t)
            {
                A[i][j] = DP_node(penalty_lt, DP_node::left_top);
            }else if(penalty_l <= penalty_lt and penalty_l <= penalty_t)
            {
                A[i][j] = DP_node(penalty_l, DP_node::left);
            }else if(penalty_t <= penalty_lt and penalty_t <= penalty_l)
            {
                A[i][j] = DP_node(penalty_t, DP_node::top);
            }else{
                // no thing
            }
        }
    }

    size_t j = m;
    size_t i = n;
    int a_count = 0, b_count = 0;
    for (; i >= 1 and j >= 1; )
    {
        char x_i = a[i-1];
        char y_j = b[j-1];
        if (A[i][j].direct == DP_node::left_top)
        {
            a_aligned.push_back(x_i);
            b_aligned.push_back(y_j);
            --j; --i;
            a_count++; b_count++;
        }
        else if (A[i][j].direct == DP_node::left)
        {
            a_aligned.push_back(x_i);
            b_aligned.push_back('-');
            --i;
            a_count++;
        }
        else if(A[i][j].direct == DP_node::top) 
        {
            a_aligned.push_back('-');
            b_aligned.push_back(y_j);
            --j;
            b_count++;
        }else{
            cerr << "Error" << endl;
        }
    }

    while (i >= 1 && j < 1)
    {
        a_aligned.push_back(a[i-1]);
        b_aligned.push_back('-');
        --i;
        a_count++;
    }
    while (j >= 1 && i < 1)
    {
        a_aligned.push_back('-');
        b_aligned.push_back(b[j-1]);
        --j;
        b_count++;
    }

    reverse(a_aligned.begin(), a_aligned.end());
    reverse(b_aligned.begin(), b_aligned.end());

    return A[n][m].penalty;
}

}
