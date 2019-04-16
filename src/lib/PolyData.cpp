#include "PolyData.hpp"
#include <glob.h>
#include <fstream>
#include "Exception.hpp"

using namespace std;

template <class T>
bool item_in_vector(vector<T> const &v, T item) {
    return find(v.begin(), v.end(), item) != v.end();
}

vector<string> GlobVector(const string &pattern) {
    glob_t glob_result;
    glob(pattern.c_str(), GLOB_TILDE, NULL, &glob_result);
    vector<string> files;
    for (unsigned int i = 0; i < glob_result.gl_pathc; ++i) {
        files.emplace_back(string(glob_result.gl_pathv[i]));
    }
    globfree(&glob_result);
    return files;
}

PolyData::PolyData(CodonSequenceAlignment *from_alignment, string const &ali_path)
    : SampleSize(from_alignment->GetNtaxa(), 1) {
    Alignment = from_alignment;

    // PolyData will search for every files (.vcf) that are in the same
    // directory (folder) as ali_path. The files (.vcf)  will be matched
    // to taxa (leaves of the tree) if the name of the taxon is included
    // in the file name (e.g. *HomoSapiens*.vcf)

    size_t find_last_sep = ali_path.rfind('.');
    string root_dir = "./";
    if (find_last_sep != string::npos) {
        root_dir = ali_path.substr(0, find_last_sep);
        cerr << "Searching for files (.vcf) in " << root_dir << endl;
    } else {
        cerr << "Searching for files (.vcf) in the current directory" << endl;
    }

    for (int taxon{0}; taxon < Alignment->GetNtaxa(); taxon++) {
        string pattern = root_dir + '*' + Alignment->GetTaxonSet()->GetTaxon(taxon) + "*.vcf";
        vector<string> taxon_files = GlobVector(pattern);

        if (taxon_files.empty()) {
            continue;
        } else if (taxon_files.size() > 1) {
            cerr << "Found " << taxon_files.size() << " files (.vcf) for taxon "
                 << Alignment->GetTaxonSet()->GetTaxon(taxon) << ":" << endl;
            for (auto const &t : taxon_files) { cerr << "\t" << t << endl; }
            cerr << "Only " << *taxon_files.begin() << " will be used" << endl;
        }

        ifstream vcf_stream(*taxon_files.begin());

        Nvcf++;
        string vcf_line;
        vector<string> header_names;

        unsigned sample_size = 0;

        while (getline(vcf_stream, vcf_line)) {
            if (vcf_line.size() > 1 and vcf_line[0] == '#') {
                //

                if (vcf_line[1] != '#') {
                    istringstream vcf_line_stream(vcf_line);
                    string word;

                    while (getline(vcf_line_stream, word, '\t')) { header_names.push_back(word); }

                    // Assert the line contains the information necessary
                    assert(item_in_vector<string>(header_names, "POS"));
                    assert(item_in_vector<string>(header_names, "REF"));
                    assert(item_in_vector<string>(header_names, "ALT"));
                    assert(item_in_vector<string>(header_names, "INFO"));

                } else {
                    string tag = "numberGenotypes";
                    size_t find_sample_size = vcf_line.find(tag);
                    if (find_sample_size != string::npos) {
                        string value = vcf_line.substr(find_sample_size + tag.size() + 1);
                        sample_size = static_cast<unsigned>(stoi(value));
                        SampleSize.at(taxon) = sample_size;
                    }
                }


            } else if (vcf_line.size() > 1 and vcf_line[0] != '#') {
                map<int, unsigned> codon_state_to_count;

                unsigned word_index = 0;
                int nuc_site = -1;
                char ref_nuc = -1;
                char alt_nuc = -1;

                unsigned alt_count = 0;
                string ref_codon{};
                string alt_codon{};

                istringstream vcf_line_stream(vcf_line);
                string word;

                while (getline(vcf_line_stream, word, '\t')) {
                    if (header_names.at(word_index) == "POS") {
                        nuc_site = stoi(word);
                    } else if (header_names.at(word_index) == "REF") {
                        ref_nuc = word[0];
                    } else if (header_names.at(word_index) == "ALT") {
                        alt_nuc = word[0];
                    } else if (header_names.at(word_index) == "INFO") {
                        istringstream word_stream(word);
                        string info_field;
                        while (getline(word_stream, info_field, ';')) {
                            size_t find_equal = info_field.find('=');
                            if (find_equal != string::npos) {
                                if (info_field.substr(0, find_equal) == "ALTCOUNT") {
                                    string value = info_field.substr(find_equal + 1);
                                    alt_count = static_cast<unsigned>(stoi(value));
                                }
                                if (info_field.substr(0, find_equal) == "REFCODON") {
                                    ref_codon = info_field.substr(find_equal + 1);
                                }
                                if (info_field.substr(0, find_equal) == "ALTCODON") {
                                    alt_codon = info_field.substr(find_equal + 1);
                                }
                            }
                        }
                    }
                    word_index++;
                }

                // Assert position has been found
                assert(nuc_site != -1);
                // Assert reference nucleotide has been found
                assert(ref_nuc != -1);
                // Assert alternate nucleotide has been found
                assert(alt_nuc != -1);

                // Assert the alternate count is within the sample size
                // VL removed "alt_count >= 0" because it was always true
                assert(alt_count <= sample_size);
                unsigned ref_count = sample_size - alt_count;

                int site = nuc_site / 3;
                // Assert codon position is not greater than the alignment size
                assert(site < Alignment->GetNsite());

                int ref_state = Alignment->GetState(taxon, site);
                // Assert the reference codon could be found in the alignment
                assert(ref_state >= 0);
                assert(ref_state < Alignment->GetCodonStateSpace()->GetNstate());

                if (not ref_codon.empty()) {
                    string ref_codon_from_ali =
                        Alignment->GetCodonStateSpace()->GetState(ref_state);
                    // Assert that the reference codon in the .vcf is the same as in .ali
                    assert(ref_codon == ref_codon_from_ali);
                }
                ref_codon = Alignment->GetCodonStateSpace()->GetState(ref_state);
                // Assert that the reference nucleotide in the .vcf is the same as in .ali
                assert(ref_codon[nuc_site % 3] == ref_nuc);

                // Assert that the alternate codon in the .vcf is right according to the alternate
                // nucleotide
                if (not alt_codon.empty()) {
                    assert(alt_codon[nuc_site % 3] == alt_nuc);
                    if ((alt_codon[(nuc_site + 1) % 3] != ref_codon[(nuc_site + 1) % 3]) or
                        (alt_codon[(nuc_site + 2) % 3] != ref_codon[(nuc_site + 2) % 3])) {
                        cerr << "Double mutation observed at site " << site << " for taxon "
                             << Alignment->GetTaxonSet()->GetTaxon(taxon) << "." << endl;
                        alt_codon[(nuc_site + 1) % 3] = ref_codon[(nuc_site + 1) % 3];
                        alt_codon[(nuc_site + 2) % 3] = ref_codon[(nuc_site + 2) % 3];
                    }
                }
                alt_codon[nuc_site % 3] = alt_nuc;

                int alt_state = Alignment->GetCodonStateSpace()->GetState(alt_codon);
                // Assert the alternate codon is in the genetic code (and also not a stop codon)
                if (alt_state < 0 or alt_state >= Alignment->GetCodonStateSpace()->GetNstate()) {
                    cerr << "The mutation lead to a stop codon at site " << site << " for taxon "
                         << Alignment->GetTaxonSet()->GetTaxon(taxon) << "." << endl;
                    continue;
                }

                codon_state_to_count[ref_state] = ref_count;
                codon_state_to_count[alt_state] = alt_count;

                if (Data.count(taxon) == 1 and Data.at(taxon).count(site) == 1) {
                    cerr << "There is already a SNP at site " << site << " for taxon "
                         << Alignment->GetTaxonSet()->GetTaxon(taxon) << "." << endl;
                } else {
                    if (alt_count == sample_size) {
                        Alignment->SetState(taxon, site, alt_state);
                    } else {
                        Data[taxon][site] = codon_state_to_count;
                        SiteSampleSize[taxon][site] = sample_size;
                    }
                }
            }
        }
        assert(sample_size != 0);
    }
    if (Nvcf == 0) {
        cerr << "No file (.vcf) found." << endl;
        cerr << "The files (.vcf) must contains the taxon name: e.g '"
             << '*' + Alignment->GetTaxonSet()->GetTaxon(0) + "*.vcf'" << endl;
    } else {
        cerr << "Found " << Nvcf << " files (.vcf) out of " << Alignment->GetNtaxa() << " taxa"
             << endl;
    }
};

int PolyData::GetNstate() const { return Alignment->GetNstate(); }

int PolyData::GetNtaxa() const { return Alignment->GetNtaxa(); }

int PolyData::GetNsite() const { return Alignment->GetNsite(); }

int PolyData::GetNvcf() const { return Nvcf; }

unsigned PolyData::GetCount(int taxon, int site, int state) const {
    if (Data.count(taxon) == 1 and Data.at(taxon).count(site) == 1 and
        Data.at(taxon).at(site).count(state) == 1) {
        return Data.at(taxon).at(site).at(state);
    } else {
        if (Alignment->GetState(taxon, site) == state) {
            return SampleSize.at(taxon);
        } else {
            return 0;
        }
    }
}

unsigned PolyData::GetSampleSize(int taxon, int site) const {
    if (SiteSampleSize.count(taxon) == 1 and SiteSampleSize.at(taxon).count(site) == 1) {
        return SiteSampleSize.at(taxon).at(site);
    } else {
        return SampleSize.at(taxon);
    }
}