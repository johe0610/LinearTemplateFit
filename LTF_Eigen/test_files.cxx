#include <iostream>
#include "TFile.h"
#include "TH1D.h"

// Example usage:
// TH1D* tmp_data = GetHistogram("data.root", "data_alt.root", "my_histogram");

TH1D* GetHistogram(const char* datafile, const char* datafile_alt, const char* name) {
    TFile* file1 = TFile::Open(datafile);
    if (!file1 || file1->IsZombie()) {
        std::cout << "Could not open file: " << datafile << std::endl;
        return nullptr;
    }

    TH1D* tmp_data = dynamic_cast<TH1D*>(file1->Get(name));
    if (tmp_data) {
        std::cout << "Found histogram " << name << " in file " << datafile << std::endl;
        return tmp_data; // Note: caller must manage the ownership!
    }

    std::cout << "Could not find histogram " << name << " in file " << datafile
              << ", trying alternative file " << datafile_alt << std::endl;

    TFile* file2 = TFile::Open(datafile_alt);
    if (!file2 || file2->IsZombie()) {
        std::cout << "Could not open alternative file: " << datafile_alt << std::endl;
        file1->Close();
        delete file1;
        return nullptr;
    }

    tmp_data = dynamic_cast<TH1D*>(file2->Get(name));
    if (tmp_data) {
        std::cout << "Found histogram " << name << " in alternative file " << datafile_alt << std::endl;
        file1->Close();
        delete file1;
        return tmp_data;
    }

    std::cout << "Could not find histogram " << name << " in either file." << std::endl;
    file1->Close();
    delete file1;
    file2->Close();
    delete file2;
    return nullptr;
}

// Example usage in main:
int test_files() {
    TH1D* tmp_data = GetHistogram(
				  "files/unfolding_SR_Whad_Final_l_Whad_particle_TUnfoldStandalone_OptionA_data_nonClosureAlternative_paper.root",
				  "files/unfolding_SR_Whad_Final_l_Whad_particle_TUnfoldStandalone_OptionA_data_nonClosureAlternative_extra.root",
				  "unfolding_rbj_direct_TOTAL__1up");

    if (!tmp_data) {
        std::cerr << "Histogram not found in either file. Exiting with code 1." << std::endl;
        return 1;
    }

    // Use tmp_data here...
    // Remember: TFile::Get() returns a pointer to an object owned by the file.
    // If you want to keep the histogram after closing the file, you must clone it:
    TH1D* histogram = static_cast<TH1D*>(tmp_data->Clone("cloned_histogram"));
    histogram->Print("All");
    // Now you can use histogram even after closing the file.

    // Clean up
    delete histogram;

    return 0;
}
