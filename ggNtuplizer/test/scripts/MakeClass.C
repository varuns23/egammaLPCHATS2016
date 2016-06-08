void MakeClass()                                                                                                                                                               
{
    cout<<"Opening the input file."<<endl;

    //TString Inputfile = "root://cmseos.fnal.gov//store/user/leptonjets/varun/13TeV/Ntuples/Data_76X/ggNtuples_25ns_data_980.root";
    TString Inputfile = "../hgg_ggtree_mc.root";

    TFile *file = TFile::Open(Inputfile);

    cout<< " Input root file opended" << endl;
    // TTree *ttr = (TTree*)file->Get("myEvent");
    TTree *ttr = (TTree*)file->Get(Inputfile+":/ggNtuplizer/EventTree");
    //  ttr->MakeSelector("myEvent");
    ttr->MakeClass("getPhotonIsoEfficency");
}

