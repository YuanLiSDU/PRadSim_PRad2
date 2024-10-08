//
// HyCalDigitization.cc
// Developer : Chao Peng, Chao Gu
// History:
//   Aug 2012, C. Peng, Original version.
//   Jan 2017, C. Gu, Rewrite to separate the digitization.
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "HyCalDigitization.hh"

#include "PRadClusterProfile.h"
#include "PRadEventStruct.h"
#include "PRadHyCalSystem.h"
#include "PRadHyCalModule.h"
#include "ConfigParser.h"

#include "TROOT.h"
#include "TError.h"
#include "TObject.h"
#include "TChain.h"
#include "TMath.h"
#include "TRandom2.h"
#include "TF1.h"

#include <ctime>
#include <string>
#include <vector>
#include <algorithm>

Double_t asyGaussian(Double_t* x, Double_t* par)
{
    Double_t t = x[0];
    if (t < par[0]){
        return exp(-pow((t-par[0])/(sqrt(2.)*par[1]), 2));
    }
    else return exp(-pow((t-par[0])/(sqrt(2.)*par[2]), 2));
}

static TF1* asyGaussFunc = new TF1("asyGaussFunc", asyGaussian, -1000, 1000., 3);

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

static TRandom2 *RandGen = new TRandom2();

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

bool SortForID(PRadHyCalModule* a, PRadHyCalModule* b){
    return a->GetID() < b->GetID();
}


HyCalDigitization::HyCalDigitization(const std::string &abbrev, const std::string &path, double ebeam, int smearMode, int pedData) 
: StandardDigiBase(abbrev), fDMethod(0), fBeamEnergy(ebeam), fSmearMode(smearMode), fPedData(pedData)
{
    //RandGen->SetSeed((UInt_t)time(NULL));
    gRandom->SetSeed((UInt_t)time(NULL));
    RandGen->SetSeed(1);
    asyGaussFunc->SetParameters(1., 1., 1.);
    asyGaussFunc->SetNpx(500);

    fHyCal = new PRadHyCalSystem(path);

    fModuleList = fHyCal->GetModuleList();
    std::sort(fModuleList.begin(), fModuleList.end(), SortForID);

    if (fModuleList.size() != NModules)
        std::cout << "ERROR: number of modules do not match" << std::endl;

    fModuleHitList.clear();

    for (auto &module : fModuleList){
	//std::cout<<module->GetID()<<std::endl;
        //fModuleHitList.push_back(ModuleHit(module->GetID(), module->GetGeometry(), module->GetLayout(), 0, false));
        fModuleHitList.emplace_back(module, module->GetID(), 0, false);
    }
    //fProfile = &PRadClusterProfile::Instance();

    fTotalEdep = 0;

    for (int i = 0; i < NModules; i++) {
        fModuleEdep[i] = 0;
        fModuleTrackL[i] = 0;
    }
    
    LoadMCCaliConst();
    
    if (fPedData){
        if (ebeam > 1500.){
            fPedFile = new TFile("./database/norm_ped_1443.root" ,"READ");
        }else{
            fPedFile = new TFile("./database/norm_ped_1288.root" ,"READ");
        }
        
        fPedTree = (TTree*)fPedFile->Get("T");
        fPedTree->SetBranchAddress("PedEvent", fModulePed);
        fPedEntry = gRandom->Integer(fPedTree->GetEntries() - 1);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HyCalDigitization::~HyCalDigitization()
{
    //
    if (fPedData){
        fPedFile->Close();
        delete fPedFile;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void HyCalDigitization::LoadMCCaliConst()
{
    ConfigParser parser;
    std::string fileName = "./database/calibration/new_mc_cali_const.dat";
    
    if (!parser.ReadFile(fileName)){
        std::cout<<"cannot find mc calibration file, using default value 1 and sigma 0"<<std::endl;
        for (int i=0; i<T_BLOCKS; i++){
            fMCCaliConst[i] = 1.;
            fMCCaliSigma[i] = 0.;
            fResoPar[i][0] = 0.;
            fResoPar[i][1] = 0.;
            fResoPar[i][2] = 0.;
        }
        return;
    }
    int count = 0;
    while (parser.ParseLine()){
        double input[8];
        for (int i=0; i<8; i++) input[i] = parser.TakeFirst().Double();
        if (fBeamEnergy < 1500){
            fMCCaliConst[count] = fSmearMode == 0 ? input[0] : input[2];
            fMCCaliSigma[count] = fSmearMode == 0 ? input[1] : input[3];
        }
        else{
            fMCCaliConst[count] = fSmearMode == 0 ? input[4] : input[6];
            fMCCaliSigma[count] = fSmearMode == 0 ? input[5] : input[7];
        }
        count++;
        
    }
    //parser.CloseFile();
    
    if (!parser.ReadFile("./database/hycal_resolution_curve_2terms.dat")){
        std::cout<<"cannot find hycal_resolution_curve.dat"<<std::endl;
        exit(0);
    }
    
    while(parser.ParseLine()){
        std::string channelName = parser.TakeFirst();
        double input[3];
        for (int i=0; i<3; i++) input[i] = parser.TakeFirst().Double();
        
        int id = -1;
        
        if (channelName[0] == 'W' || channelName[0] == 'G'){
            id = std::stoi(channelName.substr(1, channelName.length() - 1));
            if (channelName[0] == 'W') id += 1000;
        }else{
            id = -1;
        }
        
        if (id <= 0 ) continue;
        
        for (int i=0; i<3; i++) fResoPar[id-1][i] = input[i];
        

    }
    
    //parser.CloseFile();
    
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HyCalDigitization::RegisterData(TChain *t)
{
    StandardDigiBase::RegisterData(t);

    t->SetBranchAddress(Form("%s.TotalEdep", fAbbrev), &fTotalEdep);

    if (fDMethod != 1) {
        t->SetBranchAddress(Form("%s.ModuleEdep", fAbbrev), fModuleEdep);
        t->SetBranchAddress(Form("%s.ModuleTrackL", fAbbrev), fModuleTrackL);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int HyCalDigitization::PreStart(uint32_t *buffer, int base_index)
{
    int index = 0;

    for (int roc_id = 6; roc_id >= 4; --roc_id)
        index += addRocData(&buffer[index + base_index], roc_id, index + base_index);

    return index;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

bool HyCalDigitization::ProcessEvent(uint32_t *buffer)
{
    if (fTotalEdep < TRIGGER_THRESHOLD) {
        std::cout<<fTotalEdep<<std::endl;
        Clear();
        return false;
    }

    if (fDMethod == 1) UpdateEnergy();
    
    if (fPedData){
        fPedTree->GetEntry(fPedEntry);
        fPedEntry++;
        if (fPedEntry == fPedTree->GetEntries() - 1) fPedEntry = 0;
        
        //std::cout<<"event"<<std::endl;
        //for (int i=0; i<2156; i++) std::cout<<fModulePed[i]<<std::endl;
    }

    for (int i = 0; i < NModules; i++){
        if (fModuleList[i]->IsLeadGlass()){
            FillBuffer(buffer, *(fModuleList[i]), fModuleTrackL[i]);
        }else{
            FillBuffer(buffer, *(fModuleList[i]), fModuleEdep[i]);
        }
    }
    Clear();
    return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HyCalDigitization::Clear()
{
    StandardDigiBase::Clear();

    fTotalEdep = 0;

    for (int i = 0; i < NModules; i++) {
        fModuleEdep[i] = 0;
        fModuleTrackL[i] = 0;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void HyCalDigitization::UpdateEnergy()
{
    /*auto det = fHyCal->GetDetector();

    for (int i = 0; i < fN; i++) {
        float fx = fX[i];
        float fy = fY[i];

        double fracsum = 0;
        double frac[NModules];

        int sid = det->GetSectorID(fx, fy);
        int type = det->GetSectorInfo().at(sid).mtype;

        for (int j = 0; j < NModules; j++) {
            auto &hit = fModuleHitList[j];
            double dist = det->QuantizedDist(fx, fy, sid, hit->GetX(), hit->GetY(), hit->GetSectorID());
            auto profile = fProfile->GetProfile(type, dist, fMomentum[i]);

            if (profile.frac > 0)
                frac[j] = RandGen->Gaus(profile.frac, profile.err);
            else
                frac[j] = 0;

            if (frac[j] < 0) frac[j] = 0;

            fracsum += frac[j];
        }

        for (int j = 0; j < NModules; j++)
            fModuleEdep[j] += fMomentum[i] * frac[j] / fracsum;
    }*/

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int HyCalDigitization::addRocData(uint32_t *buffer, int roc_id, int base_index)
{
    int index = 0;
    int nslot, slot[25];

    switch (roc_id) {
    default:
        return 0;

    case 4:
        nslot = 10;

        for (int i = 0; i < nslot; ++i)
            slot[i] = 22 - 2 * i;

        break;

    case 5:
    case 6:
        nslot = 10;

        for (int i = 0; i < nslot; ++i)
            slot[i] = 23 - 2 * i;

        break;
    }

    // add roc header
    buffer[index++] = 0x00000000;
    buffer[index++] = (roc_id << 16) | 0x1001;

    // add TI bank 11 words
    buffer[index++] = 0x0000000a; // 10 + 1 words in total
    buffer[index++] = 0xe10a0100; // TI bank header
    buffer[index + 2] = 2 << 24; // only 2nd word matters, it defines trigger type, here it is total sum
    index += 9; // TI bank expects 9 words

    buffer[index++] = 0x00000000;
    buffer[index++] = 0xe1200100; // Fastbus bank header
    // roc id and board number
    buffer[index++] = 0xdc0adc00 | ((roc_id & 0xf) << 20) | (nslot & 0xff);

    for (int i = 0; i < nslot; ++i) {
        buffer[index++] = (slot[i] << 27) | 65;
        data_index[(6 - roc_id) * 10 + i] = index + base_index;

        for (int ch = 0; ch < 64; ++ch)
            buffer[index++] = (slot[i] << 27) | (ch << 17);
    }

    buffer[index++] = 0xfabc0005; // end word
    buffer[0] = index - 1; // roc bank size
    buffer[13] = index - 14; // fastbus bank size

    return index;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HyCalDigitization::FillBuffer(uint32_t *buffer, const PRadHyCalModule &module, double edep)
{
    int crate = module.GetChannel()->GetAddress().crate;
    int slot = module.GetChannel()->GetAddress().slot;
    int channel = module.GetChannel()->GetAddress().channel;

    int pos = (6 - crate) * 10 + ((23 - slot) / 2);

    int index = data_index[pos] + channel;

    double ped = RandGen->Gaus(module.GetChannel()->GetPedestal().mean, module.GetChannel()->GetPedestal().sigma);
    
    if (fPedData) ped = fModulePed[module.GetID()-1];
    
    //std::cout<<module.GetID()<<" "<<ped<<std::endl;
    
    unsigned short val = 0;
    
    double mcConst = fMCCaliConst[module.GetID()-1];
    double mcSigma = fMCCaliSigma[module.GetID()-1];

    //ped = module.GetChannel()->GetPedestal().mean; //test
    
    if (!module.GetChannel()->IsDead() && edep > 0 && edep < 10000) {
        /*if (module.IsLeadTungstate()) {
            double reso = (0.026*TMath::Sqrt(0.73)/TMath::Sqrt(edep/1000.)+mcSigma);
            if (reso < 0.) reso = 0.;
            if (module.GetID() == 1491) std::cout<<edep<<" "<<reso<<std::endl;
            val = ped + (RandGen->Gaus(edep,  edep * reso )) *  (mcConst / module.GetCalibrationFactor());
        } else if (module.IsLeadGlass()) {
            double reso = (0.053*TMath::Sqrt(0.73)/TMath::Sqrt(edep/1000.)+mcSigma);
            if (reso < 0.) reso = 0.;
            
            val = ped + (RandGen->Gaus(edep,  edep * reso )) *  (mcConst / module.GetCalibrationFactor());
        }
        */
        double reso = TMath::Sqrt(0.76)*TMath::Sqrt( pow(fResoPar[module.GetID()-1][0]/TMath::Sqrt(edep/1000.), 2) + 
                                                     pow(fResoPar[module.GetID()-1][1]/(edep/1000.), 2)            + 
                                                     pow(fResoPar[module.GetID()-1][2], 2) );
                
        reso += mcSigma;

        if (reso < 0.) reso = 0.;
        
        //reso = 0.; //test

        //val = ped + (RandGen->Gaus(edep,  edep * reso )) *  (mcConst / module.GetCalibrationFactor());

        //reso = 0.;
        
        double tmpval = 0.;
        if (edep > 10. && reso > 1e-3 && (module.GetID() < 0)) {
            asyGaussFunc->SetParameters(0., edep*reso*0.5, edep*reso*1.45);
            tmpval = asyGaussFunc->GetRandom() + edep;
        }
        else tmpval = (RandGen->Gaus(edep,  edep * reso ));

        val = ped + tmpval * (mcConst / module.GetCalibrationFactor());
    } else
        val = ped;

    buffer[index] = (slot << 27) | (channel << 17) | val;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
