#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <string>
#include <map>    
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TGaxis.h"
#include "TText.h"
#include "TStyle.h"
#include "TMath.h"
#include "TTimeStamp.h"
#include "TRandom.h"
#include "TROOT.h"
#include "TF1.h"
#include "TF2.h"
#include "TString.h"
#include "TPad.h"

using namespace std;

// __________________________________________
//
// DA CAMBIARE: nome tree e file @ 197-198 
//              colore @ 296: kYellow+2
//              num label @ 275: %8
//              media sul tempo @ 209: 9000
//              titolo asse Y + u.m @ 293: He4 temperature [K]
//              range asse Y @ 295: 0.05, 0.052
//              titolo canvas @ 260: He4 evaporator bottom vs Time GMT+00 (2.5-hour interval)
//___________________________________________

// Create a custom formatter function for your timestamps
TString formatTimeStamp(double unixTime) {
    TTimeStamp ts(unixTime);
    
    // Use UInt_t instead of int
    UInt_t year, month, day, hour, min, sec;
    ts.GetDate(kFALSE, 0, &year, &month, &day);
    ts.GetTime(kFALSE, 0, &hour, &min, &sec);
    
    // Format as a shorter string
    return TString::Format("%u-%02u-%02u %02u:%02u", 
                          year, month, day, hour, min);
}

void DCS_TPO4(){

    TFile* eventFile = new TFile("histDCS-newdata.root", "READ");
    TTree* eventTree = (TTree*)eventFile->Get("OutTree");

    // Variables to read from 2023 eventTree
    int eventNum, nVertices, nBeam;
    double unixTime, timeInSpill;

    eventTree->SetBranchAddress("Nevent", &eventNum);
    eventTree->SetBranchAddress("UnixSeconds", &unixTime);
    eventTree->SetBranchAddress("TiS", &timeInSpill);
    eventTree->SetBranchAddress("NVrtx", &nVertices);
    eventTree->SetBranchAddress("Nbeam", &nBeam);

    const int timeWindow = 9000;  // 2.5 hours
    int nEvents = eventTree->GetEntries();

   // cout << "n entries in eventTree " << nEvents << endl;
    
    
    // Get the first entry to initialize min/max
    eventTree->GetEntry(0);
    double minTimeEV = unixTime;
    eventTree->GetEntry(nEvents-1);
    double maxTimeEV = unixTime;
    
    
    // Now we know our time range, calculate number of windows
    const int numWindowsEV = int((maxTimeEV - minTimeEV) / timeWindow) + 1;
    
    // Create arrays once we know the size
    int* sumVertices = new int[numWindowsEV]();  
    int* sumBeam = new int[numWindowsEV]();
    bool* windowHasRatioData = new bool[numWindowsEV]();
    double* meanVrtx = new double[numWindowsEV]();
    double* meanBeam = new double[numWindowsEV]();
    double* meanRatio = new double[numWindowsEV]();
    double* windowTime = new double[numWindowsEV]();
    vector<TString> labels(numWindowsEV); // Use TString to properly manage memory

    // counters for number of vrtx/beam in each time window
    int* VerticesPerWindow = new int[numWindowsEV](); 
    int* BeamPerWindow = new int[numWindowsEV](); 
    
    
    // Initialize window times and labels (only once)
    for (int i = 0; i < numWindowsEV; i++) {
        windowTime[i] = minTimeEV + (i * timeWindow) + (timeWindow / 2);
        labels[i] = TTimeStamp(windowTime[i]).AsString();
    }
    
    // Now collect the data in a single pass through the events
    for (int i = 0; i < nEvents; i++) {
        eventTree->GetEntry(i);
        
        int windowIndex = int((unixTime - minTimeEV) / timeWindow);
        
        if (windowIndex >= 0 && windowIndex < numWindowsEV && nBeam > 0) {
            sumVertices[windowIndex] += nVertices;
            if(nVertices > 0 )VerticesPerWindow[windowIndex] ++;
            sumBeam[windowIndex] += nBeam;
            if(nBeam>0) BeamPerWindow[windowIndex] ++;
            windowHasRatioData[windowIndex] = true;

        }
        
    }
    
    // Calculate ratios
    int validPoints = 0;
    for (int i = 0; i < numWindowsEV; i++) {
            meanVrtx[i] = (double)sumVertices[i]/VerticesPerWindow[i];
            meanBeam[i] = (double)sumBeam[i]/BeamPerWindow[i];
            meanRatio[i] = (meanVrtx[i] /meanBeam[i]);
            validPoints++;

            // cout << "Window " << i << ": sum vrtx " << sumVertices[i] << " num vrtx " << VerticesPerWindow[i] << " sum beam " << sumBeam[i] << " num beam " << BeamPerWindow[i] << endl;
            // cout << " mean ratio " << meanRatio[i] << endl;
    }
    
    // Create graph with the correct number of points
    TGraph* gRatio = new TGraph(numWindowsEV, windowTime, meanRatio);
    gRatio->SetName("gRatio");
    gRatio->SetTitle("Ratio nVertices/nBeam vs Time GMT+00 (2.5-hour interval)");
    
    TCanvas* c1 = new TCanvas("c1", "TPO4 vs Time", 1200, 800);
    c1->Divide(1, 2);
    c1->cd(2);
    
    // Set axis properties
    gRatio->GetXaxis()->SetTimeDisplay(1);
    
    // Apply labels (much more efficient approach)
    TAxis* xaxis = gRatio->GetXaxis();
    for (int i = 0; i < numWindowsEV; i++) {
        // Find the bin that contains this x value
        int bin = xaxis->FindBin(windowTime[i]);

        if (i % 8 == 0) { // Only show every 8th label
            xaxis->SetBinLabel(bin, formatTimeStamp(windowTime[i]).Data());
        }
        else {
            xaxis->SetBinLabel(bin, ""); // Empty string for other bins
        }

        //xaxis->SetBinLabel(bin, labels[i].Data());
    }
    
    gRatio->GetXaxis()->SetLabelSize(0.05);
    gRatio->GetXaxis()->SetLabelOffset(0.01);
    gRatio->GetXaxis()->LabelsOption("h");
    gRatio->GetXaxis()->SetNdivisions(510);
    
    gRatio->SetMarkerStyle(20);
    gRatio->SetMarkerSize(0.9);
    gRatio->GetYaxis()->SetRangeUser(0.9985, 1.0001);
    gRatio->GetXaxis()->SetTitle("Date");
    gRatio->GetYaxis()->SetTitle("Ratio");
    gRatio->SetMarkerColor(kRed);

   // gRatio->GetYaxis()->SetRangeUser(0.95, 1);
    gPad->SetGrid(1,1);
    gRatio->Draw("AP");
    c1->Update();


    c1->SaveAs("/eos/home-g/gmeinard/2023DCS-correlation/DCS_TPO4_vs_time.png");

    // Clean up
    delete[] sumVertices;
    delete[] sumBeam;
    delete[] windowHasRatioData;
    delete[] meanVrtx;
    delete[] meanBeam;
    delete[] meanRatio;
    delete[] windowTime;
    delete[] VerticesPerWindow;
    delete[] BeamPerWindow;
    eventFile->Close();

    // __________________________________________________________ TPO4 _____________________________________________________________________


    TFile *hfileTPO4  = new TFile("/eos/home-g/gmeinard/root_MCfiles/DCS_data/He4_evaporator_bottom_(TPO4_temperature)_2023-05-19_00_00_00_2023-06-26_00_00_00.root", "READ");
    TTree *treeTPO4 = (TTree*)hfileTPO4->Get("He4 evaporator bottom (TPO4 temperature)");
     
    // Variables to read from DCS TTrees

    float TPO4 = 0;
    TTimeStamp* time_TPO4 = nullptr;

    // trees branches
    treeTPO4->SetBranchAddress("timestamp", &time_TPO4);
    treeTPO4->SetBranchAddress("fValue",&TPO4);

    const int timeWindow_TPO4 = 9000;  // 2.5 hours
    int nTPO4 = treeTPO4->GetEntries();
    // cout << " n entries in TPO4 " << nTPO4 << endl;

    // Get the first entry to initialize min/max
    treeTPO4->GetEntry(0);
    double minTimeTPO4 = time_TPO4->AsDouble();
    // treeTPO4->GetEntry(nTPO4-1);
    // double maxTimeTPO4 = time_TPO4->AsDouble();

    // Now we know our time range, calculate number of windows
    const int numWindowsTPO4 = int((maxTimeEV - minTimeTPO4) / timeWindow_TPO4) + 1;
    
    // Create arrays once we know the size
    float* sumTPO4 = new float[numWindowsTPO4](); 
    double* windowTime_TPO4 = new double[numWindowsTPO4]();
    vector<TString> labels_TPO4(numWindowsTPO4); // Use TString to properly manage memory
    // counters for TPO4 values in each time window
    int* TPO4PerWindow = new int[numWindowsTPO4]();
    double* meanTPO4 = new double[numWindowsTPO4]();
    
    // Initialize window times and labels (only once)
    for (int i = 0; i < numWindowsTPO4; i++) {
        windowTime_TPO4[i] = minTimeTPO4 + (i * timeWindow_TPO4) + (timeWindow_TPO4 / 2);
        labels_TPO4[i] = TTimeStamp(windowTime_TPO4[i]).AsString();                    // time_TPO4->AsString();
    }
    
    // Now collect the data in a single pass through the events
    for (int i = 0; i < nTPO4; i++) {
        treeTPO4->GetEntry(i);
        
        int windowIndex_TPO4 = int((time_TPO4->AsDouble() - minTimeTPO4) / timeWindow_TPO4);
        
        if (windowIndex_TPO4 >= 0 && windowIndex_TPO4 < numWindowsTPO4) {
            sumTPO4[windowIndex_TPO4] += TPO4;
            TPO4PerWindow[windowIndex_TPO4] ++;
            // cout << " windows index " << windowIndex_TPO4 <<" tpo4 " << sumTPO4[windowIndex_TPO4] << endl;
        }
    }
    
    // Calculate mean
    int validPointsTPO4 = 0;
    for (int i = 0; i < numWindowsTPO4; i++) {
            meanTPO4[i] = (sumTPO4[i]/(double)(TPO4PerWindow[i]));
            validPointsTPO4++;
    }


    // Create graph with the correct number of points
    TGraph* gTPO4 = new TGraph(numWindowsTPO4, windowTime_TPO4, meanTPO4);
    gTPO4->SetName("gTPO4");
    gTPO4->SetTitle("He4 evaporator bottom vs Time GMT+00 (2.5-hour interval)");
    
    // TCanvas* c1 = new TCanvas("c1", "TPO4 vs Time", 1200, 800);
   // c1->Divide(1, 2);
    c1->cd(1);
    
    // Set axis properties
    gTPO4->GetXaxis()->SetTimeDisplay(1);
    
    // Apply labels (much more efficient approach)
    TAxis* xaxis_TPO4 = gTPO4->GetXaxis();
    for (int i = 0; i < numWindowsTPO4; i++) {
        // Find the bin that contains this x value
        int bin_TPO4 = xaxis_TPO4->FindBin(windowTime_TPO4[i]);

        if (i % 8 == 0) { // Only show every 8th label
            xaxis_TPO4->SetBinLabel(bin_TPO4, formatTimeStamp(windowTime_TPO4[i]).Data());
        }
        else {
            xaxis_TPO4->SetBinLabel(bin_TPO4, ""); // Empty string for other bins
        }

        //xaxis->SetBinLabel(bin, labels[i].Data());
    }
    
    gTPO4->GetXaxis()->SetLabelSize(0.05);
    gTPO4->GetXaxis()->SetLabelOffset(0.01);
    gTPO4->GetXaxis()->LabelsOption("h");
    gTPO4->GetXaxis()->SetNdivisions(510);
    
    gTPO4->SetMarkerStyle(20);
    gTPO4->SetMarkerSize(0.9);
    gTPO4->GetXaxis()->SetTitle("Date");
    gTPO4->GetYaxis()->SetTitle("He4 temperature [K]");
    gTPO4->GetYaxis()->LabelsOption("h");
    gTPO4->GetYaxis()->SetRangeUser(0.05, 0.052);
    gTPO4->SetMarkerColor(kYellow+2);

    gPad->SetGrid(1,1);
    gTPO4->Draw("AP");

    // Clean up
    delete[] sumTPO4;
    delete[] meanTPO4;
    delete[] TPO4PerWindow;
    delete[] windowTime_TPO4;
    
    c1->Update();
    c1->SaveAs("/eos/home-g/gmeinard/2023DCS-correlation/DCS_TPO4_vs_time.png");
    
    hfileTPO4->Close();
    
}
