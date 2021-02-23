// Generator 2pi2pi0
// Everything should be in GeV
#include "Cmd3Generator2pi2pi0_ke.h"
#include "TClassHolder.h"
#include <TGenPhaseSpace.h>
#include <iostream>
#include <gsl/gsl_rng.h>
#include "globals.hh"
#include <cmath>
#include <stdlib.h>
#include <stdio.h>
#include <complex>
#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include "Randomize.hh"
#include <TRandom1.h>
#include <TRandom.h>
#include <TRandom3.h>
#include <TVector3.h>
#include "TLorentzVectorC.h"
#include "TMath.h"
#include "TF1.h"
#include "TF2.h"
#include "TGraph2D.h"
#include <math.h>
#include <fstream>



typedef HepMC::FourVector HepLorentzVector;

//--------------------------------------------------------------------------------------------------------------------------------------------------
static TClassHolder<Cmd3Generator2pi2pi0_ke, Cmd3PrimaryGenerator> fpipi2pi0_ke_gen("pipi2pi0_ke_gen");
typedef HepMC::FourVector HepLorentzVector;
//**************************************************************************************************************************************************

Cmd3Generator2pi2pi0_ke::Cmd3Generator2pi2pi0_ke(){}


double Cmd3Generator2pi2pi0_ke::factor_ampl(double enrgy, int number){
ifstream stream1(results.c_str());
TString btv;
double energy;
double energy0;
double data[100];
double data0[100];
double res = -10000;
bool check = false;
int npar = number;
while(stream1.eof()==0){
        int i = -1;
        stream1 >> energy;
        if(stream1.eof()==1)break;
        while(btv!="end"){
            if(i>-1){
                data[i] = btv.Atof();
            }
            stream1 >> btv;
            if(stream1.eof()==1)break;
            i++;
        }
        btv = "as";
        if(enrgy <= energy){
            res = data[npar];
            if(check == true){
                res = data[npar]*(enrgy - energy0)/(energy - energy0) + data0[npar]*(- enrgy + energy)/(energy - energy0);
                if(number%4==2)res = data0[npar];
            }
            break;
        }
        check = true;
        for(int i = 0; i < 100; i++){data0[i] = data[i];}
        energy0 = energy;
  }
  if(res < -1000)res = data[npar];
  return res;
}


double Cmd3Generator2pi2pi0_ke::factor_ampl1(double enrgy, int number){

double res = 0;
double width = 0.09;
double step = 0.00005;
double en = enrgy - width/2.;
while(en < (enrgy + width/2.)){

        res = res+step*factor_ampl(en, number);
        en = en + step;
}

return res/width;
}

double Cmd3Generator2pi2pi0_ke::cross_norm(double *e, double *par){

    ifstream streamd(("/home/eakozyrev/diskD/4pi_new/mygenerator/data/new/"+scross_norm).c_str());
    int i = 0;double norm[1000];
    while(streamd.eof()==0){streamd >> norm[i]; i++; }
    double factor;
    if(e[0] <= 0.625)factor = norm[0];
    else if(e[0] >= 1.975)factor = norm[27];
    else {
        int index = (int)((e[0] - 0.625)/0.05);
        double e1 = 0.625 + index*0.05;
        double e2 = 0.625 + (index+1.)*0.05;
        factor = norm[index]*(e2-e[0])/(e2-e1) + norm[index+1]*(-e1+e[0])/(e2-e1);
    }
    return factor;
}

void Cmd3Generator2pi2pi0_ke::change(int num, double value){

        if(num == 1)omegapi0_ = value;
        if(num == 2)a1pi_rhopi_ = value;
        if(num == 3)a1pi_rhopi_phase = value;
        if(num == 4)a1pi_sigmapi_ = value;
        if(num == 5)a1pi_sigmapi_phase = value;
        if(num == 6)rhof0_ = value;
        if(num == 7)rhof0phase_ = value;
        if(num == 8)rhoprhom_ = value;
        if(num == 9)rhoprhom_phase = value;
        if(num == 10)rhosigma_ = value;
        if(num == 11)rhosigmaphase_ = value;
        if(num == 12)phsp_ = value;
        if(num == 13)phspphase_ = value;
        if(num == 14)a2pi_ = value;
        if(num == 15)a2pi_ph = value;
        if(num == 16)h1pi = value;
        if(num == 17)h1pi_ph = value;
        if(num == 18)a1pi_rhopi_II = value;
        if(num == 19)a1pi_rhopi_phase_II = value;
}


void Cmd3Generator2pi2pi0_ke::fill_parameters(){

    omegapi0_ = factor_ompi0();
    omegapi0phase_ = 0.;
    
    a1pi_rhopi_ = factor_a1pi_rhopi();
    a1pi_rhopi_phase = factor_a1pi_rhopi_ph();

    a1pi_sigmapi_ = factor_a1pi_sigmapi();
    a1pi_sigmapi_phase = factor_a1pi_sigmapi_ph();

    rhof0_ = factor_rho_f0();
    rhof0phase_ = factor_rho_f0_ph();

    rhosigma_ = factor_rho_sigma();
    rhosigmaphase_ = factor_rho_sigma_ph();

    phsp_ = factor_phsp();
    phspphase_ = factor_phsp_ph();

    rhoprhom_ = factor_rhop_rhom();
    rhoprhom_phase = factor_rhop_rhom_ph();

    h1pi = factor_h1_rhopi();
    h1pi_ph = factor_h1_rhopi_ph();
    
    a2pi_ = factor_a2pi();
    a2pi_ph = factor_a2pi_ph();
    a1pi_rhopi_II =  factor_ampl(Energy,36);
    a1pi_rhopi_phase_II = factor_ampl(Energy,38);

    rhof0phase_inter = factor_ampl(Energy,42);
    rhof0_inter =  factor_ampl(Energy,40);

    rhof2 = factor_ampl(Energy,44);
    rhof2_phase = factor_ampl(Energy,46);
    rhoprhom_II = factor_ampl(Energy,48);
    rhoprhom_phase_II = factor_ampl(Energy,50);
    
    double total = omegapi0_+a1pi_rhopi_+a1pi_sigmapi_+rhof0_+rhosigma_+phsp_;
    total = total + rhoprhom_ + h1pi + a2pi_ + a1pi_rhopi_II;
    
    fompi0 = omegapi0_/total;dompi0 = factor_ampl(Energy,3);
    fa1pi_rhopi = a1pi_rhopi_/total; da1pi_rhopi = factor_ampl(Energy,5);
    fa1pi_sigmapi = a1pi_sigmapi_/total; da1pi_sigmapi = factor_ampl(Energy,9);
    frhoprhom = rhoprhom_/total; drhoprhom = factor_ampl(Energy,17);
    frhof0 = rhof0_/total; drhof0 = factor_ampl(Energy,13);
    frhosigma = rhosigma_/total; drhosigma = factor_ampl(Energy,21);
    fh1pi=h1pi/total; dh1pi= factor_ampl(Energy,33);
    
    
    
}


void Cmd3Generator2pi2pi0_ke::print(){
        cout << "============================================================" << endl;
        cout << "Energy = " << Energy << endl;
        cout << "isISR = " << isISR << endl; 

        cout << "omegapi0_ = " << omegapi0_ << endl;

        cout << "a1pi_rhopi_ = " << a1pi_rhopi_ << endl;
        cout << "a1pi_rhopi_phase = " << a1pi_rhopi_phase << endl;

        cout << "a1pi_sigmapi_ = " << a1pi_sigmapi_ << endl;
        cout << "a1pi_sigmapi_phase = " << a1pi_sigmapi_phase << endl;

        cout << "rhoprhom_ = " << rhoprhom_ << endl;
        cout << "rhoprhom_phase = " << rhoprhom_phase << endl;
        cout << "rhoprhom_II  = " << rhoprhom_II << endl;
        cout << "rhoprhom_phase_II = " << rhoprhom_phase_II << endl;

        cout << "rhof0_ = " << rhof0_ << endl;
        cout << "rhof0phase_ = " <<rhof0phase_ << endl;

        cout << "rhof0_inter = " << rhof0_inter << endl;
        cout << "rhof0phase_inter = " <<rhof0phase_inter << endl;

        cout << "rhosigma_ = " << rhosigma_ << endl;
        cout << "rhosigmaphase_ = " << rhosigmaphase_ << endl;

        cout << "rhof2 = " << rhof2 << endl;
        cout << "rhof2_phase = " << rhof2_phase << endl;

        cout << "h1pi = " << h1pi << endl;
        cout << "h1pi_ph = " << h1pi_ph << endl;

        cout << "a1pi_rhopi_II = " << a1pi_rhopi_II << endl;
        cout << "a1pi_rhopi_phase_II = " << a1pi_rhopi_phase_II << endl;

        cout << "phsp_ = " << phsp_ << endl;
        cout << "phspphase_ = " << phspphase_ << endl;

        cout << "a2pi_ = " << a2pi_ << endl;
        cout << "a2pi_ph = " << a2pi_ph << endl;
        cout << "============================================================" << endl;
}

double Cmd3Generator2pi2pi0_ke::levicivita(int A, int B, int C, int D){
        if(A == 1 && B == 2 && C == 3 && D == 4)return 1;
        else if(A == 1 && B == 2 && C == 4 && D == 3)return -1;
        else if(A == 1 && B == 3 && C == 2 && D == 4)return 1;
        else if(A == 1 && B == 3 && C == 4 && D == 2)return -1;
        else if(A == 1 && B == 4 && C == 2 && D == 3)return 1;
        else if(A == 1 && B == 4 && C == 3 && D == 2)return -1;
        else if(A == 2 && B == 3 && C == 4 && D == 1)return 1;
        else if(A == 2 && B == 4 && C == 3 && D == 1)return -1;
        else if(A == 3 && B == 2 && C == 4 && D == 1)return 1;
        else if(A == 3 && B == 4 && C == 2 && D == 1)return -1;
        else if(A == 4 && B == 2 && C == 3 && D == 1)return 1;
        else if(A == 4 && B == 3 && C == 2 && D == 1)return -1;
        else if(A == 3 && B == 4 && C == 1 && D == 2)return 1;
        else if(A == 4 && B == 3 && C == 1 && D == 2)return -1;
        else if(A == 2 && B == 4 && C == 1 && D == 3)return 1;
        else if(A == 4 && B == 2 && C == 1 && D == 3)return -1;
        else if(A == 2 && B == 3 && C == 1 && D == 4)return 1;
        else if(A == 3 && B == 2 && C == 1 && D == 4)return -1;
        else if(A == 4 && B == 1 && C == 2 && D == 3)return 1;
        else if(A == 3 && B == 1 && C == 2 && D == 4)return -1;
        else if(A == 4 && B == 1 && C == 3 && D == 2)return 1;
        else if(A == 2 && B == 1 && C == 3 && D == 4)return -1;
        else if(A == 3 && B == 1 && C == 4 && D == 2)return 1;
        else if(A == 2 && B == 1 && C == 4 && D == 3)return -1;
        else return 0.;
}

double Cmd3Generator2pi2pi0_ke::TLorentzVectorC_index(TLorentzVectorC P, int index){
        
        if(index == 0)return P.E().real();
        else if(index == 1)return P.X().real();
        else if(index == 2)return P.Y().real();
        else if(index == 3)return P.Z().real();
        else return 0.;
        
}

double Cmd3Generator2pi2pi0_ke::svertka_leviciv(TLorentzVectorC P1, TLorentzVectorC P2, TLorentzVectorC P3, TLorentzVectorC P4){
        
        double res = 0.;
        for(int i = 0; i < 4; i++){
                for(int j = 0; j < 4; j++){
                        for(int k = 0; k < 4; k++){
                                for(int l = 0; l < 4; l++){
                                
                                        res = res + levicivita(i,j,k,l)*TLorentzVectorC_index(P1,i)*TLorentzVectorC_index(P2,j)*TLorentzVectorC_index(P3,k)*TLorentzVectorC_index(P4,l);
                                        
                                }
                        }
                }
        }
        return res;
}

TLorentzVectorC Cmd3Generator2pi2pi0_ke::sver_vec_leviciv(TLorentzVectorC P2, TLorentzVectorC P3, TLorentzVectorC P4){
        
        double res0 = 0.;
        double res1 = 0.;
        double res2 = 0.;
        double res3 = 0.;
        for(int j = 0; j < 4; j++){
                for(int k = 0; k < 4; k++){
                        for(int l = 0; l < 4; l++){
                                
                                res0 = res0 + levicivita(1,j,k,l)*TLorentzVectorC_index(P2,j)*TLorentzVectorC_index(P3,k)*TLorentzVectorC_index(P4,l);
                                res1 = res1 + levicivita(2,j,k,l)*TLorentzVectorC_index(P2,j)*TLorentzVectorC_index(P3,k)*TLorentzVectorC_index(P4,l);
                                res2 = res2 + levicivita(3,j,k,l)*TLorentzVectorC_index(P2,j)*TLorentzVectorC_index(P3,k)*TLorentzVectorC_index(P4,l);
                                res3 = res3 + levicivita(4,j,k,l)*TLorentzVectorC_index(P2,j)*TLorentzVectorC_index(P3,k)*TLorentzVectorC_index(P4,l);

                        }
                }
        }
        
        return TLorentzVectorC(res0, res1, res2, res3);
}



// ------------------------------- PHASE SPACE------------------------------------

double Cmd3Generator2pi2pi0_ke::cross_ph_sp(double *e, double *par){
    double norm[] = {0.027171, 0.0258879, 0.0247892, 0.0238614, 0.0230847, 0.0223951, 0.021841, 0.0211563, 0.0207638, 0.0203469, 0.0200074, 0.0196746, 0.0192915, 0.0190267, 0.0187609, 0.0185032, 0.0182173, 0.0180387, 0.0177281, 0.0174934, 0.017325, 0.0172337, 0.0170029, 0.0168753, 0.0167488, 0.0165608, 0.0164036, 0.0162988};
    double factor;
    if(e[0] <= 0.625)factor = norm[0];
    else if(e[0] >= 1.975)factor = norm[27];
    else {
        int index = (int)((e[0] - 0.625)/0.05);
        double e1 = 0.625 + index*0.05;
        double e2 = 0.625 + (index+1.)*0.05;
        factor = norm[index]*(e2-e[0])/(e2-e1) + norm[index+1]*(-e1+e[0])/(e2-e1);
    }
    return factor;
}

TLorentzVectorC  Cmd3Generator2pi2pi0_ke::H_ph_sp_beforepi0sim(TLorentzVectorC P, TLorentzVectorC P1, TLorentzVectorC P2, TLorentzVectorC P3, TLorentzVectorC P4){

    TLorentzVectorC result = (P3-P4);//(0,1./sqrt(2),1./sqrt(2),0);
    //TLorentzVectorC result = result0*(P1+P2).Dot(P3);
    //result = result - result0*(P1+P2).Dot(P4);
    double ener = P.E().real();
    return result/sqrt(cross_ph_sp(&ener, &ener));


}



TLorentzVectorC  Cmd3Generator2pi2pi0_ke::H_ph_sp(TLorentzVectorC P, TLorentzVectorC P1, TLorentzVectorC P2, TLorentzVectorC P3, TLorentzVectorC P4){

    return H_ph_sp_beforepi0sim(P,P1,P2,P3,P4) + H_ph_sp_beforepi0sim(P,P2,P1,P3,P4) - H_ph_sp_beforepi0sim(P,P1,P2,P4,P3) - H_ph_sp_beforepi0sim(P,P2,P1,P4,P3);

}


double Cmd3Generator2pi2pi0_ke::factor_phsp(){
    return factor_ampl(Energy,24);
}

double Cmd3Generator2pi2pi0_ke::factor_phsp_ph(){
    return factor_ampl(Energy,26);
}


// ------------------------- OMEGA PI0 -----------------------------------------------------

double Cmd3Generator2pi2pi0_ke::cross_omega_pi0(double *e, double *par){
    double norm[] = {1.87049e-07, 1.56049e-06, 7.99364e-06, 3.37705e-05, 0.000135718, 0.000596526, 0.0107153, 0.100773, 0.216074, 0.309316, 0.401808, 0.4623,0.538993, 0.560054, 0.602989, 0.667697, 0.667369, 0.785723, 0.806614, 0.835757, 0.854405, 0.862221, 0.87398, 0.800103, 0.758419, 0.752406, 0.783933, 0.925019};
    double factor;
    if(e[0] <= 0.625)factor = norm[0];
    else if(e[0] >= 1.975)factor = norm[27];
    else {
        int index = (int)((e[0] - 0.625)/0.05);
        double e1 = 0.625 + index*0.05;
        double e2 = 0.625 + (index+1.)*0.05;
        factor = norm[index]*(e2-e[0])/(e2-e1) + norm[index+1]*(-e1+e[0])/(e2-e1);
    }
    return factor;
}

TLorentzVectorC  Cmd3Generator2pi2pi0_ke::H_omega_pi0_init(TLorentzVectorC P, TLorentzVectorC P1, TLorentzVectorC P2, TLorentzVectorC P3, TLorentzVectorC P4){


    TLorentzVectorC W = P2 + P3 + P4;
    TLorentzVectorC k = P3 + P4;
    TLorentzVectorC result = (P3-P4)*P.Dot(k)*W.Dot(W) + W*P.Dot(P3-P4)*W.Dot(k) + k*P.Dot(W)*W.Dot(P3-P4);
    result = result - W*P.Dot(k)*W.Dot(P3-P4) - k*P.Dot(P3-P4)*W.Dot(W) - (P3-P4)*P.Dot(W)*W.Dot(k);
    result = result*propagator_omega(W.M()*W.M())*propagator_rho(k.M()*k.M());
    return result;
}

TLorentzVectorC Cmd3Generator2pi2pi0_ke::H_omega_pi0_beforepi0sim(TLorentzVectorC P, TLorentzVectorC P1, TLorentzVectorC P2, TLorentzVectorC P3, TLorentzVectorC P4){
    return - H_omega_pi0_init(P,P1,P2,P3,P4) + H_omega_pi0_init(P,P1,P4,P3,P2) + H_omega_pi0_init(P,P1,P3,P2,P4);
}

TLorentzVectorC  Cmd3Generator2pi2pi0_ke::H_omega_pi0(TLorentzVectorC P, TLorentzVectorC P1, TLorentzVectorC P2, TLorentzVectorC P3, TLorentzVectorC P4){
    TLorentzVectorC a = H_omega_pi0_beforepi0sim(P,P1,P2,P3,P4) + H_omega_pi0_beforepi0sim(P,P2,P1,P3,P4);
    double ener = P.E().real();
    return  a/sqrt(cross_omega_pi0(&ener, &ener));
}


double Cmd3Generator2pi2pi0_ke::factor_ompi0(){
    return factor_ampl(Energy,2);
}
double Cmd3Generator2pi2pi0_ke::factor_ompi0_ph(){
    return -1.;
    return factor_ampl(Energy,2);
}

//============================================================
// ------------A1 (rho pi) PI first part----------------------
//============================================================
double Cmd3Generator2pi2pi0_ke::cross_a1_rhopi_pi(double *e, double *par){
    double norm[] = {9.38941e-06, 3.37591e-05, 0.000104477, 0.000286863, 0.000767181, 0.00200169, 0.00523157, 0.0145385, 0.0411876, 0.114388, 0.2936, 0.639556, 1.20124, 2.06013, 3.09992, 4.44872, 6.03592, 7.64013, 9.14649, 10.8322, 11.9934, 12.988, 14.3327, 14.3168, 15.4297, 16.1528, 16.5091, 17.1639};
    double factor;
    if(e[0] <= 0.625)factor = norm[0];
    else if(e[0] >= 1.975)factor = norm[27];
    else {
        int index = (int)((e[0] - 0.625)/0.05);
        double e1 = 0.625 + index*0.05;
        double e2 = 0.625 + (index+1.)*0.05;
        factor = norm[index]*(e2-e[0])/(e2-e1) + norm[index+1]*(-e1+e[0])/(e2-e1);
    }
    return factor;
}

TLorentzVectorC  Cmd3Generator2pi2pi0_ke::H_a1_rhopi_pi_beforepi0sim(TLorentzVectorC P, TLorentzVectorC P1, TLorentzVectorC P2, TLorentzVectorC P3, TLorentzVectorC P4){
    TLorentzVectorC W = P1 + P2 + P3;
    TLorentzVectorC k = P2 + P3;
    TLorentzVectorC result(0.,0.,0.,0.);
    result = result + (P3-P2)*P.Dot(W)*W.Dot(k);
    result = result - W*P.Dot(P3-P2)*W.Dot(k);
    result = result + W*P.Dot(k)*W.Dot(P3-P2);
    result = result - k*W.Dot(P3-P2)*P.Dot(W);
    double ener = P.E().real();
    double ph = 0.;//1.*(W.M()*W.M()/Ma1/Ma1-1.);
    complex<double>  FF = complex<double> (cos(ph),sin(ph));
    return -result*propagator_rho(k.M()*k.M())*propagator_a1(W.M()*W.M())/sqrt(cross_a1_rhopi_pi(&ener, &ener))*FF; 
}

TLorentzVectorC  Cmd3Generator2pi2pi0_ke::H_a1_rhopi_pi(TLorentzVectorC P, TLorentzVectorC P1, TLorentzVectorC P2, TLorentzVectorC P3, TLorentzVectorC P4){
        
    TLorentzVectorC result(0.,0.,0.,0.);
    result = result + H_a1_rhopi_pi_beforepi0sim(P,P1,P2,P3,P4);
    result = result + H_a1_rhopi_pi_beforepi0sim(P,P2,P1,P3,P4);
    result = result - H_a1_rhopi_pi_beforepi0sim(P,P1,P2,P4,P3);
    result = result - H_a1_rhopi_pi_beforepi0sim(P,P2,P1,P4,P3);

    return result;

}


double Cmd3Generator2pi2pi0_ke::factor_a1pi_rhopi(){
    return factor_ampl(Energy,4);
}
double Cmd3Generator2pi2pi0_ke::factor_a1pi_rhopi_ph(){
    return factor_ampl(Energy,6);
}

//=============================================
// ------A1 (rho pi) PI second part -----------
//=============================================

double Cmd3Generator2pi2pi0_ke::cross_a1_rhopi_pi_II(double *e, double *par){
    double norm[] = {1.89775e-06, 6.93526e-06, 2.02051e-05, 5.56619e-05, 0.000144712, 0.000371023, 0.000920294, 0.00242596, 0.00627131, 0.0179908, 0.0418408, 0.0831727, 0.159426, 0.26784, 0.408592, 0.63246, 0.920412, 1.38438, 2.03832, 3.08834, 4.68597, 7.24819, 12.4944, 19.9858, 32.6581, 45.0234, 56.3054, 64.2716};
    double factor;
    if(e[0] <= 0.625)factor = norm[0];
    else if(e[0] >= 1.975)factor = norm[27];
    else {
        int index = (int)((e[0] - 0.625)/0.05);
        double e1 = 0.625 + index*0.05;
        double e2 = 0.625 + (index+1.)*0.05;
        factor = norm[index]*(e2-e[0])/(e2-e1) + norm[index+1]*(-e1+e[0])/(e2-e1);
    }
    return factor;
}

TLorentzVectorC  Cmd3Generator2pi2pi0_ke::H_a1_rhopi_pi_beforepi0sim_II(TLorentzVectorC P, TLorentzVectorC P1, TLorentzVectorC P2, TLorentzVectorC P3, TLorentzVectorC P4){

    double ener = P.E().real();

    TLorentzVectorC W = P1 + P2;
    TLorentzVectorC e = P3 - P4;
    W = P1 + P2 + P4;
    TLorentzVectorC k = P1 + P4;
    e = P4 - P1;
    TLorentzVectorC result = (W*P.Dot(P3) - P3*W.Dot(P))*(e.Dot(P2)*k.Dot(W)-e.Dot(W)*k.Dot(P2));
    result = result*propagator_pi1300(W.Dot(W))*propagator_rho(k.Dot(k));
 
    return result/sqrt(cross_a1_rhopi_pi_II(&ener, &ener));

}


TLorentzVectorC  Cmd3Generator2pi2pi0_ke::H_a1_rhopi_pi_II(TLorentzVectorC P, TLorentzVectorC P1, TLorentzVectorC P2, TLorentzVectorC P3, TLorentzVectorC P4){
        
    TLorentzVectorC result(0.,0.,0.,0.);
    result = result + H_a1_rhopi_pi_beforepi0sim_II(P,P1,P2,P3,P4);
    //result = result + H_a1_rhopi_pi_beforepi0sim_II(P,P2,P1,P3,P4);
    //result = result - H_a1_rhopi_pi_beforepi0sim_II(P,P1,P2,P4,P3);
    //result = result - H_a1_rhopi_pi_beforepi0sim_II(P,P2,P1,P4,P3);
    return result;
}

double Cmd3Generator2pi2pi0_ke::factor_a1pi_rhopi_II(){
    return factor_ampl(Energy,36);
}
double Cmd3Generator2pi2pi0_ke::factor_a1pi_rhopi_ph_II(){
    return factor_ampl(Energy,38);
}


//==================================================
// ---------A1 (sigma pi) PI -----------------------
//==================================================

double Cmd3Generator2pi2pi0_ke::cross_a1_sigmapi_pi(double *e, double *par){
    double norm[] = {2.35448e-05, 7.50863e-05, 0.000206688, 0.000509195, 0.00120981, 0.00277813, 0.00631554, 0.0142642, 0.0309426, 0.0652028, 0.135337, 0.268977, 0.499254, 0.88869, 1.45045, 2.3462, 3.49934, 4.77602, 6.22402, 7.7155, 9.27811, 10.428, 11.9285, 13.1289, 14.2201, 15.8061, 16.4864, 18.3851};
    double factor;
    if(e[0] <= 0.625)factor = norm[0];
    else if(e[0] >= 1.975)factor = norm[27];
    else {
        int index = (int)((e[0] - 0.625)/0.05);
        double e1 = 0.625 + index*0.05;
        double e2 = 0.625 + (index+1.)*0.05;
        factor = norm[index]*(e2-e[0])/(e2-e1) + norm[index+1]*(-e1+e[0])/(e2-e1);
    }
    return factor;
}

TLorentzVectorC  Cmd3Generator2pi2pi0_ke::H_a1_sigmapi_pi_beforepi0sim(TLorentzVectorC P, TLorentzVectorC P1, TLorentzVectorC P2, TLorentzVectorC P3, TLorentzVectorC P4){
        TLorentzVectorC W = P1 + P2 + P3;
        TLorentzVectorC k = P1 + P2;
        double ener = P.E().real();
        TLorentzVectorC result(0.,0.,0.,0.);
        result = result + P3*P.Dot(W)*k.Dot(W);
        result = result - k*P.Dot(W)*P3.Dot(W);
        result = result - W*P3.Dot(P)*k.Dot(W);
        result = result + W*P3.Dot(W)*k.Dot(P);
        double ph = 0.;//3.*(W.M()*W.M()/Ma1/Ma1-1.);
        complex<double>  FF = complex<double> (cos(ph),sin(ph));
        result = result*propagator_a1(W.M()*W.M())*propagator_sigma(k.M()*k.M())*FF;
        return -result/sqrt(cross_a1_sigmapi_pi(&ener, &ener));
}

TLorentzVectorC  Cmd3Generator2pi2pi0_ke::H_a1_sigmapi_pi(TLorentzVectorC P, TLorentzVectorC P1, TLorentzVectorC P2, TLorentzVectorC P3, TLorentzVectorC P4){
    TLorentzVectorC res(0.,0.,0.,0.);
    res = res + H_a1_sigmapi_pi_beforepi0sim(P,P1,P2,P3,P4) - H_a1_sigmapi_pi_beforepi0sim(P,P1,P2,P4,P3);
    res = res + H_a1_sigmapi_pi_beforepi0sim(P,P2,P1,P3,P4) - H_a1_sigmapi_pi_beforepi0sim(P,P2,P1,P4,P3);
    return res;
}

double Cmd3Generator2pi2pi0_ke::factor_a1pi_sigmapi(){
    return factor_ampl(Energy,8);
}
double Cmd3Generator2pi2pi0_ke::factor_a1pi_sigmapi_ph(){
    return factor_ampl(Energy,10);
}


//===============================================
//------------------rho f2-----------------------
//===============================================

double Cmd3Generator2pi2pi0_ke::cross_rho_f2(double *e, double *par){
    double norm[] = {0.00203014, 0.00305518, 0.00340913, 0.00319095, 0.00271788, 0.00212643, 0.00159291, 0.00111007, 0.000769828, 0.000456646, 0.000340798, 0.000284551, 0.000250472, 0.000240615, 0.00024844, 0.000257863, 0.000258196, 0.000253757, 0.000262872, 0.000250844, 0.000238593, 0.000216461, 0.000181339, 0.000164848, 0.000141188, 0.000147802, 0.000161769, 0.000176241};
    double factor;
    if(e[0] <= 0.625)factor = norm[0];
    else if(e[0] >= 1.975)factor = norm[27];
    else {
        int index = (int)((e[0] - 0.625)/0.05);
        double e1 = 0.625 + index*0.05;
        double e2 = 0.625 + (index+1.)*0.05;
        factor = norm[index]*(e2-e[0])/(e2-e1) + norm[index+1]*(-e1+e[0])/(e2-e1);
    }
    return factor;
}


TLorentzVectorC  Cmd3Generator2pi2pi0_ke::H_rho_f2(TLorentzVectorC P, TLorentzVectorC P1, TLorentzVectorC P2, TLorentzVectorC P3, TLorentzVectorC P4){

    double ener = P.E().real();
    TLorentzVectorC W = P1 + P2;
    TLorentzVectorC e = P3 - P4;

    double nomin = W.Dot(W).real();
    TLorentzVectorC result = (P2 - W*W.Dot(P2)/nomin)*(e.Dot(P1) - e.Dot(W)*P1.Dot(W)/nomin)/2.;
    result = result + (P1 - W*W.Dot(P1)/nomin)*(e.Dot(P2) - e.Dot(W)*P2.Dot(W)/nomin)/2.;
    result = result - (e-W*W.Dot(e)/nomin)*(P1.Dot(P2) - W.Dot(P1)*W.Dot(P2)/nomin)/3.;
    result = result - P*result.Dot(P)/P.Dot(P).real();
    result = result*propagator_rho((P3+P4).Dot(P3+P4))*propagator_f2(W.Dot(W));   

    return result/sqrt(cross_rho_f2(&ener, &ener));

}


// =======================================================
//----------------------- A2 PI---------------------------
// =======================================================
double Cmd3Generator2pi2pi0_ke::cross_a2_pi(double *e, double *par){
    double norm[] = {5.09254e-14, 1.12193e-12, 1.03475e-11, 6.11362e-11, 2.81282e-10, 1.07244e-09, 3.61526e-09, 1.11212e-08, 3.23475e-08, 8.83757e-08, 2.39315e-07, 6.36937e-07, 1.67296e-06, 4.37807e-06, 1.14225e-05, 2.95183e-05, 7.61204e-05, 0.000201503, 0.00052664, 0.00121244, 0.00256362, 0.0049281, 0.00854306, 0.0138288, 0.0221141, 0.0323682, 0.0459657, 0.0659924};
    double factor;
    if(e[0] <= 0.625)factor = norm[0];
    else if(e[0] >= 1.975)factor = norm[27];
    else {
        int index = (int)((e[0] - 0.625)/0.05);
        double e1 = 0.625 + index*0.05;
        double e2 = 0.625 + (index+1.)*0.05;
        factor = norm[index]*(e2-e[0])/(e2-e1) + norm[index+1]*(-e1+e[0])/(e2-e1);
    }
    return factor;
}

TLorentzVectorC  Cmd3Generator2pi2pi0_ke::H_a2_pi_beforepi0sim(TLorentzVectorC P, TLorentzVectorC P1, TLorentzVectorC P2, TLorentzVectorC P3, TLorentzVectorC P4){

    TLorentzVectorC W = P1 + P2 + P3;
    double ener = P.E().real();
    TLorentzVectorC result = - sver_vec_leviciv(W,P4,P1)*svertka_leviciv(W, P3-P2, P1, P4);
    TLorentzVectorC btv = (P3-P2)*W.Dot(W)*P4.Dot(P1) + P1*P4.Dot(W)*(P3-P2).Dot(W);
    btv = btv + W*W.Dot(P1)*P4.Dot(P3-P3) - P1*W.Dot(W)*P4.Dot(P3-P3);
    btv = btv - W*P4.Dot(P1)*(P3-P2).Dot(W) - (P3-P2)*W.Dot(P1)*W.Dot(P4);
    result = result + btv*(P4.Dot(P1) - W.Dot(P4)*W.Dot(P1)/W.Dot(W));
    result = result*propagator_a2(W.Dot(W))*propagator_rho((P3+P2).Dot(P3+P2));
    return result/sqrt(cross_a2_pi(&ener, &ener));

}


TLorentzVectorC  Cmd3Generator2pi2pi0_ke::H_a2_pi(TLorentzVectorC P, TLorentzVectorC P1, TLorentzVectorC P2, TLorentzVectorC P3, TLorentzVectorC P4){

    TLorentzVectorC res = H_a2_pi_beforepi0sim(P,P1,P2,P3,P4)  + H_a2_pi_beforepi0sim(P,P2,P1,P3,P4);
    return res;
    
}

double Cmd3Generator2pi2pi0_ke::factor_a2pi(){
    return factor_ampl(Energy,28);
}
double Cmd3Generator2pi2pi0_ke::factor_a2pi_ph(){
    return factor_ampl(Energy,30);
}



//----------------------------------------------------------------------------------
// -----------------------------rho sigma-------------------------------------------
//----------------------------------------------------------------------------------

double Cmd3Generator2pi2pi0_ke::cross_rho_sigma(double *e, double *par){

    double norm[] = {0.000876075, 0.00203927, 0.00413222, 0.00768718, 0.0136843, 0.0242913, 0.0427909, 0.0782594, 0.150395, 0.282294, 0.520047, 0.805074, 1.22889, 1.65961, 2.33903, 2.99221, 3.97308, 4.4554, 5.1025, 5.70802, 6.67461, 6.67031, 7.22386, 7.26811, 7.27732, 7.9501, 7.80068, 8.19164};
    double factor;
    if(e[0] <= 0.625)factor = norm[0];
    else if(e[0] >= 1.975)factor = norm[27];
    else {
        int index = (int)((e[0] - 0.625)/0.05);
        double e1 = 0.625 + index*0.05;
        double e2 = 0.625 + (index+1.)*0.05;
        factor = norm[index]*(e2-e[0])/(e2-e1) + norm[index+1]*(-e1+e[0])/(e2-e1);
    }
    return factor;
}

TLorentzVectorC  Cmd3Generator2pi2pi0_ke::H_rho_sigma_beforesim(TLorentzVectorC P, TLorentzVectorC P1, TLorentzVectorC P2, TLorentzVectorC P3, TLorentzVectorC P4){
    TLorentzVectorC W = P3 + P4;
    double ener = P.E().real();
    double ph = 1.5;//2.98311 -2.57009*ener + 0.378539*ener*ener;
    complex<double>  FF = complex<double> (cos(ph),sin(ph));
    TLorentzVectorC result = (P3-P4)*P.Dot(W) - W*P.Dot(P3-P4);
    result = result*propagator_rho(W.Dot(W))*FF*propagator_sigma((P1+P2).Dot(P1+P2));
    return -result/sqrt(cross_rho_sigma(&ener, &ener));
}

TLorentzVectorC  Cmd3Generator2pi2pi0_ke::H_rho_sigma(TLorentzVectorC P, TLorentzVectorC P1, TLorentzVectorC P2, TLorentzVectorC P3, TLorentzVectorC P4){
    TLorentzVectorC res(0.,0.,0.,0.);
    res = res + H_rho_sigma_beforesim(P,P1,P2,P3,P4);
    return res;
}

double Cmd3Generator2pi2pi0_ke::factor_rho_sigma(){
    return factor_ampl(Energy,20);
}

double Cmd3Generator2pi2pi0_ke::factor_rho_sigma_ph(){
    return factor_ampl(Energy,22);
}


//---------------------------------------------------------------------------
// -----------------------------------RHO F0---------------------------------
//---------------------------------------------------------------------------
double Cmd3Generator2pi2pi0_ke::cross_rho_f0(double *e, double *par){

    double norm[] = {0.000120128, 0.0002924, 0.000618387, 0.00117923, 0.00212543, 0.00378536, 0.0067353, 0.0122104, 0.0233578, 0.0438315, 0.0820602, 0.127734, 0.194391, 0.260065, 0.361597, 0.46138, 0.626127, 0.746501, 0.949543, 1.23827, 1.72136, 2.1918, 3.04271, 3.61597, 4.39168, 4.85355, 5.56185, 5.57303};
    double factor;
    if(e[0] <= 0.625)factor = norm[0];
    else if(e[0] >= 1.975)factor = norm[27];
    else {
        int index = (int)((e[0] - 0.625)/0.05);
        double e1 = 0.625 + index*0.05;
        double e2 = 0.625 + (index+1.)*0.05;
        factor = norm[index]*(e2-e[0])/(e2-e1) + norm[index+1]*(-e1+e[0])/(e2-e1);
    }
    return factor;
}

TLorentzVectorC  Cmd3Generator2pi2pi0_ke::H_rho_f0_before_sim(TLorentzVectorC P, TLorentzVectorC P1, TLorentzVectorC P2, TLorentzVectorC P3, TLorentzVectorC P4){
    TLorentzVectorC W = P3 + P4;
    double ener = P.E().real();
    TLorentzVectorC result(0.,0.,0.,0.);
    result = (P3-P4)*P.Dot(W) - W*P.Dot(P3-P4);
    double ph = 0.;
    complex<double>  FF = complex<double> (cos(ph),sin(ph));
    result = result*propagator_rho(W.Dot(W))*propagator_f0((P1+P2).Dot(P1+P2))*FF;
    return result/sqrt(cross_rho_f0(&ener, &ener));
}

TLorentzVectorC  Cmd3Generator2pi2pi0_ke::H_rho_f0(TLorentzVectorC P, TLorentzVectorC P1, TLorentzVectorC P2, TLorentzVectorC P3, TLorentzVectorC P4){
    TLorentzVectorC res(0.,0.,0.,0.);
    res = H_rho_f0_before_sim(P,P1,P2,P3,P4);
    return res;
}


double Cmd3Generator2pi2pi0_ke::factor_rho_f0(){
    return factor_ampl(Energy,12);
}

double Cmd3Generator2pi2pi0_ke::factor_rho_f0_ph(){
    return factor_ampl(Energy,14);
}

//--------------------------------------------------------------------------
// -----------------------------------RHO f0_inter--------------------------
//--------------------------------------------------------------------------
double Cmd3Generator2pi2pi0_ke::cross_rho_f0_inter(double *e, double *par){

    double norm[] = {0.000120128, 0.0002924, 0.000618387, 0.00117923, 0.00212543, 0.00378536, 0.0067353, 0.0122104, 0.0233578, 0.0438315, 0.0820602, 0.127734, 0.194391, 0.260065, 0.361597, 0.46138, 0.626127, 0.746501, 0.949543, 1.23827, 1.72136, 2.1918, 3.04271, 3.61597, 4.39168, 4.85355, 5.56185, 5.57303};
    double factor;
    if(e[0] <= 0.625)factor = norm[0];
    else if(e[0] >= 1.975)factor = norm[27];
    else {
        int index = (int)((e[0] - 0.625)/0.05);
        double e1 = 0.625 + index*0.05;
        double e2 = 0.625 + (index+1.)*0.05;
        factor = norm[index]*(e2-e[0])/(e2-e1) + norm[index+1]*(-e1+e[0])/(e2-e1);
    }
    return factor;
}

TLorentzVectorC  Cmd3Generator2pi2pi0_ke::H_rho_f0_inter_before_sim(TLorentzVectorC P, TLorentzVectorC P1, TLorentzVectorC P2, TLorentzVectorC P3, TLorentzVectorC P4){
    TLorentzVectorC W = P3 + P4;
    double ener = P.E().real();
    TLorentzVectorC result(0.,0.,0.,0.);
    result = (P3-P4)*P.Dot(W) - W*P.Dot(P3-P4);
    double ph = 0.;//-2.33713 + 6.89022*ener -3.23062*ener*ener;
    complex<double>  FF = complex<double> (cos(ph),sin(ph));
    result = result*propagator_rho(W.Dot(W))*propagator_f0_inter((P1+P2).Dot(P1+P2))*FF;
    return -result/sqrt(cross_rho_f0_inter(&ener, &ener));
}

TLorentzVectorC  Cmd3Generator2pi2pi0_ke::H_rho_f0_inter(TLorentzVectorC P, TLorentzVectorC P1, TLorentzVectorC P2, TLorentzVectorC P3, TLorentzVectorC P4){
    TLorentzVectorC res(0.,0.,0.,0.);
    res = H_rho_f0_inter_before_sim(P,P1,P2,P3,P4);
    return res;
}

double Cmd3Generator2pi2pi0_ke::factor_rho_f0_inter(){
    return factor_ampl(Energy,12);
}

double Cmd3Generator2pi2pi0_ke::factor_rho_f0_inter_ph(){
    return factor_ampl(Energy,14);
}


//==================================================
//---------------------rho+ rho- -------------------
//==================================================

double Cmd3Generator2pi2pi0_ke::cross_rhop_rhom(double *e, double *par){

    double norm[] = {1.98223e-07, 1.43368e-06, 6.12057e-06, 2.02618e-05, 5.76425e-05, 0.000145236, 0.000351006, 0.000812736, 0.00177885, 0.00404068, 0.0088983, 0.0193122, 0.0398636, 0.0832108, 0.162197, 0.304643, 0.569157, 1.01843, 1.75296, 2.84571, 4.69834, 6.51253, 9.54061, 12.1267, 15.4348, 19.176, 22.3248, 27.298};
    double factor;
    if(e[0] <= 0.625)factor = norm[0];
    else if(e[0] >= 1.975)factor = norm[27];
    else {
        int index = (int)((e[0] - 0.625)/0.05);
        double e1 = 0.625 + index*0.05;
        double e2 = 0.625 + (index+1.)*0.05;
        factor = norm[index]*(e2-e[0])/(e2-e1) + norm[index+1]*(-e1+e[0])/(e2-e1);
    }
    return factor;
}


TLorentzVectorC  Cmd3Generator2pi2pi0_ke::H_rhop_rhom_beforpi0sim(TLorentzVectorC P, TLorentzVectorC P1, TLorentzVectorC P2, TLorentzVectorC P3, TLorentzVectorC P4){
    TLorentzVectorC rhop = P1 + P3;
    TLorentzVectorC rhom = P2 + P4;
    TLorentzVectorC omega = P1 - P3;
    TLorentzVectorC eee = P2 - P4;
    double ener = P.E().real();
    TLorentzVectorC result(0.,0.,0.,0.);
    result = (P2-P4)*(P1+P3-P2-P4).Dot(P1-P3) + (P1-P3)*(P1+P3-P2-P4).Dot(P2-P4);
    double ph1 = 0.;double ph2 = 0.;
    complex<double>  FF = complex<double> (cos(ph1),sin(ph1))*complex<double> (cos(ph2),sin(ph2));
    result = -result*propagator_rho(rhop.Dot(rhop))*propagator_rho(rhom.Dot(rhom))*FF;
    
    return result/sqrt(cross_rhop_rhom(&ener, &ener));
}

TLorentzVectorC  Cmd3Generator2pi2pi0_ke::H_rhop_rhom(TLorentzVectorC P, TLorentzVectorC P1, TLorentzVectorC P2, TLorentzVectorC P3, TLorentzVectorC P4){

    TLorentzVectorC result = H_rhop_rhom_beforpi0sim(P,P1,P2,P3,P4) + H_rhop_rhom_beforpi0sim(P,P2,P1,P3,P4);
    result = result - H_rhop_rhom_beforpi0sim(P,P1,P2,P4,P3) - H_rhop_rhom_beforpi0sim(P,P2,P1,P4,P3);
    return result;

}


double Cmd3Generator2pi2pi0_ke::factor_rhop_rhom(){
    return factor_ampl(Energy,16);
}

double Cmd3Generator2pi2pi0_ke::factor_rhop_rhom_ph(){
    return factor_ampl(Energy,18);
}


//==========================================
//---------------rho+rho-_II ---------------
//==========================================

double Cmd3Generator2pi2pi0_ke::cross_rhop_rhom_II(double *e, double *par){
    double norm[] = {5.09254e-14, 1.12193e-12, 1.03475e-11, 6.11362e-11, 2.81282e-10, 1.07244e-09, 3.61526e-09, 1.11212e-08, 3.23475e-08, 8.83757e-08, 2.39315e-07, 6.36937e-07, 1.67296e-06, 4.37807e-06, 1.14225e-05, 2.95183e-05, 7.61204e-05, 0.000201503, 0.00052664, 0.00121244, 0.00256362, 0.0049281, 0.00854306, 0.0138288, 0.0221141, 0.0323682, 0.0459657, 0.0659924};
    double factor;
    if(e[0] <= 0.625)factor = norm[0];
    else if(e[0] >= 1.975)factor = norm[27];
    else {
        int index = (int)((e[0] - 0.625)/0.05);
        double e1 = 0.625 + index*0.05;
        double e2 = 0.625 + (index+1.)*0.05;
        factor = norm[index]*(e2-e[0])/(e2-e1) + norm[index+1]*(-e1+e[0])/(e2-e1);
    }
    return factor;
}

TLorentzVectorC  Cmd3Generator2pi2pi0_ke::H_rhop_rhom_II(TLorentzVectorC P, TLorentzVectorC P1, TLorentzVectorC P2, TLorentzVectorC P3, TLorentzVectorC P4){

    TLorentzVectorC rhop = P1 + P3;
    TLorentzVectorC rhom = P2 + P4;
    TLorentzVectorC omega = P1 - P3;
    TLorentzVectorC eee = P2 - P4;
    double ener = P.E().real();
    TLorentzVectorC result(0.,0.,0.,0.);
    result = (P1 + P3 - P2 - P4)*(P2-P4).Dot(P1-P3);
    double ph1 = 0.; double ph2 = 0.;
    complex<double>  FF = complex<double> (cos(ph1),sin(ph1))*complex<double> (cos(ph2),sin(ph2));
    result = result*propagator_rho(rhop.Dot(rhop))*propagator_rho(rhom.Dot(rhom));
    
    return result/sqrt(cross_rhop_rhom_II(&ener, &ener));

}






//==================================================================
// ------------------h1 (rho pi) pi --------------------------------
//==================================================================

double Cmd3Generator2pi2pi0_ke::cross_h1_rhopi(double *e, double *par){
    double norm[] = {1.2186e-05, 4.41931e-05, 0.00013166, 0.000352621, 0.000941573, 0.00232855, 0.00587275, 0.0156255, 0.0420362, 0.117464, 0.298101, 0.656967, 1.31859,  2.26348, 3.4833, 4.9666, 6.33923, 7.45706, 8.57597, 9.34138, 10.635, 10.9817, 11.8845, 11.787, 12.3091, 12.6336, 12.9628, 13.0517};
    double factor;
    if(e[0] <= 0.625)factor = norm[0];
    else if(e[0] >= 1.975)factor = norm[27];
    else {
        int index = (int)((e[0] - 0.625)/0.05);
        double e1 = 0.625 + index*0.05;
        double e2 = 0.625 + (index+1.)*0.05;
        factor = norm[index]*(e2-e[0])/(e2-e1) + norm[index+1]*(-e1+e[0])/(e2-e1);
    }
    return factor;
}

TLorentzVectorC  Cmd3Generator2pi2pi0_ke::H_h1_rhopi_beforepi0sim(TLorentzVectorC P, TLorentzVectorC P1, TLorentzVectorC P2, TLorentzVectorC P3, TLorentzVectorC P4){
    TLorentzVectorC W = P2 + P3 + P4;
    TLorentzVectorC k = P3 + P4;
    TLorentzVectorC result = (P3-P4)*P.Dot(W)*W.Dot(k);
    result = result - W*P.Dot(P3-P4)*W.Dot(k);
    result = result + W*P.Dot(k)*W.Dot(P3-P4);
    result = result - k*W.Dot(P3-P4)*P.Dot(W);
    double ener = P.E().real();
    return result*propagator_rho(k.Dot(k))*propagator_h1(W.Dot(W))/sqrt(cross_h1_rhopi(&ener, &ener));
}

TLorentzVectorC  Cmd3Generator2pi2pi0_ke::H_h1_rhopi_bef(TLorentzVectorC P, TLorentzVectorC P1, TLorentzVectorC P2, TLorentzVectorC P3, TLorentzVectorC P4){
        
    TLorentzVectorC result(0.,0.,0.,0.);
    result = result + H_h1_rhopi_beforepi0sim(P,P1,P4,P3,P2);
    result = result - H_h1_rhopi_beforepi0sim(P,P1,P2,P3,P4);
    result = result + H_h1_rhopi_beforepi0sim(P,P1,P3,P2,P4);
    return result;
}

TLorentzVectorC  Cmd3Generator2pi2pi0_ke::H_h1_rhopi(TLorentzVectorC P, TLorentzVectorC P1, TLorentzVectorC P2, TLorentzVectorC P3, TLorentzVectorC P4){
        
    TLorentzVectorC result(0.,0.,0.,0.);
    result = result + H_h1_rhopi_bef(P,P1,P2,P3,P4);
    result = result + H_h1_rhopi_bef(P,P2,P1,P3,P4);
    return result;
}

double Cmd3Generator2pi2pi0_ke::factor_h1_rhopi(){
    return factor_ampl(Energy,32);
}
double Cmd3Generator2pi2pi0_ke::factor_h1_rhopi_ph(){
    return factor_ampl(Energy,34);
}


// ----------- Cross section function -------------------------------------
double  Cmd3Generator2pi2pi0_ke::cross2pi2pi0(double E){

    if(E < 0.7) return 0.;
    double en[100],cr[100];
    en[0] = 0.7; cr[0] = 0.0005;
    en[1] = 0.8; cr[1] = 0.1;
    en[2] = 0.8400000; cr[2] = 0.2;
    en[3] = 0.8750000; cr[3] = 0.3;
    en[4] = 0.9000000; cr[4] = 0.4;
    en[5] =  0.9250000; cr[5] = 1.184512;
    en[6] =   0.9500001; cr[6] = 2.972132;
    en[7] =  0.9750000; cr[7] =  4.066643;
    en[8] =  1.000000; cr[8] =  5.724451;
    en[9] =  1.025000; cr[9] =  8.523726;
    en[10] =  1.050000; cr[10] =  8.702085;
    en[11] =  1.075000; cr[11] =  11.97795;
    en[12] =  1.100000; cr[12] =  13.26839;
    en[13] = 1.125000; cr[13] =  15.06046;
    en[14] = 1.150000; cr[14] =  16.71859;
    en[15] = 1.175000; cr[15] =  17.34131;
    en[16] = 1.200000; cr[16] =  19.80570;
    en[17] = 1.225000; cr[17] =  20.18546;
    en[18] = 1.250000; cr[18] =  21.54455;
    en[19] = 1.275000; cr[19] =  23.38247;
    en[20] = 1.300000; cr[20] =  25.77932;
    en[21] = 1.325000; cr[21] =  26.29928;
    en[22] = 1.350000; cr[22] =  26.91924;
    en[23] = 1.375000; cr[23] =  28.40883;
    en[24] = 1.400000; cr[24] =  30.18806;
    en[25] = 1.425000; cr[25] =  31.97551;
    en[26] = 1.450000; cr[26] =  32.10345;
    en[27] = 1.475000; cr[27] =  32.26585;
    en[28] = 1.500000; cr[28] =  31.69104;
    en[29] = 1.525000; cr[29] =  29.54573;
    en[30] = 1.550000; cr[30] =  29.83364;
    en[31] = 1.575000; cr[31] =  27.92464;
    en[32] = 1.600000; cr[32] =  27.31240;
    en[33] = 1.625000; cr[33] =  25.89993;
    en[34] = 1.650000; cr[34] =  24.90640;
    en[35] = 1.675000; cr[35] =  23.58823;
    en[36] = 1.700000; cr[36] =  22.66760;
    en[37] = 1.725000; cr[37] =  20.37718;
    en[38] = 1.750000; cr[38] =  19.49763;
    en[39] = 1.775000; cr[39] =  16.94212;
    en[40] = 1.800000; cr[40] =  15.12758;
    en[41] = 1.825000; cr[41] =  12.53257;
    en[42] = 1.850000; cr[42] =  11.96361;
    en[43] = 1.875000; cr[43] =  10.39149;
    en[44] = 1.900000; cr[44] =  9.308997;
    en[45] = 1.925000; cr[45] =  9.284678;
    en[46] = 1.950000; cr[46] =  8.576995;
    en[47] = 1.975000; cr[47] =  8.684048;
    en[48] = 2.000000; cr[48] =  9.254687;
    en[49] = 2.025000; cr[49] =  9.181911;
    if(E > en[49])return cr[49];
    int i = 0;
    while(E > en[i]){i++;}
    if(i < 2){return cr[i];}
    return (E-en[i-1])/(en[i]- en[i-1])*(cr[i]-cr[i-1]) + cr[i-1];
}

double  Cmd3Generator2pi2pi0_ke::cross4pi(double E){
    if(E < 0.6125) return 0.;
    double en[100],cr[100];
    en[0] = 0.6125; cr[0] = 0.02;
    en[1] = 0.6375; cr[1] = 0.04;
    en[2] = 0.6625; cr[2] = 0.02;
    en[3] = 0.6875; cr[3] = 0.01;
    en[4] = 0.7125; cr[4] = 0.02;
    en[5] = 0.7375; cr[5] = 0.03;
    en[6] = 0.7625; cr[6] = 0.05;
    en[7] = 0.7875; cr[7] = 0.11;
    en[8] = 0.8125; cr[8] = 0.11;
    en[9] = 0.8375; cr[9] = 0.12;
    en[10] = 0.8625; cr[10] = 0.17;
    en[11] = 0.8875; cr[11] = 0.26;
    en[12] = 0.9125; cr[12] = 0.33;
    en[13] = 0.9375; cr[13] = 0.57;
    en[14] = 0.9625; cr[14] = 0.71;
    en[15] = 0.9875; cr[15] = 0.89;
    en[16] = 1.0125; cr[16] = 1.2;
    en[17] = 1.0375; cr[17] = 1.61;
    en[18] = 1.0625; cr[18] = 2.17;
    en[19] = 1.0875; cr[19] = 3.29;
    en[20] = 1.1125; cr[20] = 4.49;
    en[21] = 1.1375; cr[21] = 5.95;
    en[22] = 1.1625; cr[22] = 7.37;
    en[23] = 1.1875; cr[23] = 8.84;
    en[24] = 1.2125; cr[24] = 10.79;
    en[25] = 1.2375; cr[25] = 12.62;
    en[26] = 1.2625; cr[26] = 14.56;
    en[27] = 1.2875; cr[27] = 16.39;
    en[28] = 1.3125; cr[28] = 19.06;
    en[29] = 1.3375; cr[29] = 21.14;
    en[30] = 1.3625; cr[30] = 23.37;
    en[31] = 1.3875; cr[31] = 25.76;
    en[32] = 1.4125; cr[32] = 27.53;
    en[33] = 1.4375; cr[33] = 29.95;
    en[34] = 1.4625; cr[34] = 30.32;
    en[35] = 1.4875; cr[35] = 32.04;
    en[36] = 1.5125; cr[36] = 30.98;
    en[37] = 1.5375; cr[37] = 30.11;
    en[38] = 1.5625; cr[38] = 28.26;
    en[39] = 1.5875; cr[39] = 26.81;
    en[40] = 1.6125; cr[40] = 24.66;
    en[41] = 1.6375; cr[41] = 22.69;
    en[42] = 1.6625; cr[42] = 20.95;
    en[43] = 1.6875; cr[43] = 18.78;
    en[44] = 1.7125; cr[44] = 17.25;
    en[45] = 1.7375; cr[45] = 15.33;
    en[46] = 1.7625; cr[46] = 13.37;
    en[47] = 1.7875; cr[47] = 11.61;
    en[48] = 1.8125; cr[48] = 10.23;
    en[49] = 1.8375; cr[49] = 8.87;
    en[50] = 1.8625; cr[50] = 7.67;
    en[51] = 1.8875; cr[51] = 7.29;
    en[52] = 1.9125; cr[52] = 7.17;
    en[53] = 1.9375; cr[53] = 6.93;
    en[54] = 1.9625; cr[54] = 6.54;
    en[55] = 1.9875; cr[55] = 6.04;
    en[56] = 2.0125; cr[56] = 6.18;
    en[57] = 2.0375; cr[57] = 5.66;
    en[58] = 2.0625; cr[58] = 5.68;
    en[59] = 2.0875; cr[59] = 5.34;
    en[60] = 2.1125; cr[60] = 4.92;
    en[61] = 2.1375; cr[61] = 4.83;
    en[62] = 2.1625; cr[62] = 4.59;
    en[63] = 2.1875; cr[63] = 4.28;
    if(E > en[63])return cr[63];
    int i = 0;
    while(E > en[i]){i++;}
    if(i < 2){return cr[i];}
    return (E-en[i-1])/(en[i]- en[i-1])*(cr[i]-cr[i-1]) + cr[i-1];
}




//------------------ Generate events ------------------------------------------
double radiatorW(double s, double x){
    double L = 2.*TMath::Log(sqrt(s)/me);
    //return 2.*alpha/3.141592*(L-1.)*(1.-x+x*x/2.)/x;
    double dzeta3 = 1.2020569;
    double dzeta2 = 1.64493407;
    double beta = 2.*alpha/3.141592*(L-1.);
    double delta2 = (9./8. - 2*dzeta2)*L*L - (45./16.-11./2.*dzeta2-3.*dzeta3)*L-6./5.*dzeta2*dzeta2-9./2.*dzeta3-6.*dzeta2*TMath::Log(2.)+3./8.*dzeta2+57./12.;
    double delta = 1. + alpha/3.141592*(1.5*L+1./3.*3.141592*3.141592-2.) + alpha*alpha/3.141592/3.141592*delta2;
    //cout << delta*beta*pow(x,beta-1.) << " " << beta/2.*(2.-x) << " " << beta*beta/8.*((2.-x)*(3.*TMath::Log(1.-x)-4.*TMath::Log(x))-4.*TMath::Log(1.-x)/x-6.+x) << endl;
    return delta*beta*pow(x,beta-1.)-beta/2.*(2.-x)+beta*beta/8.*((2.-x)*(3.*TMath::Log(1.-x)-4.*TMath::Log(x))-4.*TMath::Log(1.-x)/x-6.+x);
}

double radiatorW_th(double E, double x, double costh){
    double ss = E*E;
    double s = sqrt(1. - costh*costh);
    double c = costh;
    double one = s*s;
    double two = - x*x*s*s*s*s/2./(x*x - 2.*x + 2.);
    double three = - me*me/ss*4.*((1.-2.*x)*s*s-x*x*c*c*c*c)/(x*x - 2.*x + 2.);
    double P = (one + two + three)/pow(s*s+4.*me*me/ss*c*c,2.);
    return P;
}

double radiatorW_th_integral(double E, double x){
    double L = 2.*TMath::Log(E/me);
    return L - 1.;
}

double radiatorW_delta_x_beta(double s, double x){
    double L = 2.*TMath::Log(sqrt(s)/me);
    double dzeta3 = 1.2020569;
    double dzeta2 = 1.64493407;
    double beta = 2.*alpha/3.141592*(L-1.);
    double delta2 = (9./8. - 2*dzeta2)*L*L - (45./16.-11./2.*dzeta2-3.*dzeta3)*L-6./5.*dzeta2*dzeta2-9./2.*dzeta3-6.*dzeta2*TMath::Log(2.)+3./8.*dzeta2+57./12.;
    double delta = 1. + alpha/3.141592*(1.5*L+1./3.*3.141592*3.141592-2.) + alpha*alpha/3.141592/3.141592*delta2;
    return delta*pow(x,beta) - 0.5*beta*(2*x-0.5*x*x) + beta*beta/8.*0.0612457; // actually, it is valid just for x = 0.001
}
//-------------- Visible cross section ---------------------------------
double Cmd3Generator2pi2pi0_ke::cross(double E){
    //cout << E << " " << cross2pi2pi0(E) << endl;
    if(Mode2pi2pi0 == 1)return cross2pi2pi0(E);
    return cross4pi(E);
}

double Cmd3Generator2pi2pi0_ke::visible_cross_section(double E){
    double s = E*E;
    double result = 0.;
    double x = 0.001;
    double stepx = 0.00001;
    result = cross(E)*radiatorW_delta_x_beta(s,x);
    double xA = x;
    double xB = xA;
    while(sqrt(fabs(s*(1.-x - stepx))) > 4.*Mpi && (x+stepx) < 1.){
        xA = xB;
        xB = xA + stepx;
        double prirost = (cross(sqrt(fabs(s*(1.-xA))))*radiatorW(s,xA) + cross(sqrt(fabs(s*(1.-0.5*xA-0.5*xB))))*radiatorW(s,0.5*xA + 0.5*xB)*4. + cross(sqrt(fabs(s*(1.-xB))))*radiatorW(s,xB))/6.*stepx;
        stepx = stepx/prirost*result*0.001;
        if(sqrt(fabs(s*(1.-x - stepx))) < 4.*Mpi)break;
        result = result + prirost;
        x = x + stepx;
    }
    return result;
}

// -----------  generate x ---------------------------------------
double Cmd3Generator2pi2pi0_ke::generate_photon_energy(double E){
    double radiatorW_delta_x_beta(double, double);
    double radiatorW(double, double);
    double xnew = gRandom->Rndm();
    double s = E*E;
    double Eprime = E;
    double crvis_x = 0.;
    double crvis = visible_cr_sec;
    double x = 0.001;
    double stepx = 0.00001;
    crvis_x = cross(E)*radiatorW_delta_x_beta(s,x);
    //cout << "cross(E) = " << cross(E) << "    radiatorW_delta_x_beta(s,x) = " << radiatorW_delta_x_beta(s,x) << "   visible_cross_section(E) = " << visible_cross_section(E) << endl;
    if(xnew < crvis_x/crvis){return 0.;}
    double xA = x;
    double xB = xA;
    while(xnew >= crvis_x/crvis){
        xA = xB;
        xB = xA + stepx;
        double prirost = (cross(sqrt(fabs(s*(1.-xA))))*radiatorW(s,xA) + cross(sqrt(fabs(s*(1.-0.5*xA-0.5*xB))))*radiatorW(s,0.5*xA + 0.5*xB)*4. + cross(sqrt(fabs(s*(1.-xB))))*radiatorW(s,xB))/6.*stepx;
        stepx = stepx/prirost*crvis_x*0.001;
        crvis_x = crvis_x + prirost;
        x = x + stepx;
        if(sqrt(fabs(s*(1.-x))) < (4.*Mpi + 0.03)){x = x - stepx;break;}
        //cout<<"x,step,cross ="<<","<<x<<", "<<stepx<< " , " << crvis_x <<endl;
    }
    double xmax = (s-16.*0.140*0.140)/s;
    if(x >= xmax)x = xmax;
    return x*E/2.;
}


double Cmd3Generator2pi2pi0_ke::generate_photon_angle(double E, double x){

    double xnew = gRandom->Rndm();
    xnew = xnew*radiatorW_th_integral(E, x);
    double result = 0.;
    int nrun = 10000000;
    double step = (double)2./nrun;
    double costh = -1.;
    int i = 0;
    double xA = costh;
    double xB = xA;
    while(result < xnew){
        xA = xB;
        xB = xA + step;
        costh = costh + step;
        result = result + 2.*(radiatorW_th(E,x,xA)+4.*radiatorW_th(E,x,0.5*xA+0.5*xB)+radiatorW_th(E,x,xB))/6.*step;
        if(fabs(costh) < 0.996)step = (double)2./1000;
        i++;
        //cout << result << " " << costh << " " << step << endl;
    }
    xnew = gRandom->Rndm();
    if(xnew > 0.5){costh = -costh;}
    //cout << "costh = " << costh << endl;
    return costh;
}

// -----------------------------------------------------------------------
TLorentzVectorC Cmd3Generator2pi2pi0_ke::matrix(TLorentzVectorC P, TLorentzVectorC P1, TLorentzVectorC P2, TLorentzVectorC P3, TLorentzVectorC P4){

    TLorentzVectorC H = H_omega_pi0(P,P1,P2,P3,P4)*omegapi0_*complex<double>(cos(omegapi0phase_),sin(omegapi0phase_));
    complex<double> fact = a1pi_rhopi_*complex<double>(cos(a1pi_rhopi_phase),sin(a1pi_rhopi_phase))/(1.+a1pi_sigmapi_);
    if(a1pi_rhopi_!=0.)H = H + (H_a1_rhopi_pi(P,P1,P2,P3,P4) + H_a1_sigmapi_pi(P,P1,P2,P3,P4)*a1pi_sigmapi_*complex<double>(cos(a1pi_sigmapi_phase),sin(a1pi_sigmapi_phase)))*fact;
    if(rhof0_!=0.)H = H + H_rho_f0(P,P1,P2,P3,P4)*rhof0_*complex<double>(cos(rhof0phase_),sin(rhof0phase_));
    if(rhof0_inter!=0.)H = H + H_rho_f0_inter(P,P1,P2,P3,P4)*rhof0_inter*complex<double>(cos(rhof0phase_inter),sin(rhof0phase_inter));
    if(rhosigma_!=0.)H = H + H_rho_sigma(P,P1,P2,P3,P4)*rhosigma_*complex<double>(cos(rhosigmaphase_),sin(rhosigmaphase_));
    if(rhoprhom_II!=0.)H = H + (H_rhop_rhom(P,P1,P2,P3,P4)*rhoprhom_*complex<double>(cos(rhoprhom_phase),sin(rhoprhom_phase))+H_rhop_rhom_II(P,P1,P2,P3,P4))/(1.+rhoprhom_)*rhoprhom_II*complex<double>(cos(rhoprhom_phase_II),sin(rhoprhom_phase_II));
    if(a2pi_!=0.)H = H + H_a2_pi(P,P1,P2,P3,P4)*a2pi_*complex<double>(cos(a2pi_ph),sin(a2pi_ph));
    if(phsp_!=0.)H = H + H_ph_sp(P,P1,P2,P3,P4)*phsp_*complex<double>(cos(phspphase_),sin(phspphase_));
    if(h1pi!=0.)H = H + H_h1_rhopi(P,P1,P2,P3,P4)*h1pi*complex<double>(cos(h1pi_ph),sin(h1pi_ph));
    if(a1pi_rhopi_II!=0.)H = H + H_a1_rhopi_pi_II(P,P1,P2,P3,P4)*a1pi_rhopi_II*complex<double>(cos(a1pi_rhopi_phase_II),sin(a1pi_rhopi_phase_II));
    if(rhof2!=0.)H = H + H_rho_f2(P,P1,P2,P3,P4)*rhof2*complex<double>(cos(rhof2_phase),sin(rhof2_phase));
    
    return H;

}


double Cmd3Generator2pi2pi0_ke::matrix_squared(TLorentzVectorC P, TLorentzVectorC P1, TLorentzVectorC P2, TLorentzVectorC P3, TLorentzVectorC P4){

    TLorentzVectorC H = matrix(P,P1,P2,P3,P4); 

    return H.X().real()*H.X().real() + H.X().imag()*H.X().imag() + H.Y().real()*H.Y().real() + H.Y().imag()*H.Y().imag();

}

double Cmd3Generator2pi2pi0_ke::matrix_squared_c(TLorentzVectorC P, TLorentzVectorC P1, TLorentzVectorC P2, TLorentzVectorC P3, TLorentzVectorC P4){

    TLorentzVectorC H = matrix(P,P1,P2,P3,P4) + matrix(P,P3,P2,P1,P4); 
    H = H + matrix(P,P1,P4,P3,P2) + matrix(P,P3,P4,P1,P2);

    return H.X().real()*H.X().real() + H.X().imag()*H.X().imag() + H.Y().real()*H.Y().real() + H.Y().imag()*H.Y().imag();
}

double Cmd3Generator2pi2pi0_ke::matrix_squared1(TLorentzVector P, TLorentzVector P1, TLorentzVector P2, TLorentzVector P3, TLorentzVector P4){

    if(isPhsp == 1)return 1;
    TLorentzVectorC Pc(complex<double>(P.X(),0.),complex<double>(P.Y(),0.),complex<double>(P.Z(),0.),complex<double>(P.E(),0.));
    TLorentzVectorC P1c(complex<double>(P1.X(),0.),complex<double>(P1.Y(),0.),complex<double>(P1.Z(),0.),complex<double>(P1.E(),0.));
    TLorentzVectorC P2c(complex<double>(P2.X(),0.),complex<double>(P2.Y(),0.),complex<double>(P2.Z(),0.),complex<double>(P2.E(),0.));
    TLorentzVectorC P3c(complex<double>(P3.X(),0.),complex<double>(P3.Y(),0.),complex<double>(P3.Z(),0.),complex<double>(P3.E(),0.));
    TLorentzVectorC P4c(complex<double>(P4.X(),0.),complex<double>(P4.Y(),0.),complex<double>(P4.Z(),0.),complex<double>(P4.E(),0.));

    return matrix_squared(Pc,P1c,P2c,P3c,P4c);
}

double Cmd3Generator2pi2pi0_ke::matrix_squared1_c(TLorentzVector P, TLorentzVector P1, TLorentzVector P2, TLorentzVector P3, TLorentzVector P4){

    if(isPhsp == 1)return 1;
    TLorentzVectorC Pc(complex<double>(P.X(),0.),complex<double>(P.Y(),0.),complex<double>(P.Z(),0.),complex<double>(P.E(),0.));
    TLorentzVectorC P1c(complex<double>(P1.X(),0.),complex<double>(P1.Y(),0.),complex<double>(P1.Z(),0.),complex<double>(P1.E(),0.));
    TLorentzVectorC P2c(complex<double>(P2.X(),0.),complex<double>(P2.Y(),0.),complex<double>(P2.Z(),0.),complex<double>(P2.E(),0.));
    TLorentzVectorC P3c(complex<double>(P3.X(),0.),complex<double>(P3.Y(),0.),complex<double>(P3.Z(),0.),complex<double>(P3.E(),0.));
    TLorentzVectorC P4c(complex<double>(P4.X(),0.),complex<double>(P4.Y(),0.),complex<double>(P4.Z(),0.),complex<double>(P4.E(),0.));

    return matrix_squared_c(Pc,P1c,P2c,P3c,P4c);
}

//**************************************************************************************************************************************************
void  Cmd3Generator2pi2pi0_ke::find_majoranta(double E){
    int i = 0;
    TGenPhaseSpace event;
    TLorentzVector P(0.,0.,0., E);
    TLorentzVector* P1 = new TLorentzVector();
    TLorentzVector* P2 = new TLorentzVector();
    TLorentzVector* P3 = new TLorentzVector();
    TLorentzVector* P4 = new TLorentzVector();
    double masses[4] = {mpiz, mpiz, Mpi, Mpi};
    if(Mode2pi2pi0 == 0)masses[0] = Mpi;
    if(Mode2pi2pi0 == 0)masses[1] = Mpi;
    event.SetDecay(P, 4, masses);
    majoranta = 0.;
    while(i < 1000000){
        double mes = event.Generate();
        P1 = event.GetDecay(0);
        P2 = event.GetDecay(1);
        P3 = event.GetDecay(2);
        P4 = event.GetDecay(3);
        double matrix_squar = -1;
        if(Mode2pi2pi0 == 1)matrix_squar = matrix_squared1(P,*P1,*P2,*P3,*P4);
        if(Mode2pi2pi0 == 0)matrix_squar = matrix_squared1_c(P,*P1,*P2,*P3,*P4);
        if(mes*matrix_squar > majoranta){majoranta = mes*matrix_squar;}
        i++;
    }
    P1->Delete();
    P2->Delete();
    P3->Delete();
    P4->Delete();
    majoranta = majoranta*1.5;
    cout << "Majoranta =  "<< majoranta << endl;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------


void  Cmd3Generator2pi2pi0_ke::Start_Of_Run(CmdParam& p){
  Energy = p.GetDouble("generator.Energy");
  isPhsp = p.GetInt("generator.isPhsp");
  isISR  = p.GetInt("generator.isISR");
  Mode2pi2pi0 =  p.GetInt("generator.Mode2pi2pi0");
  try{results = p.GetString("generator.results");}catch(...){}
  fill_parameters();
  find_majoranta(Energy);
  visible_cr_sec = visible_cross_section(Energy);
  cout << "======================= Generator parameters =================================" << endl;
  cout << "Energy = " << Energy << endl;
  cout << "isISR = " << isISR << endl; 
  cout << "isPhsp = " << isPhsp << endl;
  cout << "Mode2pi2pi0 = " << Mode2pi2pi0 << endl;    
  cout << "visible_cr_sec = " << visible_cr_sec << endl;
  cout << "amplitudes data from " << results << endl;
  print();
  cout << "============================================================" << endl;

}

HepMC::GenEvent* Cmd3Generator2pi2pi0_ke::Event(int NumEvt)
{
   double xnew = gRandom->Rndm();
   TGenPhaseSpace event;
   TLorentzVector* P1 = new TLorentzVector();
   TLorentzVector* P2 = new TLorentzVector();
   TLorentzVector* P3 = new TLorentzVector();
   TLorentzVector* P4 = new TLorentzVector();
   TLorentzVector* Ph = new TLorentzVector();
   double masses[4] = {mpiz, mpiz, Mpi, Mpi};
   if(Mode2pi2pi0 == 0)masses[0] = Mpi;
   if(Mode2pi2pi0 == 0)masses[1] = Mpi;
   int Pi1ID = 111;
   int Pi2ID = 111;
   if(Mode2pi2pi0 == 0)Pi1ID = 211;
   if(Mode2pi2pi0 == 0)Pi2ID = -211;

   HepMC::GenEvent* evt = new HepMC::GenEvent(0, NumEvt);
   HepMC::GenVertex* vertex = new HepMC::GenVertex(HepLorentzVector(0,0,0,0),0);
   evt->add_vertex(vertex);

   int proba = 0;
   double ephoton = 0;
   double costhetaphoton = 0.;
   while(proba < 50000){       
       if(isISR == 1)ephoton = generate_photon_energy(Energy);
       if(isISR == 1)costhetaphoton = generate_photon_angle(Energy, 2.*ephoton/Energy);
       double sinthetaphoton = sqrt(1. - costhetaphoton*costhetaphoton);
       Ph = new TLorentzVector(ephoton*sinthetaphoton*cos(xnew*2.*3.14),ephoton*sinthetaphoton*sin(xnew*2.*3.14),ephoton*costhetaphoton,ephoton);
       TLorentzVector P(-Ph->X(),-Ph->Y(),-Ph->Z(),Energy - ephoton);
       event.SetDecay(P, 4, masses);
       double mes = event.Generate();
       P1 = event.GetDecay(0);
       P2 = event.GetDecay(1);
       P3 = event.GetDecay(2);
       P4 = event.GetDecay(3);
       xnew = gRandom->Rndm();
       if(proba > 40000.) cout << proba << " ERROR: proba is more than 40000!";
       double matrix_squar = -1;
       if(Mode2pi2pi0 == 1)matrix_squar = matrix_squared1(P,*P1,*P2,*P3,*P4);
       if(Mode2pi2pi0 == 0)matrix_squar = matrix_squared1_c(P,*P1,*P2,*P3,*P4);
       if(mes*matrix_squar < xnew*majoranta){proba++;continue;}
       break;
       proba++;
   }


    HepLorentzVector _P1(P1->X(),P1->Y(),P1->Z(),P1->E());
    HepLorentzVector _P2(P2->X(),P2->Y(),P2->Z(),P2->E());
    HepLorentzVector _P3(P3->X(),P3->Y(),P3->Z(),P3->E());
    HepLorentzVector _P4(P4->X(),P4->Y(),P4->Z(),P4->E());
    HepLorentzVector _Ph(Ph->X(),Ph->Y(),Ph->Z(),Ph->E());

   HepMC::GenParticle* particle1=new HepMC::GenParticle(_P1,Pi1ID,1);
   vertex->add_particle_out(particle1);
   HepMC::GenParticle* particle2=new HepMC::GenParticle(_P2,Pi2ID,1);
   vertex->add_particle_out(particle2);
   HepMC::GenParticle* particle3=new HepMC::GenParticle(_P3,211,1);
   vertex->add_particle_out(particle3);
   HepMC::GenParticle* particle4=new HepMC::GenParticle(_P4,-211,1);
   vertex->add_particle_out(particle4);
   HepMC::GenParticle* particleph=new HepMC::GenParticle(_Ph,22,1);
   vertex->add_particle_out(particleph);

   return evt;
}
