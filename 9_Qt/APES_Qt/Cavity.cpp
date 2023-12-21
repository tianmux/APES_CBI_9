#include "Beam.h"
#include "Cavity.h"
#include <iostream>
#include <vector>
#include <iomanip>
#include <omp.h>

// Define the imaginary unit 'i'
// const std::complex<double> ii(0.0, 1.0);

Cavity::Cavity(const InputData& inputData, const Beam& beam) {
    
    this->cavityProps = inputData.cavityProps;
    this->iTurn = 0;
    this->vGen.resize(beam.nPar*beam.nBunch);
    this->vBeam.resize(beam.nPar*beam.nBunch);
    this->time.resize(beam.nPar*beam.nBunch);
    this->n_feed = cavityProps.h/cavityProps.feed_step;
    std::cout<<"Feed step: "<<cavityProps.feed_step<<std::endl;
    std::cout<<"h : "<<cavityProps.h<<std::endl;
    std::cout<<"n_feed: "<<n_feed<<std::endl;
    std::vector<int> tmp_feed_bucket(n_feed);
    this->feed_time.resize(n_feed);
    this->vRef.resize(n_feed);
    this->n_thread = inputData.generalProps.n_thread;
    ///////////////////////////////////////////////////////
    // Initialize the feed time
    std::cout<<"initializing feed time"<<std::endl;
    
    for(int i = 0; i < n_feed; ++i){
        tmp_feed_bucket[i] = int((i+1)*cavityProps.feed_step);
        // std::cout << "tmp_feed_bucket[i]: " << tmp_feed_bucket[i] << std::endl;
        // Shift the feedback time by 1/4 of the RF period relative to the centroids of bunches if they are in same RF bucket.
        this->feed_time[i] = (inputData.generalProps.t0)*this->cavityProps.TRF+this->cavityProps.TRF*tmp_feed_bucket[i]+this->cavityProps.TRF/4.0*1;
    }
    tmp_feed_bucket[n_feed-1] = cavityProps.h-1;
    this->feed_time[n_feed-1] = (inputData.generalProps.t0)*this->cavityProps.TRF+tmp_feed_bucket[n_feed-1]*this->cavityProps.TRF+this->cavityProps.TRF/4.0*1;
    ///////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////
    // Initialize the reference voltage of feedback
    std::cout<<"initializing reference voltage"<<std::endl;
    std::complex<double> tmp_Vref(cavityProps.Vsynch, cavityProps.Vquad);
    std::complex<double> tmp_phi = std::polar(1.0,inputData.pi/2.0);
    for(int i = 0; i < n_feed; ++i){
        this->vRef[i] = tmp_Vref*tmp_phi;
    }
    ///////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////
    // Initialize the generator current, 
    // and the initial value of Vg Ug and Vb.
    std::cout<<"initializing generator current"<<std::endl;
    // Vb = Vadd*(1/(1-np.exp(1j*cavs[0].wc*Trf*fillStep)*np.exp(-cavs[0].alpha*Trf*fillStep))-1/2)
    this->Vb = -beam.qPb*this->cavityProps.wCav*this->cavityProps.RoQ*(1.0/(1.0-std::exp(ii*this->cavityProps.wCav*this->cavityProps.TRF*double(beam.fill_step))*std::exp(-this->cavityProps.alpha*this->cavityProps.TRF*beam.fill_step))-1.0/2.0);
    // Ig = (Vc/Z-Vbphasor/Z)
    this->Ig = (this->cavityProps.Vc/this->cavityProps.Z-this->Vb/this->cavityProps.Z);
    this->Vg = this->Ig*this->cavityProps.Z;    
    this->Vgm1 = this->Vg;
    this->Vgp = this->Vg*ii*this->cavityProps.wRF;
    std::cout<<"Vg :"<<this->Vg<<std::endl;
    std::cout<<"Vgp: "<<this->Vgp<<std::endl;
    std::cout<<"1st in Ug: "<< this->cavityProps.L*this->Ig*std::exp(ii*this->cavityProps.wRF*inputData.generalProps.t0)<<std::endl;
    std::cout<<"2nd in Ug: "<< this->cavityProps.L/this->cavityProps.R*this->Vg<<std::endl;
    std::cout<<"3rd in Ug: "<< this->cavityProps.L*this->cavityProps.C*this->Vgp<<std::endl;
    this->Ug = this->cavityProps.L*(this->Ig*std::exp(ii*this->cavityProps.wRF*inputData.generalProps.t0)-1/this->cavityProps.R*this->Vg-this->cavityProps.C*this->Vgp);
    this->Ugm1 = this->Ug;
    std::cout<<"Ug :"<<this->Ug<<std::endl;
    std::cout<<"Vb: "<<this->Vb<<std::endl;
    std::cout<<"Vtot : "<<this->Vg+this->Vb<<std::endl;
    this->time = beam.time;
    this->tlast_feed = inputData.generalProps.t0*this->cavityProps.TRF;
    this->tlast_Vb = inputData.generalProps.t0*this->cavityProps.TRF;
}

void Cavity::printCavityProps() const{
    // Print CavityProperties
    std::cout << "Cavity Properties:" << std::endl;
    std::cout << "  Num Cavities: " << cavityProps.numCavities << std::endl;
    std::cout << "  Harmonic Number: " << cavityProps.h << std::endl;
    std::cout << "  R/Q: " << cavityProps.RoQ << std::endl;
    std::cout << "  Loaded Q: " << cavityProps.QL << std::endl;
    std::cout << "  Synchronous Voltage: " << cavityProps.Vsynch << std::endl;
    std::cout << "  Quadrature Voltage: " << cavityProps.Vquad << std::endl;
    std::cout << "  Frequency Shift: " << cavityProps.df << std::endl;
    std::cout << "  Cavity Frequency: " << cavityProps.wCav << std::endl;
    std::cout << "  RF Frequency: " << cavityProps.fRF << std::endl;
    std::cout << "  RF Period: " << cavityProps.TRF << std::endl;
    std::cout << "  Cavity Period: " << cavityProps.TCav << std::endl;
    std::cout << "  Cavity Inductance: " << cavityProps.L << std::endl;
    std::cout << "  Cavity Capacitance: " << cavityProps.C << std::endl;
    std::cout << "  Cavity Resistance: " << cavityProps.R << std::endl;
    std::cout << "  Cavity Damping Factor: " << cavityProps.alpha << std::endl;
    std::cout << "  Cavity Impedance: " << cavityProps.Z << std::endl;
    std::cout << "  Cavity Voltage: " << cavityProps.Vc << std::endl;
    std::cout << "  Generator Current: " << cavityProps.Ig << std::endl;
    std::cout << "  Beam DC Current: " << cavityProps.IbDC << std::endl;
    std::cout << "  Beam Reference Voltage: " << cavityProps.Vb_ref << std::endl;

    std::cout << "  Feedback Step: " << cavityProps.feed_step << std::endl;
    std::cout << "  Proportional Gain: " << cavityProps.gp << std::endl;
    std::cout << "  Integral Gain: " << cavityProps.gi << std::endl;
    std::cout << "  Derivative Gain: " << cavityProps.gd << std::endl;
    std::cout << "  Comb Gain: " << cavityProps.gc << std::endl;
    std::cout << "  Epsilon: " << cavityProps.epsilon << std::endl;
    std::cout << "  Feedback Delay: " << cavityProps.delay << std::endl;
    std::cout << "  Synchronous Phase: " << cavityProps.phis << std::endl;
    std::cout << "  Number of Feedback Time Points: " << n_feed << std::endl;
    std::cout << "  Number of Threads: " << n_thread << std::endl;

    std::cout << "  Feed Time: " << std::endl;
    for(int i = 0; i < feed_time.size(); ++i){
        std::cout << "    " << feed_time[i] << std::endl;
    }
    std::cout << "  Feed Voltage: " << std::endl;
    for(int i = 0; i < vRef.size(); ++i){
        std::cout << "    " << vRef[i] << std::endl;
    }

    std::cout << "  Vgen: " << std::endl;
    for(int i = 0; i < vGen.size(); ++i){
        std::cout << "    " << vGen[i] << std::endl;
    }
    std::cout << "  Vbeam: " << std::endl;
    for(int i = 0; i < vBeam.size(); ++i){
        std::cout << "    " << vBeam[i] << std::endl;
    }
    std::cout << "  Time: " << std::endl;
    for(int i = 0; i < time.size(); ++i){
        std::cout << "    " << time[i] << std::endl;
    }
    std::cout << "  Last Feedback Time: " << tlast_feed << std::endl;
    std::cout << "  Last Beam Voltage Update Time: " << tlast_Vb << std::endl;
    std::cout << "  Ig: " << Ig << std::endl;
    std::cout << "  Vg: " << Vg << std::endl;
    std::cout << "  Vgp: " << Vgp << std::endl;
    std::cout << "  Ug: " << Ug << std::endl;
    std::cout << "  Vb: " << Vb << std::endl;
    std::cout << "  Vgm1: " << Vgm1 << std::endl;
    std::cout << "  Ugm1: " << Ugm1 << std::endl;
    std::cout << "  Vtot: " << Vg+Vb << std::endl;
}

void Cavity::calVgen(Beam& beam, int& istart, int& iend, int& ilast_feed){
    // beam is the beam object
    // istart is the index of the first time point (of the particle) we need to calculate the generator voltage after last feedback time;
    // iend is the index of the last time point (of the particle) we need to calculate the generator voltage after last feedback time;
    // ilast_feed is the index of the last time point when we updated the generator current;
    // 
    // calculate the generated voltage at each time point.
    // But do not update the Vgm and Ugm for now. 
    /*
    std::cout<<"Calculating generator voltage"<<std::endl;
    std::cout<<"i Last feed: "<<ilast_feed<<std::endl;
    std::cout<<"Last feedback time: "<<this->tlast_feed<<std::endl;
    */
    std::vector<double> tmp_dt(iend-istart);
    double tmp_last_t = this->tlast_feed;
    for(int i = istart; i < iend; ++i){
        tmp_dt[i-istart] = this->time[i]-tmp_last_t;
    }
    //std::cout<<"Initial value of the generator voltage: "<<this->Vg<<std::endl;
    std::complex<double> A(this->cavityProps.alpha, this->cavityProps.wRF+this->cavityProps.wCav);
    std::complex<double> B(this->cavityProps.alpha, this->cavityProps.wRF-this->cavityProps.wCav);

    // This part can be parallelized since every particle is independent.
    // set the number of thread
    omp_set_num_threads(this->n_thread);
    #pragma omp parallel for
    for(int i = istart; i < iend; ++i){
        /*
        std::complex<double> tmp_exp_cav =std::exp(ii*this->cavityProps.wCav*tmp_dt[i-istart]);
        std::complex<double> tmp_exp_rf = std::exp(ii*this->cavityProps.wRF*tmp_dt[i-istart]);
        std::complex<double> tmp_exp_cav_1 = std::exp(-ii*this->cavityProps.wCav*tmp_dt[i-istart]);
        std::complex<double> tmp_exp_alpha = std::exp(-this->cavityProps.alpha*tmp_dt[i-istart]);
        
        this->vGen[i] = this->Ig/this->cavityProps.C*
                        (ii*this->cavityProps.wRF/A/B*tmp_exp_rf+
                        (this->cavityProps.alpha+ii*this->cavityProps.wCav)/(-2.0*ii*this->cavityProps.wCav*A)*tmp_exp_cav_1*tmp_exp_alpha+
                        (this->cavityProps.alpha-ii*this->cavityProps.wCav)/(2.0*ii*this->cavityProps.wCav*B)*tmp_exp_cav*tmp_exp_alpha)+
                        (this->Vgm1*std::cos(this->cavityProps.wRF*tmp_dt[i-istart])-this->cavityProps.alpha/this->cavityProps.wCav*std::sin(this->cavityProps.wRF*tmp_dt[i-istart])-
                        this->Ugm1/this->cavityProps.C/this->cavityProps.L/(this->cavityProps.wCav)*std::sin(this->cavityProps.wCav*tmp_dt[i-istart]))*tmp_exp_alpha;
        */
        
        this->vGen[i] = this->Ig/this->cavityProps.C*
                        (ii*this->cavityProps.wRF/A/B*std::exp(ii*this->cavityProps.wRF*tmp_dt[i-istart])+
                        (this->cavityProps.alpha+ii*this->cavityProps.wCav)/(-2.0*ii*this->cavityProps.wCav*A)*std::exp(-(ii*this->cavityProps.wCav+this->cavityProps.alpha)*tmp_dt[i-istart])+
                        (this->cavityProps.alpha-ii*this->cavityProps.wCav)/(2.0*ii*this->cavityProps.wCav*B)*std::exp((ii*this->cavityProps.wCav-this->cavityProps.alpha)*tmp_dt[i-istart]))+
                        (this->Vgm1*std::cos(this->cavityProps.wRF*tmp_dt[i-istart])-this->cavityProps.alpha/this->cavityProps.wCav*std::sin(this->cavityProps.wRF*tmp_dt[i-istart])-
                        this->Ugm1/this->cavityProps.C/this->cavityProps.L/(this->cavityProps.wCav)*std::sin(this->cavityProps.wCav*tmp_dt[i-istart]))*std::exp(-this->cavityProps.alpha*tmp_dt[i-istart]);
        
        //std::cout<<"Vgen at time "<<this->time[i]<<" : "<<this->vGen[i]<<std::endl;
    }
}

void Cavity::calVbeam(Beam& beam, int& istart, int& iend, int& ilast_feed){
    /*
    std::cout<<"Calculating beam voltage"<<std::endl;
    std::cout<<"i Last feed: "<<ilast_feed<<std::endl;
    std::cout<<"Last feedback time: "<<this->tlast_Vb<<std::endl;
    */
    std::vector<double> tmp_dt(iend-istart);
    
    double tmp_last_t = this->tlast_Vb;
    for(int i = istart; i < iend; ++i){
        tmp_dt[i-istart] = this->time[i]-tmp_last_t;
    }
    // This part can't be parallelized since every particle needs to know the voltage of all previous particles. 
    // first particle sees the roated and decayed Vb from last feedback point.
    /*
    std::cout<<"Initial value of the beam voltage: "<<this->Vb<<std::endl;
    std::cout<<"last time Vb was update:"<<tmp_last_t<<std::endl;
    std::cout<<"Time point of interest:"<<this->time[istart]<<std::endl;
    */
    double dt = this->time[istart]-tmp_last_t;
    //std::cout<<"dt: "<<std::scientific<<std::setprecision(15)<<dt<<std::endl;
    this->vBeam[istart] = this->Vb*std::exp((ii*this->cavityProps.wCav-this->cavityProps.alpha)*(this->time[istart]-tmp_last_t))+beam.vAdd[istart];
    //std::cout<<"Vb initial "<<std::exp((ii*this->cavityProps.wCav-this->cavityProps.alpha)*0.0)<<std::endl;
    //std::cout<<"Vb at time "<<this->time[istart]<<" : "<<this->vBeam[istart]<<std::endl;
    // The rest of the particales in this batch sees the rotated and decayed Vb from last particle infront of it.
    for(int i = istart+1; i < iend; ++i){
        this->vBeam[i] = this->vBeam[i-1]*std::exp((ii*this->cavityProps.wCav-this->cavityProps.alpha)*(this->time[i]-this->time[i-1]))+beam.vAdd[i];
        //std::cout<<"Vb at time "<<this->time[i]<<" : "<<this->vBeam[i]<<std::endl;
    }
    this->tlast_Vb = this->time[iend-1];
    this->Vb = this->vBeam[iend-1];
    //std::cout<<"Vb at the end of the batch: "<<this->Vb<<std::endl;
}

void Cavity::updateVadd(Beam& beam){
    //std::cout<<"Updating Vadd"<<std::endl;
    // Set the number of threads
    omp_set_num_threads(this->n_thread);
    #pragma omp parallel for
    for(int i = 0;i<beam.vAdd.size();++i){
        beam.vAdd[i] = -beam.q[i]*this->cavityProps.wCav*this->cavityProps.RoQ;
    }
    //std::cout<<"Vadd updated"<<std::endl;
}
void Cavity::updateVg(const int& ifeed){
    //std::cout<<"Updating Vg"<<std::endl;
    //std::cout<<"Vg before update: "<<this->Vg<<std::endl;    

    double tmp_dt = this->feed_time[ifeed]-this->tlast_feed;    
    //std::cout<<"dphi : "<<this->cavityProps.wRF*tmp_dt/2.0/M_PI<<std::endl;

    std::complex<double> A(this->cavityProps.alpha, this->cavityProps.wRF+this->cavityProps.wCav);
    std::complex<double> B(this->cavityProps.alpha, this->cavityProps.wRF-this->cavityProps.wCav);
    this->Vg = this->Ig/this->cavityProps.C*
                        (ii*this->cavityProps.wRF/A/B*std::exp(ii*this->cavityProps.wRF*tmp_dt)+
                        (this->cavityProps.alpha+ii*this->cavityProps.wCav)/(-2.0*ii*this->cavityProps.wCav*A)*std::exp(-(ii*this->cavityProps.wCav+this->cavityProps.alpha)*tmp_dt)+
                        (this->cavityProps.alpha-ii*this->cavityProps.wCav)/(2.0*ii*this->cavityProps.wCav*B)*std::exp((ii*this->cavityProps.wCav-this->cavityProps.alpha)*tmp_dt))+
                        (this->Vgm1*std::cos(this->cavityProps.wCav*tmp_dt)-this->cavityProps.alpha/this->cavityProps.wCav*std::sin(this->cavityProps.wCav*tmp_dt)-
                        this->Ugm1/this->cavityProps.C/this->cavityProps.L/(this->cavityProps.wCav)*std::sin(this->cavityProps.wCav*tmp_dt))*std::exp(-this->cavityProps.alpha*tmp_dt);
    /*
    std::cout<<"1st part : "<<this->Ig/this->cavityProps.C*
                        (ii*this->cavityProps.wRF/A/B*std::exp(ii*this->cavityProps.wRF*tmp_dt))<<std::endl;
    std::cout<<"2nd part : "<<this->Ig/this->cavityProps.C*
                        (this->cavityProps.alpha+ii*this->cavityProps.wCav)/(-2.0*ii*this->cavityProps.wCav*A)*std::exp(-(ii*this->cavityProps.wCav+this->cavityProps.alpha)*tmp_dt)<<std::endl;
    std::cout<<"3rd part : "<<this->Ig/this->cavityProps.C*
                        (this->cavityProps.alpha-ii*this->cavityProps.wCav)/(2.0*ii*this->cavityProps.wCav*B)*std::exp((ii*this->cavityProps.wCav-this->cavityProps.alpha)*tmp_dt)<<std::endl;
    std::cout<<"4th part : "<<(this->Vgm1*std::cos(this->cavityProps.wCav*tmp_dt)-this->cavityProps.alpha/this->cavityProps.wCav*std::sin(this->cavityProps.wCav*tmp_dt))*std::exp(-this->cavityProps.alpha*tmp_dt)<<std::endl;
    std::cout<<"5th part : "<<(-this->Ugm1/this->cavityProps.C/this->cavityProps.L/(this->cavityProps.wCav)*std::sin(this->cavityProps.wCav*tmp_dt))*std::exp(-this->cavityProps.alpha*tmp_dt)<<std::endl;
    */
    
    this->Vgp = this->Ig/this->cavityProps.C*
                        (ii*this->cavityProps.wRF/A/B*ii*this->cavityProps.wRF*std::exp(ii*this->cavityProps.wRF*tmp_dt)+
                        (this->cavityProps.alpha+ii*this->cavityProps.wCav)/(-2.0*ii*this->cavityProps.wCav*A)*(-(ii*this->cavityProps.wCav+this->cavityProps.alpha)*std::exp(-(ii*this->cavityProps.wCav+this->cavityProps.alpha)*tmp_dt))+
                        (this->cavityProps.alpha-ii*this->cavityProps.wCav)/(-2.0*ii*this->cavityProps.wCav*B)*(-(ii*this->cavityProps.wCav-this->cavityProps.alpha)*std::exp((ii*this->cavityProps.wCav-this->cavityProps.alpha)*tmp_dt)))+
                        this->Vgm1*std::exp(-this->cavityProps.alpha*tmp_dt)*((-2.0*this->cavityProps.alpha)*std::cos(this->cavityProps.wCav*tmp_dt)+(this->cavityProps.alpha*this->cavityProps.alpha-this->cavityProps.wCav)*std::sin(this->cavityProps.wCav*tmp_dt))-
                        this->Ugm1/this->cavityProps.C/this->cavityProps.L*std::exp(-this->cavityProps.alpha*tmp_dt)*(-this->cavityProps.alpha/this->cavityProps.wCav*std::sin(this->cavityProps.wCav*tmp_dt)+std::cos(this->cavityProps.wCav*tmp_dt));
    /*
    std::cout<<"Vg after update: "<<this->Vg<<std::endl;
    std::cout<<"Vgp after update: "<<this->Vgp<<std::endl;
    std::cout<<"Ugm1 before:"<<this->Ugm1<<std::endl;
    std::cout<<"1st in Ug: "<< this->cavityProps.L*this->Ig*std::exp(ii*this->cavityProps.wRF*tmp_dt)<<std::endl;
    std::cout<<"2nd in Ug: "<< this->cavityProps.L/this->cavityProps.R*this->Vg<<std::endl;
    std::cout<<"3rd in Ug: "<< this->cavityProps.L*this->cavityProps.C*this->Vgp<<std::endl;
    */
    this->Ug = this->cavityProps.L*(this->Ig*std::exp(ii*this->cavityProps.wRF*tmp_dt)-1/this->cavityProps.R*this->Vg-this->cavityProps.C*this->Vgp);
    //std::cout<<"Ugm1 after:"<<this->Ug<<std::endl;
    this->Vgm1 = this->Vg;
    this->Ugm1 = this->Ug;
}

void Cavity::updateVb(const int& ifeed){
    //std::cout<<"Updating Vb"<<std::endl;
    //std::cout<<"Vb before update: "<<this->Vb<<std::endl;
    double tmp_dt = this->feed_time[ifeed]-this->tlast_Vb;
    this->Vb = this->Vb*std::exp((ii*this->cavityProps.wCav-this->cavityProps.alpha)*tmp_dt);
    //std::cout<<"Vb after update: "<<this->Vb<<std::endl;
}

void Cavity::updateIg(const int& ifeed){
    // calculate the generator current at the current feedback time point based on the feedback strength.
    std::complex<double> dV = this->vRef[ifeed]-this->Vg-this->Vb;
    /*
    std::cout<<"dphi :"<<this->cavityProps.wRF*(this->feed_time[ifeed]-this->tlast_feed)/2.0/M_PI<<std::endl;
    std::cout<<"rotate: "<<std::exp(ii*this->cavityProps.wRF*(this->feed_time[ifeed]-this->tlast_feed))<<std::endl;
    std::cout<<"Gain:"<<this->cavityProps.gp<<std::endl;
    std::cout<<"Ig before update: "<<this->Ig<<std::endl;    
    std::cout<<"Ig test:"<<this->Ig*std::complex<double>(0,1)<<std::endl;
    std::cout<<"Ig test2: "<<this->Ig*std::complex<double>(4.79705e-13,1)<<std::endl;    
    std::cout<<"Ig no gain:"<<this->Ig*std::exp(ii*this->cavityProps.wRF*(this->feed_time[ifeed]-this->tlast_feed))<<std::endl;
    std::cout<<"Ig caviProps:"<<this->cavityProps.Ig<<std::endl;
    std::cout<<"Ig this: "<<this->Ig<<std::endl;
    */
    this->Ig = this->Ig*std::exp(ii*this->cavityProps.wRF*(this->feed_time[ifeed]-this->tlast_feed))+
                this->cavityProps.gp*dV/this->cavityProps.R;
    //std::cout<<"Ig after update: "<<this->Ig<<std::endl;
}

void Cavity::updateLastTimes(const int& ifeed){
    this->tlast_feed = this->feed_time[ifeed];
    this->tlast_Vb = this->feed_time[ifeed];
}

void Cavity::feedback(const int& ifeed){
    // update the generator current and Vg Vgm Ug Ugm and Vb. 
    //std::cout<<"Feedback at time: "<<this->feed_time[ifeed]<<std::endl;
    //std::cout<<"Last feedback time: "<<this->tlast_feed<<std::endl;
    // These four  must be used together, other wise the time will be wrong. 
    this->updateVg(ifeed);
    this->updateVb(ifeed);
    this->updateIg(ifeed);
    this->updateLastTimes(ifeed);
}

void Cavity::kickPar(Beam& beam){
    // kick the particles 
    // This part can be parallelized since every particle is independent.
    //std::cout<<"Kicking particles"<<std::endl;
    double kick = ((this->vGen[0]+this->vBeam[0]-beam.vAdd[0]/2.0)).real();
    
    /*
    std::cout<<"Vgen: "<<this->vGen[0]<<std::endl;
    std::cout<<"Vbeam: "<<this->vBeam[0]<<std::endl;
    std::cout<<"Vadd: "<<beam.vAdd[0]<<std::endl;        
    std::cout<<"Kick: "<<kick<<std::endl;
    std::cout<<"Total kick:"<<kick-beam.Vrad<<std::endl;
    */
   #pragma omp parallel for
    for(int i = 0; i < beam.nPar*beam.nBunch; ++i){
        kick = ((this->vGen[i]+this->vBeam[i]-beam.vAdd[i]/2.0)).real();
        beam.gamma[i] += kick/beam.E0;
    }
}
void Cavity::findIstartAndIend(const int& ifeed, int& istart, int& iend){
    // Assuming the time points are sorted.
    //std::cout<<"Finding istart and iend"<<std::endl;
    istart = iend;
    while(this->time[iend] < this->feed_time[ifeed] && iend < this->time.size()){
        //std::cout<<"Current time: "<<iend<<", "<<this->time[iend]<<std::endl;
        //std::cout<<"Feed time: "<<this->feed_time[ifeed]<<std::endl;
        iend++;
    }
    //std::cout<<"istart :"<<istart<<std::endl;
    //std::cout<<"iend :"<<iend<<std::endl;
}
void Cavity::fullMap(Beam& beam){
    // First we need to analyze the beam time structure and find out how to split the job since we need to update all time points of particles and
    // update the generator voltage at each feedback time point.
    // We'd like to split the job based on the time point of feedback since this is in general less than the number of arrival time point of particles.
    
    this->updateVadd(beam);
    this->time = beam.time;
    int istart = 0;
    int iend = 0;
    for(int i=0;i<this->n_feed;++i){
        //std::cout<<"Feedback step: "<<i<<std::endl;
        this->findIstartAndIend(i,istart,iend);
        this->calVgen(beam,istart,iend,i);
        this->calVbeam(beam,istart,iend,i);
        this->feedback(i);
        if(this->iTurn >= beam.dynamicOn){
            this->kickPar(beam);
        }
        istart = iend;
    }
    beam.sortBasedOnTime();
    /*
    std::cout<<"generator voltage of first particle: "<<this->vGen[0]<<std::endl;
    std::cout<<"beam voltage of first particle: "<<this->vBeam[0]<<std::endl;
    std::cout<<"Vtot of first particle: "<<this->vGen[0]+this->vBeam[0]<<std::endl;
    std::cout<<"Vtot of first particle: "<<this->vGen[0]+this->vBeam[0]-beam.vAdd[0]/2.0<<std::endl;
    std::cout<<"Arriving time of first particle"<<beam.time[0]-beam.bucket_centers[0]<<std::endl;
    std::cout<<"Arriving phase of first particle"<<(beam.time[0]-beam.bucket_centers[0])*this->cavityProps.wRF/M_PI*180<<std::endl;
    std::cout<<"Gamma of first particle:"<<beam.gamma[0]-beam.Gamma0<<std::endl;
    */
    this->iTurn++;
}