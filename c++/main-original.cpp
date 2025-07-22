
#include<iostream>
#include<math.h>
#include<stdlib.h>
#include<fstream>
#include<sstream>
#include<vector>
#include<ctime>
using namespace std;

int stop=0; //put =1 if wanting to stop
int speak=0; // out =1 if wanting extra info througout sim

// max Monte Carlo Sweeps
const int Tmax=10000;
int t;

// Grid size
const int L=60;

///////////////////////////
// HAMILTONIAN:
// H= \sum_{m,n} J(sigma(n),sigma(n)) + Lambda \sum_i^N (a_i-A)^2
// J = alpha if sigma(n) != sigma(m)
// J = 0 if sigma(n) = sigma(m)
// the sum is over {n,m} is over n lattice sites and m neighbours of n
// a_i is instantaneous area of cell i
///////////////////////////////

///////////////////////////////
// PARAMETERS OF HAMILTONIAN //
const int A=16;
double lambda=1;
double alpha_soft=1;
double P=0;
double alpha_stiff=1; //to be set by user
double alpha_stiffest=1;
double alpha_stiffer=1; //to be set by user
///////////////////////////////

//Heterotypic Boundary Length//

double sum_AA;
double sum_AB;
double sum_BA;
double sum_AC;
double sum_CA;
double sum_BB;
double sum_BC;
double sum_CB;
double sum_CC;
//
double sum_AD;
double sum_BD;
double sum_CD;
double sum_DD;
double sum_DA;
double sum_DB;
double sum_DC;


double agg_AA;
double agg_AB_BA;
double agg_AC_CA;
double agg_BB;
double agg_BC_CB;
double agg_CC;
//
double agg_DD;
double agg_AD_DA;
double agg_BD_DB;
double agg_CD_DC;

double r_AA;
double r_BB;
double r_CC;
double r_AB;
double r_AC;
double r_BC;
//
double r_DD;
double r_AD;
double r_BD;
double r_CD;

//Gyration tensor summations//

int num_x;
int num_y;
int num_z;
double sum_Rgyr_x;
double sum_Rgyr_y;
double sum_Rgyr_z;

double sum_aspher_x;
double sum_aspher_y;
double sum_aspher_z;


// 2-dimensional array of cell types (or spins)
int Cell[L][L];

// to keep track of PBC
int boxflag[L][L][2];

// TOTAL N OF CELLS //
const int Ncells=int(L*L*1.0/A);

//com
double COM[Ncells][2];
double COMold[Ncells][2];
// unwrapped
int COM_unwrapped[Ncells][2];
int pbcflag[Ncells][2];

//Useful arrays
double Lambda[Ncells];
double Alpha[Ncells];
int Area[Ncells];

double Type[Ncells]; // this distinguishes type A=sender from type B=receiver, C(=B'), D=(A')

double HeterotropicBoundary[Ncells][16];

////TENSOR
double GyrTens[Ncells][3];
double lambda1[Ncells];
double lambda2[Ncells];
double Rgyr[Ncells];
double aspher[Ncells];
double AreaNEW[Ncells];

//GENE Switching
double kon=0;
double koff=0;
double kon2=0;
double koff2=0;

//Energies
double Eold;
double Enew;

//Functions
double J(int a,int b,double alpha);
double ComputeEnergy(int x, int y);
void ComputeGyrationTensor();

//START PROGRAM
int main(int argc, char* argv[]){
    
    cout<< "## ATTENTION! HAVE YOU CREATED THE FOLDER IN WHICH I SHOULD WRITE??? ##" <<endl;
    cout<< "Type: 1. name output 2. cell adherence (stiff cells) 3. Fraction of B=receiver cells 4. Switch Rate ON B->C; 5. Switch Rate OFF C->B; 6. cell adherence (stiffer cells) 7. cell adherence (stiffest cells) 8.Switch rate on A->D 9. Switch rate off D->A " <<endl;
    srand(time(NULL));
    
    // CPU CLOCK //
    clock_t tcpu;
    tcpu = clock();
    
    // STIFFNESS PARAMETER //
    alpha_stiff=atof(argv[2]);
    alpha_stiffer=atof(argv[6]);
    alpha_stiffest=atof(argv[7]);
    
    // FRACTION OF B CELLS //
    double fB=atof(argv[3]);
    
    //RATES//
    kon=atof(argv[4]);
    koff=atof(argv[5]);
    kon2=atof(argv[8]);
    koff2=atof(argv[9]);
    // CREATE NAMES OF FILE AND OUTPUT //
    stringstream writeTF0;
    writeTF0<<"Alpha"<<alpha_stiff<<"/"<<argv[1]<<"_t"<<0;
    ofstream write(writeTF0.str().c_str());
    stringstream writecomF;
    writecomF <<"Alpha"<<alpha_stiff<<"/"<< "com_"<<argv[1];
    ofstream writeCOM(writecomF.str().c_str());
    
    // SET PENALTY FOR AREA CHANGE //
    for(int s=0;s<Ncells;s++){
        Lambda[s]=1;
    }
    
    // AT THE BEGINNING ALL CELLS ARE SOFT (both A and B cells) //
    for(int s=0;s<Ncells;s++)
    Alpha[s]=alpha_stiff;
    
    for(int s=0;s<Ncells;s++){
        if(rand()*1.0/RAND_MAX<fB){
            Type[s]=1;
            Alpha[s]=alpha_soft;
            
        }// B cell = receiver -> can become C (Sticky)
        else Type[s]=0;
        // A cell = sender (remains A, non sticky)
    }
    
    cout << "############## time " << 0 << " #################" << endl;
    cout << "writing on " << writeTF0.str().c_str()<<endl;
    
    /////////////////////////////////////////////////////////////////////////////////////
    //initialise cells in a "plaid" pattern
    //////////////////////////////////////////////////////////////////////////////
    int maxa=4;
    int type=0;
    for(int c1=0;c1<L/maxa;c1++){
        //row
        for(int c2=0;c2<L/maxa;c2++){
            //column
            for(int i=0;i<maxa;i++){
                for(int j=0;j<maxa;j++){
                    boxflag[i+c1*maxa][j+c2*maxa][0]=0;
                    boxflag[i+c1*maxa][j+c2*maxa][1]=0;
                    Cell[i+c1*maxa][j+c2*maxa]=type;
                    write <<  i+c1*maxa+0.5 << " " << j+c2*maxa+0.5<< " "<<Cell[i+c1*maxa][j+c2*maxa]<< " " << Alpha[Cell[i+c1*maxa][j+c2*maxa]] << " " << Type[Cell[i+c1*maxa][j+c2*maxa]] << endl;
                }
            }
            type++;
        }
    }
    cout << "write com" <<endl;
    writeCOM<<"#t" <<t <<endl;
    for(int s=0;s<Ncells;s++){
        COM[s][0]=0;
        COM[s][1]=0;
        COM_unwrapped[s][0]=0;
        COM_unwrapped[s][1]=0;
        pbcflag[s][0]=0;
        pbcflag[s][1]=0;
        
        for(int l1=0;l1<L;l1++){
            for(int l2=0;l2<L;l2++){
                if(Cell[l1][l2]==s){
                    COM[s][0]+=l1*1.0;
                    COM[s][1]+=l2*1.0;
                    COM_unwrapped[s][0]+=l1*1.0;
                    COM_unwrapped[s][1]+=l2*1.0;
                    //cout <<"cell " << s << " @ " << l1 << " " << l2 << " and summing " <<l1*1.0+boxflag[l1][l2][0]*L << " " <<l2*1.0+boxflag[l1][l2][1]*L << " with flags " << boxflag[l1][l2][0] << " & " << boxflag[l1][l2][1] << " w area " << currentArea[s] <<endl;
                }
            }
        }
        COM[s][0]=COM[s][0]/(maxa*maxa);
        COM[s][1]=COM[s][1]/(maxa*maxa);
        COM_unwrapped[s][0]=COM_unwrapped[s][0]/(maxa*maxa);
        COM_unwrapped[s][1]=COM_unwrapped[s][1]/(maxa*maxa);
        //writeCOM << s << " " << COM[s][0]+0.5 << " " << COM[s][1]+0.5 << " " << COM_unwrapped[s][0]+0.5 << " " << COM_unwrapped[s][1]+0.5<< endl;
        writeCOM << s << " " << COM_unwrapped[s][0]+0.5 << " " << COM_unwrapped[s][1]+0.5<< endl;
    }
    writeCOM<<endl;
    writeCOM<<endl;
    
    // COMPUTE INITIAL AREAS //
    for(int l1=0;l1<L;l1++){
        for(int l2=0;l2<L;l2++){
            Area[Cell[l1][l2]]++;
        }
    }
    if(speak)cout << "Initialised" <<endl;
    
    ////COUNT CELLS
    int nB=0; int nA=0; int nC=0; int nD=0;
    for(int s=0;s<Ncells;s++){
        if(Type[s]==0)nA++;
        if(Type[s]==1)nB++;
        if(Type[s]==2)nC++;
        if(Type[s]==4)nD++;
    }
    cout << "#CELLS: A=" << nA << ";  B="<<nB << "; C="<<nC<<"; D=" <<endl;
    
    
    ///////////////////////
    //END initialisation //
    ///////////////////////
    
    //////////////////////////
    // START DYNAMICS       //
    //////////////////////////
    cout << "Start time"<<endl;
    //START monte carlo sweep
    for(t=1;t<Tmax;t++){
        
        if(speak)cout << "Attempts" <<endl;
        
        //N attempts per monte carlo sweep//
        int Nattempts=L*L;
        
        int n=0;
        while(1==1){
            
            // STOP IF MADE ENOUGH SUCCESSFUL ATTEMPTS
            if(n==Nattempts)break;
            
            // ENERGIES //
            Eold=0;
            Enew=0;
            
            //////////////////////
            // Pick random site //
            //////////////////////
            int m=int(rand()*1.0/RAND_MAX*L*L);
            int mx=m%L; //x coordinate of m site
            int my=int(m*1.0/L); // y coordingate of m site
            if(speak)cout << "considering site " <<mx << " " <<my <<endl;
            
            //find neighbour of mx,my in 2D
            int nx,ny;
            int p=rand()%8;
            if(p==0){nx=mx+1;ny=my;} // 0 -> dx;
            if(p==1){nx=mx-1;ny=my;} // 1 -> sx;
            if(p==2){nx=mx;ny=my+1;} // 2 -> up;
            if(p==3){nx=mx;ny=my-1;} // 3 -> dwn;
            if(p==4){nx=mx+1;ny=my+1;} // 4 -> dx up;
            if(p==5){nx=mx+1;ny=my-1;} // 5 -> dx dwn;
            if(p==6){nx=mx-1;ny=my+1;} // 6 -> sx up;
            if(p==7){nx=mx-1;ny=my-1;} // 7 -> sx dwn;
            //reflecting bc
            //if(nx==L)nx=L-2;if(ny==L)ny=L-2;if(nx==-1)nx=2;if(ny==-1)ny=2;
            //periodic
            if(nx==L)nx=0;if(ny==L)ny=0;if(nx==-1)nx=L-1;if(ny==-1)ny=L-1;
            if(speak)cout << "neighbour " <<nx << " " <<ny <<endl;
            
            
            
            
            
            ////////////////////////////////
            ///CALC HETEROTROPIC BOUNDARY
            /////////////////////////////
            int px[8]={+1,-1,0,0,+1,+1,-1,-1};
            int py[8]={0,0,+1,-1,+1,-1,+1,-1};
            
            for(int pp=0;pp<8;pp++){
                
                int nnx=mx+px[pp];
                int nny=my+py[pp];
                //PBC
                if(nnx==L)nnx=0;if(nny==L)nny=0;if(nnx==-1)nnx=L-1;if(nny==-1)nny=L-1;
                //
                if(Cell[nnx][nny]!=Cell[mx][my]){
                    
                    if(Type[Cell[mx][my]]==0 && Type[Cell[nnx][nny]]==0) //AA
                        HeterotropicBoundary[Cell[mx][my]][0]++;
                    
                    if(Type[Cell[mx][my]]==0 && Type[Cell[nnx][nny]]==1) //AB
                        HeterotropicBoundary[Cell[mx][my]][1]++;
                    
                    if(Type[Cell[mx][my]]==0 && Type[Cell[nnx][nny]]==2) //AC
                        HeterotropicBoundary[Cell[mx][my]][2]++;
                    
                    if(Type[Cell[mx][my]]==1 && Type[Cell[nnx][nny]]==0) //BA
                        HeterotropicBoundary[Cell[mx][my]][3]++;
                    
                    if(Type[Cell[mx][my]]==1 && Type[Cell[nnx][nny]]==1) //BB
                        HeterotropicBoundary[Cell[mx][my]][4]++;
                    
                    if(Type[Cell[mx][my]]==1 && Type[Cell[nnx][nny]]==2) //BC
                        HeterotropicBoundary[Cell[mx][my]][5]++;
                    
                    if(Type[Cell[mx][my]]==2 && Type[Cell[nnx][nny]]==0) //CA
                        HeterotropicBoundary[Cell[mx][my]][6]++;
                    
                    if(Type[Cell[mx][my]]==2 && Type[Cell[nnx][nny]]==1) //CB
                        HeterotropicBoundary[Cell[mx][my]][7]++;
                    
                    if(Type[Cell[mx][my]]==2 && Type[Cell[nnx][nny]]==2) //CC
                        HeterotropicBoundary[Cell[mx][my]][8]++;
                    
                    //
                    if(Type[Cell[mx][my]]==0 && Type[Cell[nnx][nny]]==3) //AD
                        HeterotropicBoundary[Cell[mx][my]][9]++;
                    if(Type[Cell[mx][my]]==1 && Type[Cell[nnx][nny]]==3) //BD
                        HeterotropicBoundary[Cell[mx][my]][10]++;
                    if(Type[Cell[mx][my]]==2 && Type[Cell[nnx][nny]]==3) //CD
                        HeterotropicBoundary[Cell[mx][my]][11]++;
                    if(Type[Cell[mx][my]]==3 && Type[Cell[nnx][nny]]==3) //DD
                        HeterotropicBoundary[Cell[mx][my]][12]++;
                    if(Type[Cell[mx][my]]==3 && Type[Cell[nnx][nny]]==0) //DA
                        HeterotropicBoundary[Cell[mx][my]][13]++;
                    if(Type[Cell[mx][my]]==3 && Type[Cell[nnx][nny]]==1) //DB
                        HeterotropicBoundary[Cell[mx][my]][14]++;
                    if(Type[Cell[mx][my]]==3 && Type[Cell[nnx][nny]]==2) //DC
                        HeterotropicBoundary[Cell[mx][my]][15]++;
                    
                    
                }
                
            }
            //can we just write out AA, AB, AC, BB, BC, CC by symmetry? -- no: because our centre cell is at mx, my. Heterotropich boundary conditions will only be added to if we consider cells of A, B, C as centre cell.
            cout << "site " << mx << " " << my << " Cell Type " << Type[Cell[mx][my]]<< " -> " "  A-A: " <<HeterotropicBoundary[Cell[mx][my]][0]<< "  A-B: " <<HeterotropicBoundary[Cell[mx][my]][1]<< "  AC: " <<HeterotropicBoundary[Cell[mx][my]][2]<< "  AC: " <<HeterotropicBoundary[Cell[mx][my]][2]<< "  BA: " <<HeterotropicBoundary[Cell[mx][my]][3]<< "  BB: " <<HeterotropicBoundary[Cell[mx][my]][4]<< "  BC: " <<HeterotropicBoundary[Cell[mx][my]][5]<< "  CA: " <<HeterotropicBoundary[Cell[mx][my]][6]<< "  CB: " <<HeterotropicBoundary[Cell[mx][my]][7]<< "  CC: " <<HeterotropicBoundary[Cell[mx][my]][8]<<endl;
            //cin.get();
            
            
            
            
            
            
            //////////////////////////////////////////
            //IF THE SPINS ARE DIFFERENT
            //GO AHEAD WITH ATTEMPT
            //OTHERWISE TRY ANOTHER PAIR OF SPINS
            ////////////////////////////////////////
            if(Cell[mx][my]!=Cell[nx][ny]){
                
                
                ////////////////////
                // GENE SWITCHING STEP -- LOOK @ THE NEIGHBOUR
                // IF A IS NEAR B -> B can become C and C is adherent
                // IF a cell is C -> send it back to B with certain probability
                //////////////////////
                double pswitch1=rand()*1.0/RAND_MAX; //random number [0,1)
                double pswitch2=rand()*1.0/RAND_MAX; //random number [0,1)
                //ON
                //IF the site mx,my belongs to a B cell and it is near to a A cell -> change B to C with rate kon
                if(Type[Cell[mx][my]]==1 && Type[Cell[nx][ny]]==0){
                    if(pswitch1<kon){
                        Type[Cell[mx][my]]=2; // this is a C cell
                        Alpha[Cell[mx][my]]=alpha_stiffest;
                    }
                }
                //OFF
                //if it is a C cell -> return to a B cell
                else if(Type[Cell[mx][my]]==2){
                    if(pswitch2<koff){
                        Type[Cell[mx][my]]=1; // C goes to B
                        Alpha[Cell[mx][my]] = alpha_soft; // Alpha becomes soft
                    }
                }
                
                double pswitch3=rand()*1.0/RAND_MAX; //random number [0,1)
                double pswitch4=rand()*1.0/RAND_MAX; //random number [0,1)
                //ON
                //IF the site mx,my belongs to a A cell and it is near to a C cell -> change A to D with rate kon
                if(Type[Cell[mx][my]]==0 && Type[Cell[nx][ny]]==2){
                    if(pswitch3<kon2){
                        Type[Cell[mx][my]]=3; // this is a C cell
                        Alpha[Cell[mx][my]]=alpha_stiffer;
                    }
                }
                
                //OFF
                //if it is a D cell -> return to a A cell
                else if(Type[Cell[mx][my]]==3){
                    if(pswitch4<koff2){
                        Type[Cell[mx][my]]=0; // D goes to A
                        Alpha[Cell[mx][my]] = alpha_stiff; // Alpha becomes softer
                    }
                }
                
             
                
                ///////////////////////
                // END GENETIC SWITCH
                ///////////////////////
                
             
                
                
                //////////////////
                //MAKE CELLS MOVE
                //////////////////
                
                ////////////////
                // TO MAKE COMPUTATION EFFICIENT, WE CALCULATE ONLY LOCAL CHANGE OF ENERGY
                // THE REST OF THE SPINS REMAIN THE SAME SO THE ENERGY IS UNAFFECTED
                /////////////
                //compute old energy at site m
                double Eold=ComputeEnergy(mx,my);
                if(speak)cout << "Eold " <<Eold <<endl;
                if(stop)cin.get();
                if(isnan(Eold)){cout << "ALERT!! NaN!" <<endl; cin.get();}
                
                //Record value of old spin (old cell id)
                int oldspin=Cell[mx][my];
                //Change areas
                Area[Cell[mx][my]]--;
                Area[Cell[nx][ny]]++;
                //Change spin
                Cell[mx][my]=Cell[nx][ny]; //Spin @ mx,my is attempted to be changed into the one @ nx,ny
                //compute new energy at site m
                double Enew=ComputeEnergy(mx,my);
                if(speak)cout << "Enew " <<Enew <<endl;
                if(stop)cin.get();
                if(isnan(Enew)){cout << "ALERT!! NaN!" <<endl; cin.get();}
                
                //Metropolis test via change in energy Enew-Eold
                double DE=Enew-Eold;
                //Metropolis Criterion (DE is in units of k_B T)
                double MetProb=exp(-DE);
                if(speak)cout << DE << " " << MetProb <<endl;
                if(rand()*1.0/RAND_MAX<fmin(MetProb,1)){
                    //ACCEPTED:MAINTAIN CHANGE OF SPIN AND AREAS
                    if(speak)cout << "Accepted!! " << endl;
                }
                else{
                    if(speak)cout<<"rejected!"<<endl;
                    //REJECTED:RETURN TO OLD VALUES
                    Cell[mx][my]=oldspin;
                    Area[Cell[mx][my]]++;
                    Area[Cell[nx][ny]]--;
                }
                
                if(stop)cin.get();
                
                //COUNT AS ATTEMPT AND TRY ANOTHER
                n++;
            } //close if spin(m) != spin(n)
        }//end loop over attempts
        
        
        
        for(int k=0; k<Ncells; k++){
            sum_AA+=HeterotropicBoundary[k][0];
            sum_AB+=HeterotropicBoundary[k][1];
            sum_BA+=HeterotropicBoundary[k][3];
            sum_AC+=HeterotropicBoundary[k][2];
            sum_CA+=HeterotropicBoundary[k][6];
            sum_BB+=HeterotropicBoundary[k][4];
            sum_BC+=HeterotropicBoundary[k][5];
            sum_CB+=HeterotropicBoundary[k][7];
            sum_CC+=HeterotropicBoundary[k][8];
            
            sum_AD+=HeterotropicBoundary[k][9];
            sum_BD+=HeterotropicBoundary[k][10];
            sum_CD+=HeterotropicBoundary[k][11];
            sum_DD+=HeterotropicBoundary[k][12];
            sum_DA+=HeterotropicBoundary[k][13];
            sum_DB+=HeterotropicBoundary[k][14];
            sum_DC+=HeterotropicBoundary[k][15];
            
            
        }
        
        
        
        //cout << "there" << endl;
        agg_AA = sum_AA;
        agg_AB_BA= sum_AB + sum_BA;
        agg_AC_CA = sum_AC + sum_CA;
        agg_BB = sum_BB;
        agg_BC_CB = sum_BC + sum_CB;
        agg_CC = sum_CC;
        
        agg_DD = sum_DD;
        agg_AD_DA = sum_AD + sum_DA;
        agg_BD_DB = sum_BD + sum_DB;
        agg_CD_DC = sum_CD + sum_DC;
        
        
        
        
        //cout<< agg_AA << " " << agg_AB_BA<< " "<< agg_AC_CA << " " << agg_BB << " " << agg_BC_CB << " " << agg_CC << endl;
        
        r_AA = agg_AA/(agg_AA + agg_AB_BA + agg_AC_CA + agg_BB + agg_BC_CB + agg_CC + agg_DD + agg_AD_DA + agg_BD_DB + agg_CD_DC);
        r_BB = agg_BB/(agg_AA + agg_AB_BA + agg_AC_CA + agg_BB + agg_BC_CB + agg_CC + agg_DD + agg_AD_DA + agg_BD_DB + agg_CD_DC);
        r_CC = agg_CC/(agg_AA + agg_AB_BA + agg_AC_CA + agg_BB + agg_BC_CB + agg_CC + agg_DD + agg_AD_DA + agg_BD_DB + agg_CD_DC);
        r_AB = agg_AB_BA/(agg_AA + agg_AB_BA + agg_AC_CA + agg_BB + agg_BC_CB + agg_CC + agg_DD + agg_AD_DA + agg_BD_DB + agg_CD_DC);
        r_AC = agg_AC_CA/(agg_AA + agg_AB_BA + agg_AC_CA + agg_BB + agg_BC_CB + agg_CC + agg_DD + agg_AD_DA + agg_BD_DB + agg_CD_DC);
        r_BC = agg_BC_CB/(agg_AA + agg_AB_BA + agg_AC_CA + agg_BB + agg_BC_CB + agg_CC + agg_DD + agg_AD_DA + agg_BD_DB + agg_CD_DC);
        
        r_DD = agg_DD/(agg_AA + agg_AB_BA + agg_AC_CA + agg_BB + agg_BC_CB + agg_CC + agg_DD + agg_AD_DA + agg_BD_DB + agg_CD_DC);
        r_AD = agg_AD_DA/(agg_AA + agg_AB_BA + agg_AC_CA + agg_BB + agg_BC_CB + agg_CC + agg_DD + agg_AD_DA + agg_BD_DB + agg_CD_DC);
        r_BD = agg_BD_DB/(agg_AA + agg_AB_BA + agg_AC_CA + agg_BB + agg_BC_CB + agg_CC + agg_DD + agg_AD_DA + agg_BD_DB + agg_CD_DC);
        r_CD = agg_CD_DC/(agg_AA + agg_AB_BA + agg_AC_CA + agg_BB + agg_BC_CB + agg_CC + agg_DD + agg_AD_DA + agg_BD_DB + agg_CD_DC);
        
        //cout << r_AA << " " << r_BB << " " << r_CC <<" " << r_DD << " " << r_AB << " " << r_AC << " "<< r_BC <<" " <<r_AD << " " << r_BD << " " << r_CD << endl;
        
        
        //reset the HBL array for each loop
        memset(HeterotropicBoundary, 0, sizeof(HeterotropicBoundary));
        
        
        
        sum_AA=0;
        sum_AB=0;
        sum_BA=0;
        sum_AC=0;
        sum_CA=0;
        sum_BB=0;
        sum_BC=0;
        sum_CB=0;
        sum_CC=0;
        
        agg_AA=0;
        agg_AB_BA=0;
        agg_AC_CA=0;
        agg_BB=0;
        agg_BC_CB=0;
        agg_CC=0;
        
        r_AA=0;
        r_BB=0;
        r_CC=0;
        r_AB=0;
        r_AC=0;
        r_BC=0;
        
        
        
        
        
        
        sum_AD=0;
        sum_BD=0;
        sum_CD=0;
        sum_DD=0;
        sum_DA=0;
        sum_DB=0;
        sum_DC=0;
        
        
        agg_DD=0;
        agg_AD_DA=0;
        agg_BD_DB=0;
        agg_CD_DC=0;
        
        r_DD=0;
        r_AD=0;
        r_BD=0;
        r_CD=0;
        
        
        
        
        
        
        /*
         ////////
         ComputeGyrationTensor();
         /////////
         */
        
        //EVERY 10 MC sweeps print something
        if(t%10==0){
            clock_t elapsedtime = clock() - tcpu;
            //////////WRITE CONF AT TIME T //////////////////////////////
            cout << "############## time " << t << " #################" << endl;
            stringstream writeTF;
            writeTF<<"Alpha"<<alpha_stiff<<"/"<<argv[1]<<"_t"<<t;
            ofstream writeT(writeTF.str().c_str());
           cout <<"writing on " <<writeTF.str().c_str()<<endl;
            for(int l1=0;l1<L;l1++)for(int l2=0;l2<L;l2++)writeT<< l1+0.5 <<" " <<l2+0.5 << " " << Cell[l1][l2]<< " " << Alpha[Cell[l1][l2]]<< " " << Type[Cell[l1][l2]] << endl;
            cout << "CPU Time = " << ((float)elapsedtime)/CLOCKS_PER_SEC << " seconds" << endl;
            /////////////////////////////////////////////////////////////////////////////////////
            
            
            ////FIND IF ACROSS BORDERS
            int across_x[Ncells];
            int across_y[Ncells];
            for(int s=0;s<Ncells;s++){
                across_x[s]=0;
                across_y[s]=0;
                int max_x=int(L/2.);
                int min_x=int(L/2.);
                int max_y=int(L/2.);
                int min_y=int(L/2.);
                for(int l1=0;l1<L;l1++){
                    for(int l2=0;l2<L;l2++){
                        if(Cell[l1][l2]==s){
                            //find extrema of cell in the x direction
                            if(l1<min_x)min_x=l1;
                            if(l1>max_x)max_x=l1;
                            //find extrema of cell in the y direction
                            if(l2<min_y)min_y=l2;
                            if(l2>max_y)max_y=l2;
                            //check that min and max are on the same side
                            if(min_x<L/2. && max_x>L/2. && max_x-min_x>L/2. && COMold[s][0]>L/2.) across_x[s]=1;
                            if(min_x<L/2. && max_x>L/2. && max_x-min_x>L/2. && COMold[s][0]<L/2.) across_x[s]=-1;
                            if(min_y<L/2. && max_y>L/2. && max_y-min_y>L/2. && COMold[s][1]>L/2.) across_y[s]=1;
                            if(min_y<L/2. && max_y>L/2. && max_y-min_y>L/2. && COMold[s][1]<L/2.) across_y[s]=-1;
                            //
                            
                        }
                    }
                }
            }
            
            ////COUNT CELLS
            int nB=0; int nA=0; int nC=0;int nD=0;
            for(int s=0;s<Ncells;s++){
                if(Type[s]==0)nA++;
                if(Type[s]==1)nB++;
                if(Type[s]==2)nC++;
                if(Type[s]==3)nD++;
            }
            cout << "#CELLS: A=" << nA << ";  B="<<nB << "; C="<<nC<< "; D="<<nD <<endl;
            
            ////COMPUTE COM OF CELLS AND WRITE
           // cout << "write com" <<endl;
            writeCOM<<"#t" <<t <<endl;
            int carea[Ncells];
            for(int s=0;s<Ncells;s++){
                //
                carea[s]=0;
                COM[s][0]=0;
                COM[s][1]=0;
                for(int l1=0;l1<L;l1++){
                    for(int l2=0;l2<L;l2++){
                        if(Cell[l1][l2]==s){
                            //
                            int addy=0;
                            int addx=0;
                            //
                            if(across_x[s]==1 && l1<L/2.)addx=L;
                            if(across_x[s]==-1 && l1>L/2.)addx=-L;
                            if(across_y[s]==1 && l2<L/2.)addy=L;
                            if(across_y[s]==-1 && l2>L/2.)addy=-L;
                            //
                            //cout << " cell " << s << " @ " << l1 << " " << l2 << " " << across_x[s]<<" " << across_y[s]<<  " adding " << l1+addx<< " " << l2+addy <<endl;cin.get();
                            //
                            COM[s][0]+=l1*1.0+addx;
                            COM[s][1]+=l2*1.0+addy;
                            //
                            carea[s]++;
                            //cout <<"cell " << s << " @ " << l1 << " " << l2 << " and summing " <<l1*1.0+boxflag[l1][l2][0]*L << " " <<l2*1.0+boxflag[l1][l2][1]*L << " with flags " << boxflag[l1][l2][0] << " & " << boxflag[l1][l2][1] << " w area " << currentArea[s] <<endl;
                        }
                    }
                }
                COM[s][0]=COM[s][0]/carea[s];
                COM[s][1]=COM[s][1]/carea[s];
                
                //    writeCOM << s << " " <<COM[s][0] << " " << COM[s][1]  << " " <<  COM[s][0]-floor(COM[s][0]*1.0/L)*L << " " << COM[s][1] - floor(COM[s][1]*1.0/L)*L <<endl;
                writeCOM << s << " " << COM[s][0]-floor(COM[s][0]*1.0/L)*L << " " << COM[s][1] - floor(COM[s][1]*1.0/L)*L <<endl;
                COMold[s][0]=COM[s][0];
                COMold[s][1]=COM[s][1];
            }
            writeCOM<<endl;
            writeCOM<<endl;
            /////////////////////////////////////////////////////////////////////////////////////
        }
        
        if(t%10000==0)cin.get();
    } //close loop over time
    /////////////////////////////////////////////////////////////////////////////////////
    
    
    return 0;
}






///////////////////
// FUNCTIONS
//////////////////
double J(int a,int b, double alpha){
    double j;
    if(a==b)j=0;
    if(a!=b)j=alpha;
    return j;
}

/////////////////////////////////////
/////////////////////////////////////
double ComputeEnergy(int x, int y){
    double E=0;
    
    //loop over nearest neighbours;
    int kx,ky;
    for(int nn=0;nn<8;nn++){
        if(nn==0){kx=x+1;ky=y;} //0 -> dx;
        if(nn==1){kx=x-1;ky=y;} //1 -> sx;
        if(nn==2){kx=x;ky=y+1;} //2 -> up;
        if(nn==3){kx=x;ky=y-1;} //3 -> dwn;
        if(nn==4){kx=x+1;ky=y+1;} //4 -> dx up;
        if(nn==5){kx=x+1;ky=y-1;} //5 -> dx dwn;
        if(nn==6){kx=x-1;ky=y+1;} //6 -> sx up;
        if(nn==7){kx=x-1;ky=y-1;} //7 -> sx dwn;
        //reflecting bc
        //if(kx==L)kx=L-2;if(ky==L)ky=L-2;if(kx==-1)kx=2;if(ky==-1)ky=2;
        //periodic bc
        if(kx==L)kx=0;if(ky==L)ky=0;if(kx==-1)kx=L-1;if(ky==-1)ky=L-1;
        
        if(speak)cout << x << " " << y << " " << Cell[x][y]<<endl;
        if(speak)cout << kx << " " << ky << " " << Cell[kx][ky]<<endl;
        if(speak)cout << "J " <<J(Cell[x][y],Cell[kx][ky],Alpha[Cell[x][y]])<<endl;
        E+=J(Cell[x][y],Cell[kx][ky],Alpha[Cell[x][y]]);
    }
    //cout << "Eold1 " << Eold<<endl;
    
    //add area constraint
    //loop over cells and sum lambda*(a-A);
    for(int s=0;s<Ncells;s++){
        int currentarea=Area[s];
        E+=Lambda[s]*pow((currentarea-A),2.0);
        if(speak)cout << "old cell "<< s << " current area " << currentarea<<endl;
    }
    //cout << "Eold2 " << Eold<<endl;
    
    //OPTIONAL FOR LATER (NOT DONE AT THE MOMENT)
    // ADD MOBILITY TERM
    /*
     //add mobility term
     double direction[Ncells][2];
     double theta[Ncells];
     double noise;
     //direction is rotational diffusive;
     double Dr=1.0;
     //
     for(int s=0;s<Ncells;s++){
     theta[s] =rand()*1.0/RAND_MAX*2*M_PI;
     //noise=rand()*1.0/RAND_MAX;
     //theta[s]=noise*2*Dr;
     direction[s][0]=cos(theta[s]);
     direction[s][1]=sin(theta[s]);
     //
     E-=P*((x*direction[s][0])+(y*direction[s][1])); // this is P times (x dot n)
     //if(speak)cout << "cell " << s << " direction " <<P*((x*direction[s][0])+(y*direction[s][1]))<<endl;
     }
     */
    
    return E;
}


/*
 void ComputeGyrationTensor(){
 
 //LOOP OVER CELLS
 for(int n=0;n<Ncells;n++){
 //cout << " NN " << n <<endl;
 GyrTens[n][0]=GyrTens[n][1]=GyrTens[n][2]=0;
 lambda1[n]=0;
 lambda2[n]=0;
 Rgyr[n]=0;
 AreaNEW[n]=0;
 //LOOP OVER LATTICE POINTS
 for(int i=0;i<L;i++){
 for(int j=0;j<L;j++){
 //LOOP OVER LATTICE POINTS
 for(int ii=0;ii<L;ii++){
 for(int jj=0;jj<L;jj++){
 if(Cell[i][j]==n && Cell[ii][jj]==n){
 double dx=(i-ii); if(dx>L/2.)dx=dx-L; else if(dx<-L/2.)dx=dx+L;
 double dy=(j-jj); if(dy>L/2.)dy=dy-L; else if(dy<-L/2.)dy=dy+L;
 GyrTens[n][0]+=dx*dx;
 GyrTens[n][1]+=dx*dy;
 GyrTens[n][2]+=dy*dy;
 AreaNEW[n]++;
 }
 }
 }
 }
 }
 for(int d=0;d<3;d++)GyrTens[n][d]=GyrTens[n][d]/AreaNEW[n];
 
 double b=-(GyrTens[n][0]+GyrTens[n][2]);
 double a=1;
 double c=GyrTens[n][0]*GyrTens[n][2]-GyrTens[n][1]*GyrTens[n][1];
 
 lambda1[n]=-b/(2.*a)+pow((b*b-4*a*c)/(4*a),0.5);
 lambda2[n]=-b/(2.*a)-pow((b*b-4*a*c)/(4*a),0.5);
 
 Rgyr[n]=sqrt(pow(lambda1[n],2.0)+pow(lambda2[n],2.0));
 aspher[n]=lambda2[n]/lambda1[n];
 
 
 
 
 
 
 
 for(int k=0; k<Ncells; k++){
 
 if(Type[k]==0){
 sum_Rgyr_x+=Rgyr[k];
 num_x+=1;
 }
 if(Type[k]==1){
 sum_Rgyr_y+=Rgyr[k];
 num_y+=1;
 }
 if(Type[k]==2){
 sum_Rgyr_z+=Rgyr[k];
 num_z+=1;
 }
 }
 
 double agg_Rgyr_x = sum_Rgyr_x/num_x;
 double agg_Rgyr_y = sum_Rgyr_y/num_y;
 double agg_Rgyr_z = sum_Rgyr_z/num_z;
 
 
 
 for(int k=0; k<Ncells; k++){
 
 if(Type[k]==0){
 sum_aspher_x+=aspher[k];
 num_x+=1;
 }
 if(Type[k]==1){
 sum_aspher_y+=aspher[k];
 num_y+=1;
 }
 if(Type[k]==2){
 sum_aspher_z+=aspher[k];
 num_z+=1;
 }
 }
 
 double agg_aspher_x = sum_aspher_x/num_x;
 double agg_aspher_y = sum_aspher_y/num_y;
 double agg_aspher_z = sum_aspher_z/num_z;
 
 
 
 for(int k=0; k<Ncells; k++){
 
 if(Type[k]==0){
 num_x+=1;
 }
 if(Type[k]==1){
 num_y+=1;
 }
 if(Type[k]==2){
 num_z+=1;
 }
 }
 
 cout<<num_x<< " "<< num_y<< " "<< num_z << endl;
 
 
 
 ////COUNT CELLS
 int mB=0; int mA=0; int mC=0;
 for(int s=0;s<Ncells;s++){
 if(Type[s]==0)mA++;
 if(Type[s]==1)mB++;
 if(Type[s]==2)mC++;
 }
 cout << "#CELLS: A=" << mA << ";  B="<<mB << "; C="<<mC<<endl;
 
 
 //how do we initiliase back to zero?
 
 if (num_x % 224 ==0){
 sum_Rgyr_x=0;
 sum_Rgyr_y=0;
 sum_Rgyr_z=0;
 num_x=0;
 num_y=0;
 num_z=0;
 agg_Rgyr_x=0;
 agg_Rgyr_y=0;
 agg_Rgyr_x=0;
 }
 
 
 
 
 //cout <<num_x <<" " << num_y <<" "<< num_z << endl;
 
 sum_Rgyr_x=0;
 sum_Rgyr_y=0;
 sum_Rgyr_z=0;
 num_x=0;
 num_y=0;
 num_z=0;
 
 
 
 
 
 //cout << "WHAT " <<GyrTens[n][0] << " " << GyrTens[n][1] << " " << GyrTens[n][2] <<endl;
 //cout << "WHAT2 " << a<< " " << b << " " << c << " -> " << (b*b-4*a*c)/(4*a) << endl;
 //cout << "Cell " << n << " " <<"Type " << Type[n] << " " <<Rgyr[n] <<" " <<  aspher[n] << endl;
 //cout <<"Type " << Type[n] << " " <<Rgyr[n] << endl;
 //cout << "Cell " << n << " " <<"Type " << Type[n] << " " <<  aspher[n] << endl;
 
 //cout <<agg_Rgyr_x << " "<<agg_Rgyr_y << " "<<agg_Rgyr_z << endl; //_x: cell type 0; _y: cell type 1; _z: cell type 2
 //cout <<agg_aspher_x << " "<<agg_aspher_y << " "<<agg_aspher_z << endl;
 
 
 
 
 }
 cin.get();
 
 }
 
 
 */
// Now need to see the average shape of the cells... (A, B, C)



