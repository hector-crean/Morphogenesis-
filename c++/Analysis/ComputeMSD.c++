#include<iostream>
#include<math.h>
#include<stdlib.h>
#include<fstream>
#include<sstream>
#include<vector>
#include<ctime>
using namespace std;

const int Tmax=10000;
const int Nmax=1000;
int t;

//
double COM_wrapped[Tmax][Nmax][2];
double COM[Tmax][Nmax][2];
int PBC[2];

//
double msd[Nmax][Tmax];
double msd_ave[Tmax];

//START PROGRAM
int main(int argc, char* argv[]){

PBC[0]=PBC[1]=0;
for(int tau=0; tau<Tmax; tau++)for(int nc=0; nc<Nmax; nc++)msd[nc][tau]=0;
for(int tau=0; tau<Tmax; tau++)msd_ave[tau]=0;

cout<< "Type: 1. name input 2. How many timesteps? 3. Lagtime? 4. How many cells? 5. Size Grid? 6. output" <<endl;
int Nframes=atoi(argv[2]);
int lagtime=atoi(argv[3]);
int Ncells=atoi(argv[4]);
int L=atoi(argv[5]);

//DEFINE INPUT
ifstream read;
read.open(argv[1]);

// CREATE NAMES OF FILE FOR OUTPUT //
stringstream writeMSDF;
writeMSDF <<"MSD_" << argv[6] <<".dat";
ofstream writeMSD(writeMSDF.str().c_str());

stringstream writeCOMF;
writeCOMF <<"COM_" << argv[6] <<".dat";
ofstream writeCOM(writeCOMF.str().c_str());


for(int t=0; t<Nframes; t++){
string dummy;
getline(read,dummy); //this is to read the #time
//cout <<"d"<< dummy <<endl;

//read coordinates
for(int nc=0; nc<Ncells; nc++){
int id;
double x,y;
read >> id >> x >> y;
//cout << id << " " << x << " "<< y<<endl; cin.get();
COM_wrapped[t][id][0]=x;
COM_wrapped[t][id][1]=y;
}

//read 2 empty lines
getline(read,dummy);
//cout <<"d"<< dummy <<endl;
getline(read,dummy);
//cout <<"d"<< dummy <<endl;
getline(read,dummy);
//cout <<"d"<< dummy <<endl;
}

/////////////////
// END READING //
/////////////////

/////////////////
// RESOLVE PBC //
/////////////////
//when COM jumps of a value +/-L then add -/+L and continue

for(int nc=0; nc<Ncells; nc++){
PBC[0]=0;
PBC[1]=0;
 writeCOM << "#Cell " << nc+1<<endl;
    for(int t=0; t<Nframes; t++){
        if(t==0){
        COM[t][nc][0]=COM_wrapped[t][nc][0];
        COM[t][nc][1]=COM_wrapped[t][nc][1];
        }
        if(t>0){
            if(COM_wrapped[t][nc][0]-COM_wrapped[t-1][nc][0]>L/2.)PBC[0]--;
            if(COM_wrapped[t][nc][1]-COM_wrapped[t-1][nc][1]>L/2.)PBC[1]--;
            if(COM_wrapped[t][nc][0]-COM_wrapped[t-1][nc][0]<-L/2.)PBC[0]++;
            if(COM_wrapped[t][nc][1]-COM_wrapped[t-1][nc][1]<-L/2.)PBC[1]++;
            
            COM[t][nc][0]=COM_wrapped[t][nc][0]+PBC[0]*L;
            COM[t][nc][1]=COM_wrapped[t][nc][1]+PBC[1]*L;;
            
        }
    writeCOM<< t << " " <<  nc << " " << COM_wrapped[t][nc][0] << " " << COM_wrapped[t][nc][1]<< " " << COM[t][nc][0] << " " << COM[t][nc][1]<<endl ;
    }
    
    writeCOM<<endl;
    writeCOM<<endl;
}


//////////////////////////////////////
// COMPUTE DIFFERENCE AT LAG TIME T //
//////////////////////////////////////
writeMSD<<"#AVERAGE" <<endl;
for(int tau=0; tau<Nframes; tau++){
cout << "doing tau = " << tau <<endl;
    for(int nc=0; nc<Ncells; nc++){
    int Tave=Nframes-tau;
        for(int t=0; t<Tave; t++){
            msd[nc][tau]+=pow(COM[t+tau][nc][0]-COM[t][nc][0],2.0)+pow(COM[t+tau][nc][1]-COM[t][nc][1],2.0);
            }
msd[nc][tau]=msd[nc][tau]*1.0/Tave;
msd_ave[tau]+=msd[nc][tau];
}
msd_ave[tau]=msd_ave[tau]*1.0/Ncells;
writeMSD << tau*lagtime << " " << msd_ave[tau] <<endl;
}

writeMSD <<endl;
writeMSD <<endl;

for(int nc=0; nc<Ncells; nc++){
writeMSD<<"#Cell " << nc+1 <<endl;
    for(int tau=0; tau<Nframes; tau++){
    writeMSD << tau*lagtime << " " << msd[nc][tau] <<endl;
    }
writeMSD<<endl;
writeMSD<<endl;
}



return 0;
}

