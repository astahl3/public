#include <cstdlib>
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <math.h>
using namespace std;

int readhybrid(int,float*,float*,float*,float*,float*,float*,float*,float*,float*,float*,float*,float*);

int main(int argc, char** argv) {
    int count;
    float * bx;
    float * by;
    float * bz;
    float * ex;
    float * ey;
    float * ez;
    float * lx;
    float * ly;
    float * lz;
    float x0hc,v0hc,b0hc;
    
    const int npoints = 640*576*576; //meshpoints x*y*z;
    
    lx = (float*) malloc(npoints * sizeof (float));
    ly = (float*) malloc(npoints * sizeof (float));
    lz = (float*) malloc(npoints * sizeof (float));
    bx = (float*) malloc(npoints * sizeof (float));
    by = (float*) malloc(npoints * sizeof (float));
    bz = (float*) malloc(npoints * sizeof (float));
    ex = (float*) malloc(npoints * sizeof (float));
    ey = (float*) malloc(npoints * sizeof (float));
    ez = (float*) malloc(npoints * sizeof (float));
    int dummy = readhybrid(0, &*lx, &*ly, &*lz, &*bx, &*by, &*bz, &*ex, &*ey, &*ez, &x0hc, &b0hc, &v0hc);
    
    //int numgridsim = readhybrid(1,&*lx, &*ly, &*lz, &*bx, &*by, &*bz, &*ex, &*ey, &*ez, &x0hc, &b0hc, &v0hc);
    //int check=mkdir("aikef_fields");

    ofstream auxfile("aikef_fields/auxfile.txt");
    auxfile << "b0= " << b0hc << endl << "v0= " << v0hc << endl << "x0= " <<x0hc << endl << "numpoints= " <<npoints;
    auxfile.close();
    
    ofstream lxfile("aikef_fields/lxfile.txt");
    ofstream lyfile("aikef_fields/lyfile.txt");
    ofstream lzfile("aikef_fields/lzfile.txt");
    ofstream bxfile("aikef_fields/bxfile.txt");
    ofstream byfile("aikef_fields/byfile.txt");
    ofstream bzfile("aikef_fields/bzfile.txt");
    ofstream exfile("aikef_fields/exfile.txt");
    ofstream eyfile("aikef_fields/eyfile.txt");
    ofstream ezfile("aikef_fields/ezfile.txt");
    /*ofstream bxnTfile("aikef_fields/bxnTfile.txt");
    ofstream bynTfile("aikef_fields/bynTfile.txt");
    ofstream bznTfile("aikef_fields/bznTfile.txt");
    ofstream exmVfile("aikef_fields/exmVfile.txt");
    ofstream eymVfile("aikef_fields/eymVfile.txt");
    ofstream ezmVfile("aikef_fields/ezmVfile.txt");
    ofstream lxREfile("aikef_fields/lxREfile.txt");
    ofstream lyREfile("aikef_fields/lyREfile.txt");
    ofstream lzREfile("aikef_fields/lzREfile.txt");
    ofstream Rfile("aikef_fields/Rfile.txt");
    ofstream allfiles("aikef_fields/allfiles.csv");*/
    cout << endl << "Datapoints left: " << endl;
    
    float re=1560800.;
    //how cooper "wanted" the output
    //allfiles << "Bx.nT,By.nT,Bz.nT,Ex.mV,Ey.mV,Ez.mV,Lx.RE,Ly.RE,Lz.RE,R.RE" << endl;

    for (count = 0; count < npoints; count++) {
        lxfile << lx[count] << endl;
        lyfile << ly[count] << endl;
        lzfile << lz[count] << endl;
        bxfile << bx[count] << endl;
        byfile << by[count] << endl;
        bzfile << bz[count] << endl;
        exfile << ex[count] << endl;
        eyfile << ey[count] << endl;
        ezfile << ez[count] << endl;
	/*bxnTfile << bx[count]*b0hc*1.e9 << endl;
	bynTfile << by[count]*b0hc*1.e9 << endl;
	bznTfile << bz[count]*b0hc*1.e9 << endl;
	exmVfile << ex[count]*b0hc*v0hc*1.e3 << endl;
	eymVfile << ey[count]*b0hc*v0hc*1.e3 << endl;
	ezmVfile << ez[count]*b0hc*v0hc*1.e3 << endl;
	lxREfile << lx[count]*x0hc/re << endl;
	lyREfile << lx[count]*x0hc/re << endl;
	lzREfile << lx[count]*x0hc/re << endl;
	Rfile << sqrt(pow(lx[count]*x0hc/re,2)+pow(ly[count]*x0hc/re,2)+pow(lz[count]*x0hc/re,2)) <<endl;
	allfiles << bx[count]*b0hc*1.e9 << ',' << by[count]*b0hc*1.e9 << ',' << bz[count]*b0hc*1.e9 << ',' << ex[count]*b0hc*v0hc*1.e3 << ',' << ey[count]*b0hc*v0hc*1.e3 << ',' << ez[count]*b0hc*v0hc*1.e3 << ',' << lx[count]*x0hc/re << ',' << ly[count]*x0hc/re << ',' << lz[count]*x0hc/re << ',' << sqrt(pow(lx[count]*x0hc/re,2)+pow(ly[count]*x0hc/re,2)+pow(lz[count]*x0hc/re,2)) <<endl;
	*/
	if(count%500000==0)
	  cout << npoints-count << endl;
	
    }
    lxfile.close();
    lyfile.close();
    lzfile.close();
    bxfile.close();
    byfile.close();
    bzfile.close();
    exfile.close();
    eyfile.close();
    ezfile.close();
    /*allfiles.close();
    bxnTfile.close();
    bynTfile.close();
    bznTfile.close();
    exmVfile.close();
    eymVfile.close();
    ezmVfile.close();
    lxREfile.close(); ofstream auxfile("aikef_fields/auxfile.txt");
    auxfile << "b0= " << b0hc << endl << "v0= " << v0hc << endl << "x0= " <<x0hc << endl << "numpoints= " <<npoints;
    auxfile.close();
    lyREfile.close();
    lzREfile.close();
    Rfile.close();*/
    

    free(lx);
    free(ly);
    free(lz);
    free(bx);
    free(by);
    free(bz);
    free(ex);
    free(ey);
    free(ez);

    return 0;
}

