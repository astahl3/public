#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

double** read_fields();
void get_fields(double x1, double x2, double x3);

///////////// Fundamental Constants //////////////////////
double m_p = 1.67262192*1e-27; // mass of a proton
double mu_0 = 4*M_PI*1e-7; //
double e = 1.60217663*1e-19; // elementary charge

////////////////////// Parameters ///////////////////////////

double B_vec[3] = {82.0*1e-9,0.0,-401.0*1e-9}; // background magnetic field vector
double B0 = sqrt(B_vec[0]*B_vec[0] + B_vec[1]*B_vec[1] + B_vec[2]*B_vec[2]); //background magnetic field magnitude
double n0 = 200e6; //background density
double m0 = 18.5*1.67262192*1e-27; //upstream ion mass
double VA0 = B0/sqrt(mu_0*m0*n0);
double t0 = m0/(e*B0);
double x0 = VA0*t0;
double E0 = VA0*B0;
double Rmoon = 1560.8e3;

// Geometric parameters from aikef
const double box_size[3] = {20*Rmoon/x0, 20*Rmoon/x0, 20*Rmoon/x0};
const double obstacle_loc[3] = {0.3001*box_size[0], 0.5001*box_size[1], 0.5001*box_size[2]};
const double dx = 0.705;
const double dy = 0.784;
const double dz = 0.784;
const int meshpoints[3] = {640,576,576};
//const int num_of_lines = meshpoints[0]*meshpoints[1]*meshpoints[2];

//choose the size of the output domain, must be smaller than AIKEF domain
double pos_x_bound = 5; // in moon radii (Rmoon)
double neg_x_bound = -5;
double pos_y_bound = 5;
double neg_y_bound = -5;
double pos_z_bound = 5;
double neg_z_bound = -5;

//choose resolution
double x_res = 0.1; // in moon radii (Rmoon)
double y_res = 0.1;
double z_res = 0.1; 

int num_x_steps = (pos_x_bound - neg_x_bound)/x_res;
int num_y_steps = (pos_y_bound - neg_y_bound)/y_res;
int num_z_steps = (pos_z_bound - neg_z_bound)/z_res;

double x;
double y;
double z;

double B[3];
double E[3];

double** aikef_fields;

int main()
{

  ofstream xy_bxfile;
  ofstream xy_byfile;
  ofstream xy_bzfile;
  ofstream xy_exfile;
  ofstream xy_eyfile;
  ofstream xy_ezfile;

  ofstream xz_bxfile;
  ofstream xz_byfile;
  ofstream xz_bzfile;
  ofstream xz_exfile;
  ofstream xz_eyfile;
  ofstream xz_ezfile;

  ofstream yz_bxfile;
  ofstream yz_byfile;
  ofstream yz_bzfile;
  ofstream yz_exfile;
  ofstream yz_eyfile;
  ofstream yz_ezfile;

  
  // XY Cross Section
  std::string xy_bxfilename = "XY/Bx.txt";
  xy_bxfile.open(xy_bxfilename,std::ios::app);

  std::string xy_byfilename = "XY/By.txt";
  xy_byfile.open(xy_byfilename,std::ios::app);

  std::string xy_bzfilename = "XY/Bz.txt";
  xy_bzfile.open(xy_bzfilename,std::ios::app);

  std::string xy_exfilename = "XY/Ex.txt";
  xy_exfile.open(xy_exfilename,std::ios::app);

  std::string xy_eyfilename = "XY/Ey.txt";
  xy_eyfile.open(xy_eyfilename,std::ios::app);

  std::string xy_ezfilename = "XY/Ez.txt";
  xy_ezfile.open(xy_ezfilename,std::ios::app);

  //XZ Cross Section

  std::string xz_bxfilename = "XZ/Bx.txt";
  xz_bxfile.open(xz_bxfilename,std::ios::app);

  std::string xz_byfilename = "XZ/By.txt";
  xz_byfile.open(xz_byfilename,std::ios::app);

  std::string xz_bzfilename = "XZ/Bz.txt";
  xz_bzfile.open(xz_bzfilename,std::ios::app);

  std::string xz_exfilename = "XZ/Ex.txt";
  xz_exfile.open(xz_exfilename,std::ios::app);

  std::string xz_eyfilename = "XZ/Ey.txt";
  xz_eyfile.open(xz_eyfilename,std::ios::app);

  std::string xz_ezfilename = "XZ/Ez.txt";
  xz_ezfile.open(xz_ezfilename,std::ios::app);

  //YZ Cross Section

  std::string yz_bxfilename = "YZ/Bx.txt";
  yz_bxfile.open(yz_bxfilename,std::ios::app);

  std::string yz_byfilename = "YZ/By.txt";
  yz_byfile.open(yz_byfilename,std::ios::app);

  std::string yz_bzfilename = "YZ/Bz.txt";
  yz_bzfile.open(yz_bzfilename,std::ios::app);

  std::string yz_exfilename = "YZ/Ex.txt";
  yz_exfile.open(yz_exfilename,std::ios::app);

  std::string yz_eyfilename = "YZ/Ey.txt";
  yz_eyfile.open(yz_eyfilename,std::ios::app);

  std::string yz_ezfilename = "YZ/Ez.txt";
  yz_ezfile.open(yz_ezfilename,std::ios::app);

  cout << "about to read the aikef cube" << std::endl;

  // read in the aikef cube
  aikef_fields = read_fields();

  cout << "aikef_fields read in" << std::endl;

  cout << "Starting XY cross section..." << std::endl;
  // XY cross section
  for(int i = 0; i < num_y_steps; i++)
    {
      y = (neg_y_bound + i*y_res)*Rmoon;

      for(int j = 0; j < num_x_steps; j++)
	{

	  x = (neg_x_bound + j*x_res)*Rmoon;

	  cout << "(x,y) = (" << x << "," << y << ")" << std::endl;

	  get_fields(x/x0, y/x0, 0);
	  
	  xy_bxfile << B[0]*B0 << " ";
	  xy_byfile << B[1]*B0 << " ";
	  xy_bzfile << B[2]*B0 << " ";
	  xy_exfile << E[0]*E0 << " ";
	  xy_eyfile << E[1]*E0 << " ";
	  xy_ezfile << E[2]*E0 << " ";
	}

      xy_bxfile << std::endl;
      xy_byfile << std::endl;
      xy_bzfile << std::endl;
      xy_exfile << std::endl;
      xy_eyfile << std::endl;
      xy_ezfile << std::endl;
    }

  cout << "...done" << std::endl;

  cout << "Starting XZ cross section..." << std::endl;

  // XZ cross section
  for(int ii = 0; ii < num_z_steps; ii++)
    {

      z = (neg_z_bound + ii*z_res)*Rmoon;

      for(int jj = 0; jj < num_x_steps; jj++)
	{
	  x = (neg_x_bound + jj*x_res)*Rmoon;

	  cout << "(x,z) = (" << x << "," << z << ")" << std::endl;

	  get_fields(x/x0, 0, z/x0);
	  
	  xz_bxfile << B[0]*B0 << " ";
	  xz_byfile << B[1]*B0 << " ";
	  xz_bzfile << B[2]*B0 << " ";
	  xz_exfile << E[0]*E0 << " ";
	  xz_eyfile << E[1]*E0 << " ";
	  xz_ezfile << E[2]*E0 << " ";
	}

      xz_bxfile << std::endl;
      xz_byfile << std::endl;
      xz_bzfile << std::endl;
      xz_exfile << std::endl;
      xz_eyfile << std::endl;
      xz_ezfile << std::endl;
    }

  cout << "...done" << std::endl;

  cout << "Starting YZ cross section..." << std::endl;

  // YZ cross section
  for(int iii = 0; iii < num_z_steps; iii++)
    {

      z = (neg_z_bound + iii*z_res)*Rmoon;

      for(int jjj = 0; jjj < num_y_steps; jjj++)
	{
	  y = (neg_y_bound + jjj*y_res)*Rmoon;

	  cout << "(y,z) = (" << y << "," << z << ")" << std::endl;

	  get_fields(0, y/x0, z/x0);

	  yz_bxfile << B[0]*B0 << " ";
	  yz_byfile << B[1]*B0 << " ";
	  yz_bzfile << B[2]*B0 << " ";
	  yz_exfile << E[0]*E0 << " ";
	  yz_eyfile << E[1]*E0 << " ";
	  yz_ezfile << E[2]*E0 << " ";
	}

      yz_bxfile << std::endl;
      yz_byfile << std::endl;
      yz_bzfile << std::endl;
      yz_exfile << std::endl;
      yz_eyfile << std::endl;
      yz_ezfile << std::endl;
    }
  
  cout << "...done" << std::endl;

}


//! Routine to Read in AIKEF Field Cubes
double** read_fields()
{
  std::cout << "In routines.cpp" << endl;
  //logfile << "----Reading in AIKEF Cubes----" << endl << endl;
  ifstream line_check("aikef_fields/bzfile.txt");
  string line;

  int num_of_lines = 0;

  while(getline(line_check, line)) // determine number of lines in AIKEF files
    num_of_lines++;
  line_check.close();

  double* Bx = new double[num_of_lines]; //create arrays to hold each field component
  double* By = new double[num_of_lines];
  double* Bz = new double[num_of_lines];



  double* Ex = new double[num_of_lines];
  double* Ey = new double[num_of_lines];
  double* Ez = new double[num_of_lines];

  double** fields = new double*[num_of_lines]; //create array to hold all components

  ifstream BxFile("aikef_fields/bxfile.txt"); //read in files, MAKE SURE THE AIKEF FILES ARE AT THIS PATH
  ifstream ByFile("aikef_fields/byfile.txt");
  ifstream BzFile("aikef_fields/bzfile.txt");

  ifstream ExFile("aikef_fields/exfile.txt");
  ifstream EyFile("aikef_fields/eyfile.txt");
  ifstream EzFile("aikef_fields/ezfile.txt");

  std::cout << "Ifstream commands completed successfully" << std::endl;

  if (BxFile.is_open())
    for (int i = 0; i < num_of_lines; i++)
      {
	BxFile >> Bx[i];

	fields[i] = new double[6];

	fields[i][0] = Bx[i];

      }

  BxFile.close();
  delete[] Bx;

  std::cout << "Bx read successfully" << std::endl;

 if (ByFile.is_open())
    for (int i = 0; i < num_of_lines; i++)
      {
	ByFile >> By[i];

	fields[i][1] = By[i];

      }

  ByFile.close();
  delete[] By;

  std::cout << "By read successfully" << std::endl;

 if (BzFile.is_open())
    for (int i = 0; i < num_of_lines; i++)
      {
	BzFile >> Bz[i];

	fields[i][2] = Bz[i];

      }

 BzFile.close();
 delete[] Bz;

 std::cout << "Bz read successfully" << std::endl;

 if (ExFile.is_open())
    for (int i = 0; i < num_of_lines; i++)
      {
	ExFile >> Ex[i];

	fields[i][3] = Ex[i];

      }

  ExFile.close();
  delete[] Ex;

  std::cout << "Ex read successfully" << std::endl;

 if (EyFile.is_open())
    for (int i = 0; i < num_of_lines; i++)
      {
	EyFile >> Ey[i];

	fields[i][4] = Ey[i];

      }

  EyFile.close();
  delete[] Ey;

  std::cout << "Ey read successfully" << std::endl;

 if (EzFile.is_open())
    for (int i = 0; i < num_of_lines; i++)
      {
	EzFile >> Ez[i];

	fields[i][5] = Ez[i];

      }

  EzFile.close();
  delete[] Ez;

  std::cout << "Ez read successfully" << std::endl;

  std::cout << "----Done----" << std::endl << std::endl;

  return fields;

}


void get_fields(double x1, double x2, double x3)
{ 
  //! Here we differ from the original GENTOo, which interpolates at the border of the box,
  //! using AIKEF values on the box side and homoegeneous values outside,
  //! this approach is only valid if the fields at the edge of the AIKEF box
  //! are safely homogenous, will test

  //! If particle is inside the box, use trilinear interpolation

//! now, shift coordinate system, moving origin from center of moon to lower, left, back corner of AIKEF Cube (-xmax,-ymax,-zmax)
  double pos[3] = {x1+obstacle_loc[0],x2+obstacle_loc[1],x3+obstacle_loc[2]};

  //! now, determine the cell that the particle is in

  int cell1 = pos[0]/dx; //same as x_index in Lucas' code
  int cell2 = pos[1]/dy; // y_index
  int cell3 = pos[2]/dz; // z_index
  int cell[3] = {cell1,cell2,cell3};

  double lambda = pos[0]/dx-cell1; //distance from the cell origin in the x direction
  double mu = pos[1]/dy-cell2; //distance from the cell origin in the y direction
  double nu = pos[2]/dz-cell3; //distance from the cell origin in the z direction

  int i_j_k;
  int ip1_j_k;
  int i_jp1_k;
  int i_j_kp1;
  int ip1_jp1_k;
  int ip1_j_kp1;
  int i_jp1_kp1;
  int ip1_jp1_kp1;
  
  cout << "determining indices..." << std::endl;

  i_j_k = cell1*meshpoints[1]*meshpoints[2]+cell2*meshpoints[2]+cell3; //meshpoints determined in get_cellsize()

  ip1_j_k = (1+cell1)*meshpoints[1]*meshpoints[2] +    cell2 *meshpoints[2] +    cell3;
  i_jp1_k =    cell1 *meshpoints[1]*meshpoints[2] + (1+cell2)*meshpoints[2] +    cell3;
  i_j_kp1 =    cell1 *meshpoints[1]*meshpoints[2] +    cell2 *meshpoints[2] + (1+cell3);
	      
  //! -------------------------------------------
  ip1_jp1_k = (1+cell1)*meshpoints[1]*meshpoints[2] +(1+cell2)*meshpoints[2] +    cell3;
  ip1_j_kp1 = (1+cell1)*meshpoints[1]*meshpoints[2] +   cell2 *meshpoints[2] + (1+cell3);
  i_jp1_kp1 =    cell1 *meshpoints[1]*meshpoints[2] +(1+cell2)*meshpoints[2] + (1+cell3);
	     
  //! -------------------------------------------
  ip1_jp1_kp1 = (1+cell1)*meshpoints[1]*meshpoints[2] + (1+cell2)*meshpoints[2] + (1+cell3);

  ///////////////////////////////////////////////////////

  cout << "...done" << std::endl;

  cout << "i_j_k = " << i_j_k << std::endl
       << "ip1_j_k = " << ip1_j_k << std::endl
       << "i_jp1_k = " << i_jp1_k << std::endl
       << "i_j_kp1 = " << i_j_kp1 << std::endl
       << "ip1_jp1_k = " << ip1_jp1_k << std::endl
       << "ip1_j_kp1 = " << ip1_j_kp1 << std::endl
       << "i_jp1_kp1 = " << i_jp1_kp1 << std::endl
       << "ip1_jp1_kp1 = " << ip1_jp1_kp1 << std::endl;

  B[0] = aikef_fields[i_j_k][0]*(1-lambda)*(1-mu)*(1-nu)
    + aikef_fields[ip1_j_k][0]*(lambda)*(1-mu)*(1-nu)
    + aikef_fields[i_jp1_k][0]*(1-lambda)*(mu)*(1-nu)
    + aikef_fields[i_j_kp1][0]*(1-lambda)*(1-mu)*(nu)
    + aikef_fields[ip1_jp1_k][0]*(lambda)*(mu)*(1-nu)
    + aikef_fields[ip1_j_kp1][0]*(lambda)*(1-mu)*(nu)
    + aikef_fields[i_jp1_kp1][0]*(1-lambda)*(mu)*(nu)
    + aikef_fields[ip1_jp1_kp1][0]*(lambda)*(mu)*(nu);

  B[1] = aikef_fields[i_j_k][1]*(1-lambda)*(1-mu)*(1-nu)
    + aikef_fields[ip1_j_k][1]*(lambda)*(1-mu)*(1-nu)
    + aikef_fields[i_jp1_k][1]*(1-lambda)*(mu)*(1-nu)
    + aikef_fields[i_j_kp1][1]*(1-lambda)*(1-mu)*(nu)
    + aikef_fields[ip1_jp1_k][1]*(lambda)*(mu)*(1-nu)
    + aikef_fields[ip1_j_kp1][1]*(lambda)*(1-mu)*(nu)
    + aikef_fields[i_jp1_kp1][1]*(1-lambda)*(mu)*(nu)
    + aikef_fields[ip1_jp1_kp1][1]*(lambda)*(mu)*(nu);

  B[2] = aikef_fields[i_j_k][2]*(1-lambda)*(1-mu)*(1-nu)
    + aikef_fields[ip1_j_k][2]*(lambda)*(1-mu)*(1-nu)
    + aikef_fields[i_jp1_k][2]*(1-lambda)*(mu)*(1-nu)
    + aikef_fields[i_j_kp1][2]*(1-lambda)*(1-mu)*(nu)
    + aikef_fields[ip1_jp1_k][2]*(lambda)*(mu)*(1-nu)
    + aikef_fields[ip1_j_kp1][2]*(lambda)*(1-mu)*(nu)
    + aikef_fields[i_jp1_kp1][2]*(1-lambda)*(mu)*(nu)
    + aikef_fields[ip1_jp1_kp1][2]*(lambda)*(mu)*(nu);

  E[0] = aikef_fields[i_j_k][3]*(1-lambda)*(1-mu)*(1-nu)
    + aikef_fields[ip1_j_k][3]*(lambda)*(1-mu)*(1-nu)
    + aikef_fields[i_jp1_k][3]*(1-lambda)*(mu)*(1-nu)
    + aikef_fields[i_j_kp1][3]*(1-lambda)*(1-mu)*(nu)
    + aikef_fields[ip1_jp1_k][3]*(lambda)*(mu)*(1-nu)
    + aikef_fields[ip1_j_kp1][3]*(lambda)*(1-mu)*(nu)
    + aikef_fields[i_jp1_kp1][3]*(1-lambda)*(mu)*(nu)
    + aikef_fields[ip1_jp1_kp1][3]*(lambda)*(mu)*(nu);

  E[1] = aikef_fields[i_j_k][4]*(1-lambda)*(1-mu)*(1-nu)
    + aikef_fields[ip1_j_k][4]*(lambda)*(1-mu)*(1-nu)
    + aikef_fields[i_jp1_k][4]*(1-lambda)*(mu)*(1-nu)
    + aikef_fields[i_j_kp1][4]*(1-lambda)*(1-mu)*(nu)
    + aikef_fields[ip1_jp1_k][4]*(lambda)*(mu)*(1-nu)
    + aikef_fields[ip1_j_kp1][4]*(lambda)*(1-mu)*(nu)
    + aikef_fields[i_jp1_kp1][4]*(1-lambda)*(mu)*(nu)
    + aikef_fields[ip1_jp1_kp1][4]*(lambda)*(mu)*(nu);

  E[2] = aikef_fields[i_j_k][5]*(1-lambda)*(1-mu)*(1-nu)
    + aikef_fields[ip1_j_k][5]*(lambda)*(1-mu)*(1-nu)
    + aikef_fields[i_jp1_k][5]*(1-lambda)*(mu)*(1-nu)
    + aikef_fields[i_j_kp1][5]*(1-lambda)*(1-mu)*(nu)
    + aikef_fields[ip1_jp1_k][5]*(lambda)*(mu)*(1-nu)
    + aikef_fields[ip1_j_kp1][5]*(lambda)*(1-mu)*(nu)
    + aikef_fields[i_jp1_kp1][5]*(1-lambda)*(mu)*(nu)
    + aikef_fields[ip1_jp1_kp1][5]*(lambda)*(mu)*(nu);
  
}
