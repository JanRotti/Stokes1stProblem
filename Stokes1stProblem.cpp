//--------------------------------------------------------------------------------------//
// --------------------------LATTICE BOLTZMANN ALGORITHM--------------------------------//
//																						//
// CASE DESCRIPTION: Inkompressible Stokes First Problem 2D 							//
// AUTHOR: Jan Rottmayer																//
// INSTITUTION:	TU Kaiserslautern														//
// DATE: 23.07.2021																		//
//																						//
//--------------------------------------------------------------------------------------//

#include <omp.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <string>

using namespace std;

int main()
{
	//--------------------------------------------------------------------------------------//
	// 									  CONFIGURATION										//
	//--------------------------------------------------------------------------------------//

	//------------------------------- GEOMETRIC DEFINITION ---------------------------------//
	double height = 2;										//domain height in [m]
	double length = 0.5;									//domain width in [m]
	double deltaX = 0.01;									//grid spacing in [m]
	double epsilon = 1e-8;									//geometrical tolerance
	long nx = ceil(length / deltaX + epsilon) + 1;			//number of nodes in x direction
	long ny = ceil(height / deltaX + epsilon) + 1;			//number of nodes in y direction

	//--------------------------------- FLUID PARAMETERS -----------------------------------//
	double rho = 900;										//fluid density [kg/m^3] 
	double viscosity = 3e-3;								//kinematic fluid viscosity

	//--------------------------------- FLOW PARAMETERS ------------------------------------//						 
	double wallVelocity = 1e-3;								//wall velocity at bottom plate
	double maxExpectedVelocity = wallVelocity;				//max. expected velocity
	
	//-------------------------------- TIME CONFIGURATION ----------------------------------//
	long timeSteps = 10000;									//time steps 
	long writeInterval = 500;								//writing intervall
	double deltaT = 5e-5; //5e-5;							//time step
	double finalTime = deltaT * timeSteps;					//final time


	//-------------------------------- COMPUTATION SETUP -----------------------------------//
	int cpus = 16;											//number of virtual cpu cores 
	omp_set_num_threads(cpus);								//omp computation cores

	//-------------------------------- LATTICE VARIABLES -----------------------------------//
	double xsi0 = deltaX / deltaT;							//molecular velocity
	double machNumber = sqrt(3.) * deltaT * 
						maxExpectedVelocity / deltaX;		//mach number
	double speedOfSound = xsi0 / sqrt(3.);					//isothermal speed of sound
	double omega = pow(speedOfSound, 2) * deltaT / 
				  (viscosity + 0.5 * pow(speedOfSound, 2) 
				  * deltaT); 								//relaxation parameter

	//--------------------------------- TERMINAL OUTPUT ------------------------------------//
	cout.width(40);

	cout << "DOMAIN DEFINITION" << endl;
	cout << "Domain Height:\t\t" <<  height << endl;
	cout << "Domain Length:\t\t" <<  length << endl;
	cout << "Grid Resolution:\t" << deltaX << endl;
	cout << "Grid Size:\t\t" <<  nx << "x" << ny << endl;

	cout << "\n" <<  "FLUID MODEL" << endl;
	cout << "Density:\t\t" << rho << endl;
	cout << "Viscosity:\t\t" << viscosity << endl;

	cout << "\n" << "TIME SETTING" << endl;
	cout << "Number of Timesteps:\t" << timeSteps << endl;
	cout << "Timestep:\t\t" << deltaT << endl;
	cout << "Max Simulation Time:\t" << finalTime << endl;

	cout << "\n" << "STABILITY/VALIDITY" << endl;
	cout << "Mach number:\t\t" << machNumber << endl;
	cout << "Isothermal SoS:\t\t" << speedOfSound << endl;
	cout << "Omega:\t\t\t" << omega << endl;

	cout << "\n" << "SIMULATION PROGRESS" << endl;
	
	//--------------------------------------------------------------------------------------//
	// 							DATA STRUCTURE INITIALIZATION								//
	//--------------------------------------------------------------------------------------//
	
	// FLUID DENSITY -> first index: x; second index: y
	double** fluidDensity = new double* [nx]; 
	for (int k = 0; k < nx; k++)
	{
		fluidDensity[k] = new double[ny];
		for (int l = 0; l < ny; l++)
		{
			//fixed density initialization
			fluidDensity[k][l] = rho; 
		}//end l
	}//end k

	// FLUID VELOCITY -> first index: x; second index: y;
	// third index: velocity direction [0->x, 1->y]
	double*** fluidVelocity = new double** [nx]; 
	for (int k = 0; k < nx; k++)
	{
		fluidVelocity[k] = new double* [ny];
		for (int l = 0; l < ny; l++)
		{
			fluidVelocity[k][l] = new double[2];
			//resting fluid initialization
			fluidVelocity[k][l][0] = 0.; 
			fluidVelocity[k][l][1] = 0.;
		}//end l
	}//end k

	// LATTICE WEIGHTS
	// weigths of lattice directions for equilibrium distribution
	double weigths[9] = { 1./9., 1./9., 1./9., 1./9., 1./36., 1./36., 1./36., 1./36., 4./9.}; 

	// LOCAL DISTRIBUTION -> first index: x; second index: y;
	// third index: transport(0)/collision(1); fourth index: lattice direction
	double**** distribution = new double*** [nx]; 
	for (int k = 0; k < nx; k++)
	{
		distribution[k] = new double** [ny];
		for (int l = 0; l < ny; l++)
		{
			distribution[k][l] = new double* [2];
			distribution[k][l][0] = new double[9];
			distribution[k][l][1] = new double[9];
			
			// iterate lattice directions
			for (int i = 0; i < 9; i++) 
			{
				// initialize with rho and u=0;
				distribution[k][l][0][i] = weigths[i] * rho; //transport array
				distribution[k][l][1][i] = weigths[i] * rho; //collision array
			}//end i
		}//end l
	}//end k

	// LATTICE VELOCITIES -> first array index: lattice direction; 
	// second index velocity direction (0:x-dir;1:y-dir)
	double xsi[9][2];
	xsi[0][0] = xsi0; xsi[1][0] = 0; xsi[2][0] = -xsi0; xsi[3][0] = 0; xsi[4][0] = xsi0;
	xsi[5][0] = -xsi0; xsi[6][0] = -xsi0; xsi[7][0] = xsi0; xsi[8][0] = 0;
	xsi[0][1] = 0; xsi[1][1] = xsi0; xsi[2][1] = 0; xsi[3][1] = -xsi0; xsi[4][1] = xsi0;
	xsi[5][1] = xsi0; xsi[6][1] = -xsi0; xsi[7][1] = -xsi0; xsi[8][1] = 0;


	//--------------------------------------------------------------------------------------//
	// 								WRITING OUTPUT FUNCTION									//
	//--------------------------------------------------------------------------------------//
	auto writeResults = [](long t, double deltaT, long nx, long ny, double deltaX, double** fluidDensity, double speedOfSound, double*** fluidVelocity)
	{
		string filename("results" + to_string(static_cast<long>(t)) + ".vtk");

		ofstream results(filename.c_str());

		if (!results)
		{
			cerr << "Error opening output file " << endl;
			return 0;
		}
		else
		{
			results << setiosflags(ios::left | ios::showpoint | ios::fixed) << "# vtk DataFile Version 2.0" << endl;
			results << "Lattice Boltzmann solution at timestep " << t * deltaT << endl;
			results << "ASCII" << endl;
			results << "DATASET STRUCTURED_POINTS" << endl;
			results << "DIMENSIONS " << nx << " " << ny << " " << 1 << endl;
			results << "SPACING " << deltaX << " " << deltaX << " " << 0 << endl;
			results << "ORIGIN 0 0 0" << endl;
			results << "POINT_DATA " << nx * ny << endl;
			results << "SCALARS pressure float 1" << endl;
			results << "LOOKUP_TABLE default" << endl;

			for (int l = 0; l < ny; l++)
				for (int k = 0; k < nx; k++)
					results << setprecision(8) << fluidDensity[k][l] * pow(speedOfSound, 2) << endl;

			results << "SCALARS density float 1" << endl;
			results << "LOOKUP_TABLE default" << endl;

			for (int l = 0; l < ny; l++)
				for (int k = 0; k < nx; k++)
					results << setprecision(8) << fluidDensity[k][l] << endl;

			results << "VECTORS velocity float" << endl;

			for (int l = 0; l < ny; l++)
				for (int k = 0; k < nx; k++)
					results << setprecision(8) << fluidVelocity[k][l][0] << " " << setprecision(8) << fluidVelocity[k][l][1] << " " << 0.0 << endl;
		}
		results.close();
	};//end writeResults


	//--------------------------------------------------------------------------------------//
	// 							LATTICE BOLTZMANN ALGORITHM									//
	//--------------------------------------------------------------------------------------//
	
	// error history
	fstream errorFile;
	// open error file
	remove("history.txt");
	errorFile.open("history.txt");
	
	// start time
	long t = 0;

	//------------------------------------ TIME LOOP ---------------------------------------//
	for (t = 0; t < timeSteps; t++)
	{
	
		//---------------------------------- TRANSPORT STEP --------------------------------//
		#pragma omp parallel for
		// iterate x-direction
		for (int k = 1; k < nx - 1; ++k) 
		{	
			// iterate y-direction
			for (int l = 1; l < ny - 1; ++l) 
			{
				// iterate lattice directions
				for (int i = 0; i < 9; ++i) { 
					// initialize neighbor coordinates
					int kNN = 0, lNN = 0;
					// directional push coordinates
					switch (i) {
					case 0:
					{
						kNN = k + 1;
						lNN = l;
						break;
					}
					case 1:
					{
						kNN = k;
						lNN = l + 1;
						break;
					}
					case 2:
					{
						kNN = k - 1;
						lNN = l;
						break;
					}
					case 3:
					{
						kNN = k;
						lNN = l - 1;
						break;
					}
					case 4:
					{
						kNN = k + 1;
						lNN = l + 1;
						break;
					}
					case 5:
					{
						kNN = k - 1;
						lNN = l + 1;
						break;
					}
					case 6:
					{
						kNN = k - 1;
						lNN = l - 1;
						break;
					}
					case 7:
					{
						kNN = k + 1;
						lNN = l - 1;
						break;
					}
					case 8:
					{
						kNN = k;
						lNN = l;
						break;
					}
					} //end switch
					// push transport of distributions
					distribution[kNN][lNN][0][i] = distribution[k][l][1][i]; 
				}//end i
			}//end l
		}//end k


		//-------------------------------- PERIODIC BOUNDARY -------------------------------//
		// right boundary transport
		#pragma omp parallel for
		// iterate y-direction
		for (int l = 1; l < ny - 1; ++l) { 
			// iterate lattice directions
			for (int i = 0; i < 9; ++i) { 
				// initialize neighbor coordinates
				int kNN = 0, lNN = 0;
				// fix x to right boundary
				int k = nx - 1;
				// directional push coordinates
				switch (i) {
				case 0:
				{
					kNN = 0;
					lNN = l;
					break;
				}
				case 1:
				{
					kNN = k;
					lNN = l + 1;
					break;
				}
				case 2:
				{
					kNN = k - 1;
					lNN = l;
					break;
				}
				case 3:
				{
					kNN = k;
					lNN = l - 1;
					break;
				}
				case 4:
				{
					kNN = 0;
					lNN = l + 1;
					break;
				}
				case 5:
				{
					kNN = k - 1;
					lNN = l + 1;
					break;
				}
				case 6:
				{
					kNN = k - 1;
					lNN = l - 1;
					break;
				}
				case 7:
				{
					kNN = 0;
					lNN = l - 1;
					break;
				}
				case 8:
				{
					kNN = k;
					lNN = l;
					break;
				}
				} //end switch
				// push transport of right boundary distribution
				distribution[kNN][lNN][0][i] = distribution[k][l][1][i]; 
			}//end i
		}//end l

		// left boundary transport
		#pragma omp parallel for
		// iterate y-direction
		for (int l = 1; l < ny - 1; l++) {
			// iterate lattice directions
			for (int i = 0; i < 9; ++i) { 
				// initialize neighbor coordinates
				int kNN = 0, lNN = 0;
				// fix x to left boundary
				int k = 0;
				// directional push coordinates
				switch (i) {
				case 0:
				{
					kNN = k + 1;
					lNN = l;
					break;
				}
				case 1:
				{
					kNN = k;
					lNN = l + 1;
					break;
				}
				case 2:
				{
					kNN = nx - 1;
					lNN = l;
					break;
				}
				case 3:
				{
					kNN = k;
					lNN = l - 1;
					break;
				}
				case 4:
				{
					kNN = k + 1;
					lNN = l + 1;
					break;
				}
				case 5:
				{
					kNN = nx - 1;
					lNN = l + 1;
					break;
				}
				case 6:
				{
					kNN = nx - 1;
					lNN = l - 1;
					break;
				}
				case 7:
				{
					kNN = k + 1;
					lNN = l - 1;
					break;
				}
				case 8:
				{
					kNN = k;
					lNN = l;
					break;
				}
				}//end switch
				// push transport of left boundary distribution
				distribution[kNN][lNN][0][i] = distribution[k][l][1][i]; 
			}//end i
		} //end l


		//------------------------ BOUNCE BACK WITH WALL VELOCITY --------------------------//
		// bottom wall --> half way bounce back with fixed wall velocity
		#pragma omp parallel for
		// iterate x-direction
		for (int k = 1; k < nx - 1; ++k) {
			// fix y to bottom boundary
			int l = 0; 
			// fixed wall velocity factor
			double velocityTerm = 2  * wallVelocity / pow(speedOfSound, 2);
			// bounce back of wall distributions --> negative y-direction
			distribution[k][l + 1][0][1]	 = distribution[k][l][0][3] - fluidDensity[k][l + 1] * weigths[3] * xsi[3][0] * velocityTerm;		//i = 1	
			distribution[k + 1][l + 1][0][4] = distribution[k][l][0][6] - fluidDensity[k + 1][l + 1] * weigths[6] * xsi[6][0] * velocityTerm;	//i = 4
			distribution[k - 1][l + 1][0][5] = distribution[k][l][0][7] - fluidDensity[k -1][l + 1] * weigths[7] * xsi[7][0] * velocityTerm;	//i = 5
		}

		// bounce back corner node(0,0)
		distribution[0][1][0][1]	  = distribution[0][0][0][3] - 2 * fluidDensity[0][1] * weigths[3] * xsi[3][0] * wallVelocity / pow(speedOfSound, 2);	  //i = 1	
		distribution[1][1][0][4]	  = distribution[0][0][0][6] - 2 * fluidDensity[1][1] * weigths[6] * xsi[6][0] * wallVelocity / pow(speedOfSound, 2);	  //i = 4
		distribution[nx - 1][1][0][5] = distribution[0][0][0][7] - 2 * fluidDensity[nx - 1][1] * weigths[7] * xsi[7][0] * wallVelocity / pow(speedOfSound, 2);//i = 5
		
		
		// bounce back corner node(nx-1,0)
		distribution[nx - 1][1][0][1] = distribution[nx - 1][0][0][3] - 2 * fluidDensity[nx - 1][1] * weigths[3] * xsi[3][0] * wallVelocity / pow(speedOfSound, 2);//i = 1
		distribution[0][1][0][4]	  = distribution[nx - 1][0][0][6] - 2 * fluidDensity[0][1] * weigths[6] * xsi[6][0] * wallVelocity / pow(speedOfSound, 2);	   //i = 4
		distribution[nx - 2][1][0][5] = distribution[nx - 1][0][0][7] - 2 * fluidDensity[nx - 2][1] * weigths[7] * xsi[7][0] * wallVelocity / pow(speedOfSound, 2);//i = 5
		

		//------------------------- SPECULAR REFLECTION BOUNDARY ---------------------------//
		//top wall --> specular reflection
		#pragma omp parallel for
		// iterate x-direction
		for (int k = 1; k < nx - 1; ++k) {
			// fix y to top boundary
			int l = ny - 1; 
			// specular reflection of wall distributions --> positive y-direction
			distribution[k][l - 1][0][3]	 = distribution[k][l][0][1];  //i = 3	
			distribution[k - 1][l - 1][0][6] = distribution[k][l][0][5];  //i = 6
			distribution[k + 1][l - 1][0][7] = distribution[k][l][0][4];  //i = 7
		}

		// specular reflection corner node(0,ny-1)
		distribution[0][ny - 2][0][3]	   = distribution[0][ny - 1][0][1]; //i = 3
		distribution[nx - 1][ny - 2][0][6] = distribution[0][ny - 1][0][5]; //i = 6
		distribution[1][ny - 2][0][7]	   = distribution[0][ny - 1][0][4]; //i = 7
		
		// specular reflection corner node(nx-1,ny-1)
		distribution[nx - 1][ny - 2][0][3] = distribution[nx - 1][ny - 1][0][1]; //i = 3
		distribution[nx - 2][ny - 2][0][6] = distribution[nx - 1][ny - 1][0][5]; //i = 6
		distribution[0][ny - 2][0][7]	   = distribution[nx - 1][ny - 1][0][4]; //i = 7
		

		//---------------------------- COMPUTATION OF MOMENTS ------------------------------//
		#pragma omp parallel for
		// iterate x-direction
		for (int k = 0; k <= nx - 1; k++) {
			// iterate y-direction
			for (int l = 1; l < ny - 1; l++) {
				// resetting moments
				fluidDensity[k][l] = 0;
				fluidVelocity[k][l][0] = 0;
				fluidVelocity[k][l][1] = 0;
				// iterate lattice direction
				for (int i = 0; i < 9; i++) { 
					// density
					fluidDensity[k][l] += distribution[k][l][0][i]; 
					// fluid momentum in x and y direction
					fluidVelocity[k][l][0] += xsi[i][0] * distribution[k][l][0][i];
					fluidVelocity[k][l][1] += xsi[i][1] * distribution[k][l][0][i]; 
				}//end i

				// fluid velocity from momentum
				fluidVelocity[k][l][0] = fluidVelocity[k][l][0] / fluidDensity[k][l];
				fluidVelocity[k][l][1] = fluidVelocity[k][l][1] / fluidDensity[k][l];
			}//end l
		}//end k


		//-------------------------------- COLLISION STEP ----------------------------------//
		#pragma omp parallel for
		// iterate x-direction
		for (int k = 0; k <= nx - 1; k++) {
			// iterate y-direction
			for (int l = 1; l < ny - 1; l++) { 
				// iterate lattice direction
				for (int i = 0; i < 9; i++) { 
					// compute temporary skalar product terms
					double dotProductXsiVelo = xsi[i][0] * fluidVelocity[k][l][0] + xsi[i][1] * fluidVelocity[k][l][1];
					double dotProductVelo = pow(fluidVelocity[k][l][0], 2) + pow(fluidVelocity[k][l][1], 2);
					// compute equilibrium distribution
					double eqDistribution = weigths[i] * fluidDensity[k][l] * (1 + dotProductXsiVelo / pow(speedOfSound, 2) + pow(dotProductXsiVelo, 2) / (2 * pow(speedOfSound, 4)) - dotProductVelo / (2 * pow(speedOfSound, 2)));
					// update distribution with Bhatnagar-Gross-Krook model 
					distribution[k][l][1][i] = distribution[k][l][0][i] + omega * (eqDistribution - distribution[k][l][0][i]); 
				}//end i
			}//end l
		}//end k


		//------------------------------------ OUTPUT --------------------------------------//
		// file output 
		if (t % writeInterval == 0) {
			writeResults(t, deltaT, nx, ny, deltaX, fluidDensity, speedOfSound, fluidVelocity);
		}
		// progress bar and analytical error
		float progress = t / float(timeSteps);

		//--------------------------- CONVERGENCE CRITERIA -----------------------------//
		// deviation to analytical solution
		int k = nx / 2;
		double analyticalSolution = 0;
		// compute y coordinates
		double error = 0;
		for(int l = 1; l < ny - 1;l++){
			error += abs(fluidVelocity[k][l][0] - (wallVelocity - wallVelocity * erf(deltaX * (float(l) - 0.5) / (2 * sqrt(viscosity * (t * deltaT))))));
		}
		// normalize error
		error /= (ny - 2);

		// write to error file
		errorFile << (t * deltaT) << ", " << error << endl;

		// update on every progress percent
		if (t %  100 == 0) {
			int barWidth = 40;
			cout << "[";
			int pos = barWidth * progress;
			for (int i = 0; i < barWidth; ++i) {
				if (i < pos) { 
					cout << "=";
				}
				else if (i == pos) {
					cout << ">";
				}
				else {
					cout << " ";
				}
			}//end for
			cout << "] " << int(progress * 100.0) << " %\t";
			cout << "Analytical Error:\t" << error << "      " << "\r";
			cout.flush();
		}//end if
	}//end time loop 
	
	 
	 //write final result
	writeResults(t, deltaT, nx, ny, deltaX, fluidDensity, speedOfSound, fluidVelocity);
}//end main