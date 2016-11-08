#include <iostream>
#include <time.h>

#include "../include/rand_num.h"
#include "../include/force.h"
#include "../include/data_struct.h"
#include "../include/basic.h"
#include "../include/analysis.h"
#include "../include/trajectory.h"
#include "../include/perturb.h"

using namespace std;

typedef Vec<int> Vec_INT;
typedef Vec<double> Vec_DP;
typedef Mat<int> Mat_INT;
typedef Mat<double> Mat_DP;

// the seed of random generation method
int idum = -time(0);	

void test();

int main(int argc, char** argv) {
	extern int idum;
	printTime();
	
	string path_out = "./output/2016-11-5/";
	makeDir(path_out);
	
	// ***** initiate network and force field ***** //
	const int n_node = 34;
	string fname_mat = "./input/matrix.txt";
	string fname_node = "./input/gene_symbol.txt";
	Net net(fname_mat, fname_node, n_node);
	Force force(net);
	//net.t_mat.write("./output/2016-11-4/t_mat.txt");
	//net.w_mat.write("./output/2016-11-4/w_mat.txt");
	cout << "Network and Force Field initiate OK!" << endl;
	
	
	// ***** find all solutions test ***** //
	/*
	string fname_solutions = path_out + "fp_value.txt";
	Mat_DP solutions;
	find_all_solutions(solutions, force, 1000);
	solutions.print();
	solutions.write(fname_solutions);
	*/		
	

	// ***** trajectory generation test ***** //
	/*
	string fname_traj = path_out + "traj.txt";
	Mat_DP traj(0.0, 1000, 2);
	Vec_INT column_ids(0, 2);
	column_ids[1] = 1;
	Vec_DP x_start(0.5, n_node);
	generate_traj_det(traj, x_start, column_ids, force);
	//generate_traj_sto(traj, x_start, column_ids, force);
	//traj.print();
	traj.write(fname_traj);
	*/


	// ***** landscape generation test ***** //
	/*
	string fname_land = path_out + "freq_land.txt";
	string fname_solutions = path_out + "fp_value.txt";
	Mat_INT land(0, 200, 200);
	Mat_DP fx(0.0, 200, 200);
	Mat_DP fy(0.0, 200, 200);
	int xid = 0; 
	int yid = 1;
	int n_step = (int)1e8;
	//remove_edge(net, 9, 8);
	//force = Force(net);
	freq_land(land, fx, fy, xid, yid, n_step, force);
	land.write(fname_land);
	Mat_DP solutions;
	find_all_solutions(solutions, force, 1000);
	solutions.write(fname_solutions);
	*/
	
	// ***** random perturbing the weight matrix ***** //
	char sub_path_out[128];
	int n_sample = 1000;
	Mat_DP solutions;
	for (int i=0; i<n_sample; i++) {
		sprintf(sub_path_out,"sample%d/",i+1);
		string fname_solutions = path_out + sub_path_out + "fp_value.txt";
		string fname_w_mat = path_out + sub_path_out + "w_mat.txt";
		makeDir(fname_solutions);
		rand_w_mat(net, 0, 1);
		force = Force(net);
		find_all_solutions(solutions, force, 1000);
		net.w_mat.write(fname_w_mat);
		solutions.write(fname_solutions);
		cout << "ID: " << i+1 << '\t' << "No. FPs: " << solutions.nrows() << endl;
	}
	
	printTime();
  	return 0;
}



void test() {
	int dim = 10;
	Mat_DP mat(1.0, dim, dim);
	Mat_DP unique_mat;
	for (int i=0; i<dim; i++) {
		for (int j=0; j<dim; j++) {
			if (i % 2 == 0)
				mat[i][j] = 2;
		}
	}
	
	//mat.print();
	unique_mat_row(unique_mat, mat);
	unique_mat.print();
}

