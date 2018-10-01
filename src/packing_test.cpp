/*

	Methods implementation 
	for packing class

	BY Jack Treado

*/

#include "packing.h"
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <cmath>
#include <vector>

using namespace std;

#define DEBUG_OFF

const double PI = 3.1415926;

/* 
==================================

		PRINTER FUNCTIONS		 

================================== 
*/

void packing::print_config(){
	// throw error if file not opened
	if (!configobj.is_open()){
		cout << "ERROR: config file not opened!" << endl;		
		throw "ERROR: config file not opened!";
	}

	// local variables
	int w,p,i,d;
	w = 12;
	p = 6;	

	// print basic info
	configobj << setw(w) << N << endl;
	for (d=0; d<NDIM; d++)
		configobj << setw(w) << L[d];
	configobj << endl;
	configobj << setw(w) << setprecision(p) << phi << endl;

	// print header
	configobj << right << setw(w) << "id";
	configobj << right << setw(w) << "r";
	for (d=0; d<NDIM; d++)
		configobj << right << setw(w) << "x[" << d << "]";
	configobj << right << setw(w) << "m";
	configobj << endl;
	for (i=0; i<w*(NDIM+3); i++)
		configobj << "=";
	configobj << endl;

	// print data
	for (i=0; i<N; i++){
		configobj << setw(w) << i;
		configobj << setw(w) << setprecision(p) << r[i];
		for (d=0; d<NDIM; d++)
			configobj << setw(w) << setprecision(p) << x[i][d];
		configobj << setw(w) << setprecision(p) << m[i];
		configobj << endl;
	}	
}


void packing::print_stat(){
	// throw error if file not opened
	if (!statobj.is_open()){
		cout << "ERROR: config file not opened!" << endl;		
		throw "ERROR: config file not opened!";
	}

	// local variables
	int w,p,i,d;
	w = 10;
	p = 6;	

	// print stat info
	statobj << setw(w) << "isjammed: " << isjammed << endl;
	statobj << setw(w) << "N: " << N << endl;
	statobj << setw(w) << "L: ";
	for (d=0; d<NDIM; d++)
		statobj << setw(w) << L[d];
	statobj << endl;
	statobj << setw(w) << "phi: " << phi << endl;
	statobj << setw(w) << "U: " << U << endl;
	statobj << setw(w) << "K: " << K << endl;
	statobj << setw(w) << "csum: " << this->get_c_sum() << endl;
	statobj << setw(w) << "nr: " << nr << endl;
	statobj << setw(w) << "niso: " << DOF*(N-nr)-NDIM+1 << endl;
	statobj << setw(w) << "niso max: " << DOF*N-NDIM+1 << endl;
	statobj << setw(w) << "pc: ";
	this->print_pc(statobj,w);
	statobj << "contact matrix: " << endl;
	this->print_c_mat(statobj);
}

void packing::md_monitor(int t, int nr, double phiH, double phiL){
	int nbb,niso,i;
	nbb = N - nr;
	niso = DOF*nbb - NDIM + 1;

	// print info about system
	this->monitor_header(t);
	cout << "** phi = " << phi << endl;
	cout << "** phiH = " << phiH << endl;
	cout << "** phiL = " << phiL << endl;
	cout << "** U = " << U << endl;
	cout << "** K = " << K << endl;
	cout << "** E = " << K+U << endl;
	cout << "** mean v = " << this->get_mean_vel() << endl;
	cout << "** com = " << this->get_com() << endl;
	cout << "** c sum = " << this->get_c_sum() << endl;	
	cout << "** niso = " << niso << endl;
	cout << "** niso max = " << DOF*N-NDIM+1 << endl;
	cout << "** nbb = " << nbb << endl;
	cout << "** nr = " << nr << endl;
	if (N <= 300){
		cout << "** pc = ";
		for (i=0; i<N; i++){
			if (i % 10 == 0)
				cout << endl;
			cout << setw(4) << pc[i];		
		}
	}
}

void packing::monitor_scale(double dphi, double phiH, double phiL) {
	cout << "phiH = " << phiH << endl;
	cout << "phiL = " << phiL << endl;
	cout << "dphi = " << dphi << endl;
	cout << "total contacts = " << this->get_c_sum() << endl;
	cout << "niso max = " << DOF*N-NDIM+1 << endl;
	cout << "old phi = " << phi << endl;
	this->scale_sys(dphi);
	cout << "new phi = " << phi+dphi << endl;
	cout << endl;
}

void packing::print_xyz(){
	int i,d,dd,w;
	w = 20;	

	if (!xyzobj.is_open()){
		cout << "ERROR: xyz file obj not open!" << endl;
		throw "obj not open\n";
	}

	// print .xyz header
	xyzobj << N << endl;
	xyzobj << "Lattice=\"";
	for (d=0; d<NDIM; d++){
		for(dd=0; dd<NDIM; dd++){
			if (dd==d)
				xyzobj << L[d];
			else
				xyzobj << " 0.0 ";
		}
	}
	xyzobj << "\" ";
	xyzobj << '\t';
	xyzobj << "Properties=species:S:1:pos:R:" <<  NDIM << ":radius:R:1" << endl;

	for (i=0; i<N; i++){
		if (pc[i] == 0)
			xyzobj << setw(w) << "O";	
		else	
			xyzobj << setw(w) << "C";

		for (d=0; d<NDIM; d++){
			xyzobj << setw(w) << x[i][d];
		}
		xyzobj << setw(w) << r[i];
		xyzobj << endl;
	}
}

// print contact matrix
void packing::print_c_mat(){
	int i,j,ci,cj,w;

	// print width
	w = 6;

	for(i=0; i<N; i++){
		for (j=0; j<N; j++){
			// lower t
			if (j < i){
				ci = j*N - ((j+1)*(j+2))/2 + i;
				cout << setw(w) << c[ci];
			}
			// upper t
			else if (j > i){
				// get start of upper t row in linear array space
				ci = i*N - ((i+1)*(i+2))/2 + j;
				cout << setw(w) << c[ci];
			}
			// diag
			else
				cout << setw(w) << 0;			
		}
		cout << endl;
	}
}

void packing::print_c_mat(ofstream &obj){
	int i,j,ci,cj,w;

	// print width
	w = 6;

	for(i=0; i<N; i++){
		for (j=0; j<N; j++){
			// lower t
			if (j < i){
				ci = j*N - ((j+1)*(j+2))/2 + i;
				obj << setw(w) << c[ci];
			}
			// upper t
			else if (j > i){
				// get start of upper t row in linear array space
				ci = i*N - ((i+1)*(i+2))/2 + j;
				obj << setw(w) << c[ci];
			}
			// diag
			else
				obj << setw(w) << 0;			
		}
		obj << endl;
	}
}

void packing::print_c_data(ofstream &obj){
	int i,w,UTO;
	w = 6;

	// get upper triangle only number
	UTO = (N*(N-1))/2;

	// output entries only in the upper triangle
	for (i=0; i<UTO; i++)
		obj << setw(w) << c[i];
	obj << endl;

}

void packing::print_pc(ofstream &obj, int w){
	 int i;
	 for(i=0; i<N; i++)
	 	obj << setw(w) << pc[i];
	 obj << endl;
}


void packing::print_vars(){
	int i,d;

	cout << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&" << endl;
	cout << "Printing variables" << endl;
	cout << "U = " << U << endl;
	cout << "K = " << K << endl;
	cout << "phi = " << phi << endl;
	cout << "ep = " << ep << endl;

	cout << endl << "** arrays:" << endl << endl;
	cout << "x: " << endl;
	for (i=0; i<N; i++){
		for (d=0; d<NDIM; d++)
			cout << setw(20) << setprecision(6) << "x[" << i << "][" << d << "] = " << x[i][d];
		cout << endl;
	}
	cout << endl;
	cout << "v: " << endl;
	for (i=0; i<N; i++){
		for (d=0; d<NDIM; d++)
			cout << setw(20) << setprecision(6) << "v[" << i << "][" << d << "] = " << v[i][d];
		cout << endl;
	}
	cout << endl;
	cout << "F: " << endl;
	for (i=0; i<N; i++){
		for (d=0; d<NDIM; d++)
			cout << setw(20) << setprecision(6) << "F[" << i << "][" << d << "] = " << F[i][d];
		cout << endl;
	}
	cout << endl;
	cout << "aold: " << endl;
	for (i=0; i<N; i++){
		for (d=0; d<NDIM; d++)
			cout << setw(20) << setprecision(6) << "aold[" << i << "][" << d << "] = " << aold[i][d];
		cout << endl;
	}
	cout << endl;
	cout << "r: " << endl;
	for (i=0; i<N; i++)
		cout << setw(20) << setprecision(6) << "r[" << i << "] = " << r[i] << endl;
	cout << endl;
	cout << "m: " << endl;
	for (i=0; i<N; i++)
		cout << setw(20) << setprecision(6) << "m[" << i << "] = " << m[i] << endl;
	cout << endl;
	cout << "pc: " << endl;
	for (i=0; i<N; i++)
		cout << setw(20) << setprecision(6) << "pc[" << i << "] = " << pc[i] << endl;
	cout << endl;
}



// Neighbor list info to print
void packing::print_cell(){
	cout << "Printing cell info";
	cout << ": NCELLS = " << NCELLS;
	cout << "; NCL = " << NCL;
	cout << "; NCN = " << NCN;
	cout << "; nnupdate = " << nnupdate;
	cout << endl;

	int i,n,m;
	for (i=0; i<NCELLS; i++){
		cout << "cell " << i << ": ";
		n = celln[i];
		for (m=0; m<n; m++)
			cout << setw(5) << cell[i].at(m);
		cout << endl;
	}
}

void packing::print_celln(){
	int i;
	cout << "Printing celln:" << endl;
	for (i=0; i<NCELLS; i++)
		cout << setw(10) << celln[i];
	cout << endl << endl;
}

void packing::print_clabel(){
	int i;
	cout << "Printing clabel:" << endl;
	for (i=0; i<N; i++)
		cout << setw(10) << clabel[i];
	cout << endl << endl;
}

void packing::print_cell_pos(){
	int i,d;

	for (d=0; d<NDIM; d++)
		cout << setw(10) << g[d];
	cout << endl;

	for (i=0; i<NCELLS; i++){
		cout << "cell " << i << ": ";
		for (d=0; d<NDIM; d++)
			cout << setw(10) << setprecision(3) << cellpos[i][d];
		cout << endl;
	}
}

void packing::print_cell_neighbors(){
	int i,nn;
	cout << "Printing cell neighbors:" << endl;
	for (i=0; i<NCELLS; i++){
		for (nn=0; nn<NCN; nn++)
			cout << setw(4) << cellneighbors[i][nn];
		cout << endl;
	}			
	cout << endl << endl;
}

void packing::print_neighborlist(){
	cout << "Printing neighborlist info";
	cout << endl;

	int i,n,m;
	for (i=0; i<N; i++){
		cout << "NL " << i << ": ";
		n = neighborlist[i].size();
		for (m=0; m<n; m++)
			cout << setw(5) << neighborlist[i].at(m);
		cout << endl;
	}
}

void packing::print_nl_xyz(){
	int i,d,dd,w;
	w = 20;	

	if (!xyzobj.is_open()){
		cout << "ERROR: xyz file obj not open!" << endl;
		throw "obj not open\n";
	}

	// print .xyz header
	xyzobj << N << endl;
	xyzobj << "Lattice=\"";
	for (d=0; d<NDIM; d++){
		for(dd=0; dd<NDIM; dd++){
			if (dd==d)
				xyzobj << L[d];
			else
				xyzobj << " 0.0 ";
		}
	}
	xyzobj << "\" ";
	xyzobj << '\t';
	xyzobj << "Properties=species:S:1:pos:R:" <<  NDIM << ":radius:R:1" << endl;

	// print green if in nl of 1, red else
	int nl = neighborlist[0].size();
	int l = 0;
	int n_found = 0;
	for (i=0; i<N; i++){
		if (i == 0)
			xyzobj << setw(w) << 'C';
		else{			
			n_found = 0;
			for (int n=0; n<nl; n++){
				// get label 
				l = neighborlist[0].at(n);

				// if label == i, then i is in neighbor list of particle 0!
				if (l==i){					
					n_found = 1;
					break;
				}
				// else, not in neighbor list of 0
				else
					n_found = 0;
			}
			if (n_found == 1)
				xyzobj << setw(w) << 'H';
			else
				xyzobj << setw(w) << 'N';
		}

		for (d=0; d<NDIM; d++){
			xyzobj << setw(w) << x[i][d];
		}
		xyzobj << setw(w) << r[i];
		xyzobj << endl;
	}
}


void packing::print_all_nl_xyz(){
	int i,d,dd,w;
	w = 20;	

	if (!xyzobj.is_open()){
		cout << "ERROR: xyz file obj not open!" << endl;
		throw "obj not open\n";
	}

	// print green if in nl of 1, red else
	int K,M,l,k;
	vector<int> neigh;
	for (i=0; i<N; i++){
		// get nl size
		K = neighborlist[i].size();

		// push back neighbors
		for (l=0; l<K; l++)
			neigh.push_back(neighborlist[i][l]);

		// add upper t
		for (k=0; k<i; k++){
			M = neighborlist[k].size();
			for (l=0; l<M; l++){
				if (neighborlist[k][l] == i){
					neigh.push_back(k);
					break;
				}
			}
		}

		// update number of neighbors
		K = neigh.size();

		// print new .xyz header
		xyzobj << K+1 << endl;
		xyzobj << "Lattice=\"";
		for (d=0; d<NDIM; d++){
			for(dd=0; dd<NDIM; dd++){
				if (dd==d)
					xyzobj << L[d];
				else
					xyzobj << " 0.0 ";
			}
		}
		xyzobj << "\" ";
		xyzobj << '\t';
		xyzobj << "Properties=species:S:1:pos:R:" <<  NDIM << ":radius:R:1" << endl;	

		// print info for particle i
		xyzobj << setw(w) << 'H';
		for (d=0; d<NDIM; d++)
			xyzobj << setw(w) << x[i][d];
		xyzobj << setw(w) << r[i];
		xyzobj << endl;

		// print info for neighbors		
		for (l=0; l<K; l++){
			k = neigh[l];
			xyzobj << setw(w) << 'C';
			for (d=0; d<NDIM; d++)
				xyzobj << setw(w) << x[k][d];
			xyzobj << setw(w) << r[k];
			xyzobj << endl;

			// clear tmp neighbor list
			neigh.clear();
		}
	}
}


void packing::monitor_header(int t){
	cout << endl;
	cout << "=================" << endl;
	cout << "     t = " << t << endl;
	cout << "=================" << endl;
	cout << endl;
}



