/*

	Jam N bidisperse disks in 2 dimensions

	USING RANDOM PROTOCOLS

	-- loop over P random protocols drawn from distribution, 
	-- and check packings again a master list of packing data.
	-- If iteration k is unique, add to list;
	-- if not, move to iteration k+1

	BY Jack Treado
	05/26/2019

*/

int main(int argc, char *argv[]){
	// get numerical values for input variables
	int N,seed,fskip;
	double dphi,phi0,alpha,phimin,tmp0;