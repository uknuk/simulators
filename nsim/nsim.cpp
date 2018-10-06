// Main msim program 


#include "sim.h"
#include <stdio.h>
#include <memory.h>


int main(int argc, char *argv[])
{
		Sim* sim = new Sim(argc,argv);
		delete sim;
		printf("nsim executed \n");
		return 1;
}



