#include <stdio.h>
#include <MtlOptSolver.h>
#include <Timer.h>
#include <mosek.h>
using namespace MTLOPT;

static void MSKAPI printstr(void *handle,MSKCONST char str[]){
  printf("%s",str);
}

void MtlOptSolver::optimizeLambdaDamping_mosec(){
  
  FUNC_TIMER();
  EwLambdaDamping->reset();

  fun_Ew=EwLambdaDamping->fun(dataModel->Lambda, dataModel->Damping);
  DEBUG_LOG( "inner-it-fun(lambda): "<<funValue() );

  VectorXd g;
  EwLambdaDamping->grad(g);
  const int NUMVAR = g.size();

  const double *c = &g[0];

  MSKboundkeye  bkx[NUMVAR];
  double blx[NUMVAR], bux[NUMVAR];
  for (int i = 0; i < NUMVAR/2; ++i){
    bkx[i] = MSK_BK_RA;
    blx[i] = 0.0f;
    bux[i] = 4.0f;
  }
  for (int i = NUMVAR/2; i < NUMVAR; ++i){
    bkx[i] = MSK_BK_LO;
    blx[i] = 0.0f;
    bux[i] = +MSK_INFINITY;
  }
  
  const int NUMQNZ = NUMVAR/2*3;
  MSKint32t     qsubi[NUMQNZ];
  MSKint32t     qsubj[NUMQNZ];
  double        qval[NUMQNZ];
  
  MSKint32t     i,j;
  double        xx[NUMVAR];

  MSKenv_t      env = NULL;
  MSKtask_t     task = NULL;
  MSKrescodee   r;
  
  /* Create the mosek environment. */
  r = MSK_makeenv(&env,NULL);

  if ( r==MSK_RES_OK ){ 

	/* Create the optimization task. */
	r = MSK_maketask(env,0,NUMVAR,&task);

	if ( r==MSK_RES_OK ){

	  r = MSK_linkfunctotaskstream(task,MSK_STREAM_LOG,NULL,printstr);

	  if ( r == MSK_RES_OK ) r = MSK_appendvars(task,NUMVAR);
  
	  /* Optionally add a constant term to the objective. */
	  VectorXd a(NUMVAR);
	  a.setZero();
	  const double f0 = EwLambdaDamping->fun(&a[0]);
	  if ( r ==MSK_RES_OK )	r = MSK_putcfix(task,f0);

	  for(j=0; j<NUMVAR && r == MSK_RES_OK; ++j){
		/* Set the linear term c_j in the objective.*/  
		if(r == MSK_RES_OK) r = MSK_putcj(task,j,c[j]);
		/* Set the bounds on variable j. blx[j] <= x_j <= bux[j] */
		if(r == MSK_RES_OK) r = MSK_putvarbound(task,j,bkx[j],blx[j],bux[j]);
	  }
  
	  if ( r==MSK_RES_OK ){
		/* The lower triangular part of the Q matrix in the objective is specified.*/
		EwLambdaDamping->hessianLower(qsubi, qsubj, qval, NUMQNZ);
		r = MSK_putqobj(task,NUMQNZ,qsubi,qsubj,qval);
	  }

	  if ( r==MSK_RES_OK ){

		MSKrescodee trmcode;
		/* Run optimizer */
		r = MSK_optimizetrm(task,&trmcode);

		/* Print a summary containing information about the solution for debugging purposes*/
		MSK_solutionsummary (task,MSK_STREAM_MSG);
        
		if ( r==MSK_RES_OK ){

		  MSKsolstae solsta;
		  int j;
          
		  MSK_getsolsta (task,MSK_SOL_ITR,&solsta);
          
		  switch(solsta){
		  case MSK_SOL_STA_OPTIMAL:   
		  case MSK_SOL_STA_NEAR_OPTIMAL:
			/* Request the interior solution. */
			MSK_getxx(task,MSK_SOL_ITR,xx);
			printf("Optimal primal solution\n");
			for(j=0; j<NUMVAR; ++j) printf("x[%d]: %e\n",j,xx[j]);
			break;

		  case MSK_SOL_STA_DUAL_INFEAS_CER:
		  case MSK_SOL_STA_PRIM_INFEAS_CER:
		  case MSK_SOL_STA_NEAR_DUAL_INFEAS_CER:
		  case MSK_SOL_STA_NEAR_PRIM_INFEAS_CER:  
			printf("Primal or dual infeasibility certificate found.\n");
			break;
              
		  case MSK_SOL_STA_UNKNOWN:
			printf("The status of the solution could not be determined.\n");
			break;

		  default:
			printf("Other solution status.");
			break;
		  }
		}else{
		  printf("Error while optimizing.\n");
		}
	  }
    
	  if (r != MSK_RES_OK){

		/* In case of an error print error code and description. */      
		char symname[MSK_MAX_STR_LEN];
		char desc[MSK_MAX_STR_LEN];
        
		printf("An error occurred while optimizing.\n");
		MSK_getcodedesc (r,symname,desc);
		printf("Error %s - '%s'\n",symname,desc);
	  }
	}
  }

  MSK_getxx(task,MSK_SOL_ITR,xx);
  for (int i = 0; i < NUMVAR/2; ++i){
	dataModel->Lambda[i] = xx[i];
	dataModel->Damping[i] = xx[i+NUMVAR/2];
  }
  
  MSK_deletetask(&task);
  MSK_deleteenv(&env);

  fun_Ew=EwLambdaDamping->fun(dataModel->Lambda, dataModel->Damping);
  DEBUG_LOG( "inner-it-fun(lambda): "<<funValue() );
}
