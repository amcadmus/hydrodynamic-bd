#include <stdio.h>
#include "MDSystem_interface.h"
#include "common.h"
#include "BoxGeometry.h"
#include "MDSystem.h"
#include "RandomGenerator.h"
#include "Auxiliary.h"
#include "NeighborList_interface.h"
#include"Statistic.h"
#include "Integrator_interface.h"
#include "InteractionEngine_interface.h"
#include "Reshuffle_interface.h"
#include "Displacement_interface.h"

#include "Topology.h"
#include "SystemBondedInteraction.h"

#include "BondInteraction.h"
#include "NonBondedInteraction.h"


// #define NThreadsPerBlockCell	32
// #define NThreadsPerBlockAtom	4

#define NThreadsPerBlockCell	160
#define NThreadsPerBlockAtom	96

void posiRenormalize (int * posi,
		      const int size)
{
  while ((*posi) < 0) {
    *posi += size;
  }
  while ((*posi) >= size) {
    *posi -= size;
  }
}

double calDeltaF (const MDSystem & sys,
		  const double & kk,
		  const char * liquidName,
		  double * corrDeltaF,
		  const int corrSize,
		  int * corrDeltaFPosi,
		  int * corrDeltaFNvalid)
{
  double sumFup = 0; 
  double sumFdown = 0;
  double boxc = sys.box.size.z * 0.5;
  
  for (int ii = 0; ii < sys.hdata.numAtom; ++ii){
    if (strcmp(&(sys.hdata.atomName[ii*StringSize]), liquidName) == 0) {
      if ( sys.hdata.coord[ii].z > sys.box.size.z * 0.5 ){
	sumFup += sys.hdata.recordForcx[ii] * sin( kk * (sys.hdata.coord[ii].z - boxc) );
      }
      else {
	sumFdown += sys.hdata.recordForcx[ii] * sin( kk * (sys.hdata.coord[ii].z - boxc) );
      }
    }
  }

//  printf ("%f   ", sumFup - sumFdown);
  corrDeltaF[*corrDeltaFPosi] = (sumFup - sumFdown) * 0.5;
  (*corrDeltaFPosi) ++;
  posiRenormalize (corrDeltaFPosi, corrSize);
  
  if ((*corrDeltaFNvalid) < corrSize) {
    (*corrDeltaFNvalid) ++;
  }

  return (sumFup - sumFdown) * 0.5;
}

void depositData (const double * corrDeltaF,
		  const int corrSize,
		  const int corrDeltaFPosi,
		  const int corrDeltaFNvalid,
		  double * corrSumData,
		  int * corrSumDataCount,
		  double * corrData,
		  int * corrDataCount)
{
  int start = corrDeltaFPosi - 1;
  int end = corrDeltaFPosi - corrDeltaFNvalid - 1; 
  posiRenormalize (&start, corrSize);
  posiRenormalize (&end, corrSize);
  
  int ii = start;
  corrSumData[0] += corrDeltaF[start];
//  printf (" %f\n", corrDeltaF[start]);
  corrSumDataCount[0] ++;
  corrData[0] += corrDeltaF[start] * corrDeltaF[start];
  corrDataCount[0] ++;
  ii --;
  posiRenormalize (&ii, corrSize);
  int count = 1;
  
  while (ii != end){
    corrData[count] += corrDeltaF[ii] * corrDeltaF[start];
    corrDataCount[count] ++;
    --ii;
    ++count;
    posiRenormalize(&ii, corrSize);
  }
}


int main(int argc, char * argv[])
{
  IndexType nstep = 100000;
  IndexType confFeq = 250;
  IndexType thermoFeq = 5000;
  ScalorType rcut = 2.5;
  ScalorType nlistExten = 0.5;
  ScalorType refT = 2.80;
  ScalorType tauT = 0.1;
  char * filename;
  InteractionType recordType = 1;
  IndexType numType0 = 0;
  IndexType numType1 = 0;
  IndexType numType2 = 0;
  IndexType num_w = 0;
  IndexType num_f = 0;
  double kk = 0.224214;
  ScalorType lattice_k = 900.;
  
  if (argc != 4){
    printf ("Usage:\n%s conf.gro nstep device\n", argv[0]);
    return 1;
  }
  if (argc != 1){
    nstep = atoi(argv[2]);
    filename = argv[1];
  }
  printf ("# setting device to %d\n", atoi(argv[3]));
  cudaSetDevice (atoi(argv[3]));
  checkCUDAError ("set device");

  IndexType corrFeq = 1;
  IndexType corrStart = nstep / 10;
  IndexType corrNstep = 10000 + 1;
  double * corrDeltaF = (double *) malloc (sizeof(double) * corrNstep);
  double * corrData = (double *) malloc (sizeof(double) * corrNstep);
  int * corrDataCount = (int *) malloc (sizeof(int) * corrNstep);
  double corrSumData = 0.;
  int corrSumDataCount = 0;
  int corrDeltaFPosi = 0;
  int corrDeltaFNvalid = 0;
  for (int ii = 0; ii < corrNstep; ++ii){
    corrDeltaF[ii] = 0.;
    corrData[ii] = 0.;
    corrDataCount[ii] = 0;
  }
  double sumDF = 0.;
  double sumD2F = 0.;
  int sumD2FCount = 0;

  MDSystem sys;
  sys.initConfig(filename);
  // double com = 0;
  // for (int ii = 0; ii < sys.hdata.numAtom; ++ii){
  //   com += sys.hdata.coord[ii].z;
  // }
  // com /= double (sys.hdata.numAtom);
  // com = sys.box.size.z * 0.5 - com;
  // for (int ii = 0; ii < sys.hdata.numAtom; ++ii){
  //   sys.hdata.coord[ii].z += com;
  // }  
  
  unsigned position = 0;
  for (unsigned ii = position; ii < sys.hdata.numAtom ; ++ii){
    if (strcmp(&(sys.hdata.atomName[ii*StringSize]), "tp01") == 0) {
      continue;
    }
    else {
      numType0 = ii;
      break;
    }
  }
  position = numType0;
  for (unsigned ii = position; ii < sys.hdata.numAtom ; ++ii){
    if (strcmp(&(sys.hdata.atomName[ii*StringSize]), "tp02") == 0) {
      continue;
    }
    else {
      numType1 = ii;
      break;
    }
  }
  numType1 -= position;
  numType2 = sys.hdata.numAtom - numType1 - numType0;
  num_w = numType0 + numType1;
  num_f = numType2;
  printf ("# num tp01: %d tp02: %d ty03: %d   num wall: %d num fluid: %d\n",
	  numType0, numType1, numType2,
	  num_w, num_f);

  double minz = 1e10;
  for (int ii = 0; ii < numType0; ++ii){
    if (sys.hdata.coord[ii].z < minz){
      minz = sys.hdata.coord[ii].z;
    }
  }
  double maxz = -1;
  for (int ii = numType0; ii < numType0 + numType1; ++ii) {
    if (sys.hdata.coord[ii].z > maxz){
      maxz = sys.hdata.coord[ii].z;
    }
  }
  double shift = 0.5 * sys.box.size.z - 0.5 * (maxz + minz);
  // for (int ii = 0; ii < sys.hdata.numAtom; ++ii){
  //   sys.hdata.coord[ii].z += shift;
  // }
  printf ("# minz is %f   maxz is %f   shift the system by %f\n", minz, maxz, shift);

  TypeType freezType ;
  Topology::System sysTop;
  Topology::Molecule mol_w;
  mol_w.pushAtom (Topology::Atom (1.0, 0.0, 0));
  freezType = 100;
  Topology::Molecule mol_f;
  mol_f.pushAtom (Topology::Atom (1.0, 0.0, 1));

  sysTop.addMolecules (mol_w, num_w);
  sysTop.addMolecules (mol_f, num_f);
  
  LennardJones6_12Parameter ljparam_ww;
  ljparam_ww.reinit (1.00f, 1.00f, 1.00f, 0.f, rcut);
  LennardJones6_12Parameter ljparam_ff;
  ljparam_ff.reinit (1.00f, 1.00f, 1.00f, 0.f, rcut);
  LennardJones6_12Parameter ljparam_fw;
  ljparam_fw.reinit (1.16f, 1.04f, 0.70f, 0.f, rcut);

//  sysTop.addNonBondedInteraction (Topology::NonBondedInteraction(0, 0, ljparam_ww));
  sysTop.addNonBondedInteraction (Topology::NonBondedInteraction(0, 1, ljparam_fw));
  sysTop.addNonBondedInteraction (Topology::NonBondedInteraction(1, 1, ljparam_ff));

  sys.initTopology (sysTop);
  sys.initDeviceData ();
  
  SystemNonBondedInteraction sysNbInter;
  sysNbInter.reinit (sysTop);
  ScalorType energyCorr = sysNbInter.energyCorrection ();
  ScalorType pressureCorr = sysNbInter.pressureCorrection ();
  sysNbInter.printTable ();
  
  ScalorType maxrcut = sysNbInter.maxRcut();
  ScalorType rlist = maxrcut + nlistExten;
  CellList clist (sys, rlist, NThreadsPerBlockCell, NThreadsPerBlockAtom);
  NeighborList nlist (sysNbInter, sys, rlist, NThreadsPerBlockAtom, 2.f);
  sys.normalizeDeviceData ();
  clist.rebuild (sys, NULL);
  nlist.rebuild (sys, clist, NULL);
  Displacement_max disp (sys, NThreadsPerBlockAtom);
  disp.recordCoord (sys);
  
  MDStatistic st(sys);
  TranslationalFreedomRemover tfremover (sys, NThreadsPerBlockAtom);
  InteractionEngine inter (sys, NThreadsPerBlockAtom);
  inter.registNonBondedInteraction (sysNbInter);
  // inter.applyNonBondedInteraction (sys, nlist, recordType, st, NULL, &timer);
  // inter.applyLatticeInteraction (sys, lattice_k, 0, NULL);
  
  MDTimer timer;
  unsigned i;
  ScalorType dt = 0.002;
  ScalorType seed = 1;
  RandomGenerator_MT19937::init_genrand (seed);

  VelocityVerlet inte_vv (sys, NThreadsPerBlockAtom);
  VelocityRescale inte_vr (sys, NThreadsPerBlockAtom, refT, 0.1);
  NoseHoover_Chains2 nhc;
  // nhc.reinit (sys, num_w, NThreadsPerBlockAtom, refT, tauT);
  nhc.reinit (sys, 0, NThreadsPerBlockAtom, refT, tauT);

  Reshuffle resh (sys);
  
  timer.tic(mdTimeTotal);
  if (resh.calIndexTable (clist, &timer)){
    sys.reshuffle   (resh.indexTable, sys.hdata.numAtom, &timer);
    clist.reshuffle (resh.indexTable, sys.hdata.numAtom, &timer);  
    nlist.reshuffle (resh.indexTable, sys.hdata.numAtom, &timer);  
    disp.reshuffle  (resh.indexTable, sys.hdata.numAtom, &timer);  
  }
  
  printf ("# prepare ok, start to run\n");
  sys.recoverDeviceData (&timer);
  sys.updateHostFromRecovered (&timer);
  sys.writeHostDataGro ("confstart.gro", 0, 0.f, &timer);
  printf ("# prepare ok, start to run\n");
  printf ("#*     1     2           3         4            5       6                7        8   9\n");
  printf ("#* nstep  time  nonBondedE  kineticE  temperature  totalE  NHC_Hamiltonian pressure box\n");

  try{
    sys.initWriteXtc ("traj.xtc");
    sys.initWriteTrr ("traj.trr");
    sys.recoverDeviceData (&timer);
    sys.updateHostFromRecovered (&timer);
    sys.writeHostDataXtc (0, 0*dt, &timer);
    sys.writeHostDataTrr (0, 0*dt, &timer);
    for (i = 0; i < nstep; ++i){
      if (i%10 == 0){
//	tfremover.remove (sys, &timer);
      }
      
      nhc.operator_L (0.5 * dt, sys, &timer);
      inte_vv.step1 (sys, freezType, dt, &timer);

      st.clearDevice();
      ScalorType maxdr = disp.calMaxDisplacemant (sys, &timer);
      if (maxdr > nlistExten * 0.5){
	// printf ("# Rebuild at step %09i ... ", i+1);
	// fflush(stdout);
	// rebuild
	sys.normalizeDeviceData (&timer);
	disp.recordCoord (sys);
	clist.rebuild (sys, &timer);
	nlist.rebuild (sys, clist, &timer);
	// printf ("done\n");
	// fflush(stdout);
      }
      inter.clearInteraction (sys);
      inter.applyNonBondedInteraction (sys, nlist, recordType, st, NULL, &timer);
      inter.applyLatticeInteraction (sys, lattice_k, 0, NULL);

      inte_vv.step2 (sys, freezType, dt, &timer);
      if ((i+1) % thermoFeq == 0){	
	nhc.operator_L (0.5 * dt, sys, st, &timer);
      }
      else {
	nhc.operator_L (0.5 * dt, sys, &timer);	
      }      

      if ((i+1) % thermoFeq == 0){
	st.updateHost ();
	printf ("%09d %07e %.7e %.7e %.7e %.7e %.7e %.7e %.7e %.2e %.2e %.2e %.2e\n",
		(i+1),  
		(i+1) * dt, 
		st.nonBondedEnergy(),
		st.kineticEnergy(),
		st.kineticEnergy() * 2. / 3. / (double (num_f) - 3.),
		st.nonBondedEnergy() +
		st.kineticEnergy(),
		st.nonBondedEnergy() +
		st.kineticEnergy() +
		nhc.HamiltonianContribution (),
		st.pressure(sys.box),
		sys.box.size.x,
		nhc.xi1,
		nhc.vxi1,
		nhc.xi2,
		nhc.vxi1
	    );
	fflush(stdout);
      }

      if (i+1 >= corrStart && (i+1) % confFeq == 0){
      	sys.recoverDeviceData (&timer);
      	sys.updateHostFromRecovered (&timer);
      	sys.writeHostDataXtc (i+1, (i+1)*dt, &timer);
      	sys.writeHostDataTrr (i+1, (i+1)*dt, &timer);
      }
      // if (i+1 >= corrStart && (i+1) % corrFeq == 0){
      // 	sys.recoverDeviceData (&timer);
      // 	sys.updateHostFromRecovered (&timer);
      // 	double DF = calDeltaF (sys, kk, "tp03",
      // 			       corrDeltaF, corrNstep,
      // 			       &corrDeltaFPosi, &corrDeltaFNvalid);
      // 	sumDF += DF;
      // 	sumD2F += DF * DF;
      // 	sumD2FCount++;
      // 	depositData (corrDeltaF, corrNstep, corrDeltaFPosi, corrDeltaFNvalid,
      // 		     &corrSumData, &corrSumDataCount,
      // 		     corrData, corrDataCount);
      // }

      if ((i+1) % 100 == 0){
	if (resh.calIndexTable (clist, &timer)){
	  sys.reshuffle   (resh.indexTable, sys.hdata.numAtom, &timer);
	  clist.reshuffle (resh.indexTable, sys.hdata.numAtom, &timer);  
	  nlist.reshuffle (resh.indexTable, sys.hdata.numAtom, &timer);  
	  disp.reshuffle  (resh.indexTable, sys.hdata.numAtom, &timer);  
	}
      }
    }
    sys.endWriteXtc();
    sys.endWriteTrr();
    sys.recoverDeviceData (&timer);
    sys.updateHostFromRecovered (&timer);
    sys.writeHostDataGro ("confout.gro", nstep, nstep*dt, &timer);
    timer.toc(mdTimeTotal);
    timer.printRecord (stderr);
  }
  catch (MDExcptCuda & e){
    // resh.recoverMDDataToHost (sys, &timer);
    // sys.writeHostDataXtc (i+1, (i+1)*dt, &timer);
    timer.toc(mdTimeTotal);
    timer.printRecord (stderr);
    return 1;
  }
  catch (MDException &e){
    fprintf (stderr, "%s\n", e.what());
    return 1;
  }


  FILE * fp;
  if ((fp = fopen("corr.out","w")) == NULL) {
    printf ("cannot open file %s\n", "corr.out");
  }
  for (int ii = 0; ii < corrNstep; ++ii){
    if (corrDataCount[ii] != 0){
      double tmpa = corrData[ii] / double(corrDataCount[ii]);
      double tmpb = corrSumData  / double(corrSumDataCount);
      fprintf (fp, "%f %e   %e %e\n", ii * dt * corrFeq, tmpa - tmpb * tmpb, tmpa, tmpb);
    }
    else {
      fprintf (fp, "%f %e\n", ii * dt * corrFeq, 0.);
    }
  }
  free (corrDeltaF);
  free (corrData);
  free (corrDataCount);

  printf ("avg sumD2F: %f, avg sumDF: %f\n",
	  sumD2F / double(sumD2FCount),
	  sumDF / double (sumD2FCount)
      );
  
  return 0;
}

  
