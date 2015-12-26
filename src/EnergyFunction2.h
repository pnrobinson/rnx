#ifndef ENERGY_FUNCTION_2_H
#define ENERGY_FUNCTION_2_H

/**
 * \defgroup Folding RNA Folding Algorithms
 */


/**
 * \class EnergyFunction2
 *
 * \ingroup Folding
 * 
 *
 * \brief calculate the folding free energy change of a structure. 
 * Based on e2f2 from the mfold package.
 *
 * Prototyping. Input files from mfold dat folder
 + //open the files using the C++ method for reading
  ml1.open(miscloop);
  lo1.open(loop2);#### DONE (loop.dat)
  st1.open(stackf);
  th1.open(tstackh);
  ti1.open(tstacki);
  tl1.open(tloop);
  da1.open(danglef);
  in1.open(int22);
  in2.open(int21);
  tri.open(triloop);
  co1.open(coax);
  co2.open(tstackcoax);
  co3.open(coaxstack);
  st2.open(tstack);
  tsm.open(tstackm);
  i11.open(int11);
  * List of input files
  * - loop.dat (lo1.open(loop2)). DESTABILIZING ENERGIES BY SIZE OF LOOP (corresponds to inter_, bulge_, and hairpin_)

 *
 * \note Still prototyping.
 *
 * \author Peter Robinson
 *
 * \version  0.0.2
 *
 * \date 25 December 2015
 *
 * Contact: peter.robinson@charite.de
 *
 * Created on: 26 December 2015
 *
 *
 */

#include <string>
#if !defined(DEFINES_H)
#define DEFINES_H
#define maxfil 100    //maximum length of file names
//#define infinity   //an arbitrary value given to infinity
#define maxtloop 100 //maximum tetraloops allowed (info read from tloop)
#define maxstructures 1010 //maximum number of structures in ct file
#define maxbases 10000   //maximum number of bases in a structure
#define ctheaderlength 125 //max length of string containing info on sequence
#define ga_bonus -10 //bonus value for "almost coaxial stacking" in efn2
#define amax 400 //this is a maximum line length for void linout (below)
#define col 80  //this is the number of columns in an output file
#define numlen 8  //maximum digits in a number
#define maxforce 600 //maximum number of bases that can be forced single
#define maxgu 5 //maximum number of u's in gu pair
#endif





class Datatable {
  static const int s_infinity = 9999999;
  /**  The path to the directory containing the thermodynamic datafiles. */
  char * dirpath_;
  int poppen [5],maxpen,eparam[11],dangle[6][6][6][3];
  /** DESTABILIZING ENERGIES BY SIZE OF LOOP (INTERNAL). Loop sizes from 1-30.
   * Note that index 0 is not used. */
  int inter_[31];
  /** DESTABILIZING ENERGIES BY SIZE OF LOOP (BULGE). Loop sizes from 1-30.
   * Note that index 0 is not used. */
  int bulge_[31];
  /** DESTABILIZING ENERGIES BY SIZE OF LOOP (HAIRPIN). Loop sizes from 1-30.
   * Note that index 0 is not used. */
  int hairpin_[31];
  int stack[6][6][6][6],tstkh[6][6][6][6],tstki[6][6][6][6],
    tloop[maxtloop+1][2],numoftloops,iloop22[6][6][6][6][6][6][6][6],
    iloop21[6][6][6][6][6][6][6],iloop11[6][6][6][6][6][6],
    coax[6][6][6][6],tstackcoax[6][6][6][6],coaxstack[6][6][6][6],
    tstack[6][6][6][6],tstkm[6][6][6][6],auend,gubonus,cint,cslope,c3,
    efn2a,efn2b,efn2c,triloop[maxtloop+1][2],numoftriloops,init,gail;
  float prelog;
  /**
   * @param directory_path The path to the directory containing the thermodynamic datafiles.
   */
public:
  Datatable(const char* directory_path);
  ~Datatable();

  const char* get_data_dir() const;

  int get_destabilizing_energy_internal_loop(unsigned int loop_size);
  int get_destabilizing_energy_bulge_loop(unsigned int loop_size);
  int get_destabilizing_energy_hairpin_loop(unsigned int loop_size);

private:
  void input_data();
  void input_loop_dat(std::string &path);
};





#endif
/* eof */
