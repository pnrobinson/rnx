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
 * List of input files
 * - loop.dat (lo1.open(loop2)). DESTABILIZING ENERGIES BY SIZE OF LOOP (corresponds to inter_, bulge_, and hairpin_)
 * - stack.dat (st1.open(stackf)). stack energies
 * - tstackh.dat th1.open(tstackh); STACKING ENERGIES : TERMINAL MISMATCHES AND BASE-PAIRS
 * - tstacki.dat  ti1.open(tstacki); STACKING ENERGIES : TERMINAL MISMATCHES AND BASE-PAIRS (?what is distinction to tstackh?)
 * - tloop.dat  tl1.open(tloop)
 * - miscloop.data  ml1.open(miscloop); Miscellaneous free energy rules Extrapolation for large loops based on polymer theory, internal, bulge or hairpin loops > 30: dS(T)=dS(30)+param*ln(n/30) 
 * - dangle.dat  da1.open(danglef);
 * - int22.dat  in1.open(int22); Data tables for symetric interior loops of size 4 
 * - int21.dat
 * - coaxial.dat co1.open(coax);
 * - triloop.dat  tri.open(triloop);
 * - tstackcoax
 * - coaxstack.dat
 * - tstack.dat st2.open(tstack);
 * - tstackm.dat   tsm.open(tstackm);
 * - i11.open(int11);
 *
 * The program first reads in all of these files in the indicated order.
 *
 * \note Still prototyping.
 *
 * \author Peter Robinson
 *
 * \version  0.0.2
 *
 * \date 28 December 2015
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
//#define maxtloop 100 //maximum tetraloops allowed (info read from tloop)
#define maxstructures 1010 //maximum number of structures in ct file
#define maxbases 10000   //maximum number of bases in a structure
//#define ctheaderlength 125 //max length of string containing info on sequence
#define ga_bonus -10 //bonus value for "almost coaxial stacking" in efn2
#define amax 400 //this is a maximum line length for void linout (below)
#define col 80  //this is the number of columns in an output file
#define numlen 8  //maximum digits in a number
#define maxforce 600 //maximum number of bases that can be forced single
#define maxgu 5 //maximum number of u's in gu pair
#endif


class RNAStructure;


class Datatable {
  /** A large number (infinity for the minimisation operations). */
  static const int s_infinity = 9999999;
  /** maximum tetraloops allowed (info read from tloop). */
  static const int s_maxtloop = 100;
  /**  The path to the directory containing the thermodynamic datafiles. */
  char * dirpath_;
  /** the f(m) array (see Ninio for details)  */
  int poppen_[5];
  /** asymmetric internal loops: the ninio equation the maximum correction (miscloop.dat) */
  int maxpen_;
  /** stores various values from miscloop.dat */
  int eparam_[11];
  /** stores values from dangle.dat */
  int dangle_[6][6][6][3];
  /** DESTABILIZING ENERGIES BY SIZE OF LOOP (INTERNAL). Loop sizes from 1-30.
   * Note that index 0 is not used. */
  int inter_[31];
  /** DESTABILIZING ENERGIES BY SIZE OF LOOP (BULGE). Loop sizes from 1-30.
   * Note that index 0 is not used. */
  int bulge_[31];
  /** DESTABILIZING ENERGIES BY SIZE OF LOOP (HAIRPIN). Loop sizes from 1-30.
   * Note that index 0 is not used. */
  int hairpin_[31];
  /** stacking energy */
  int stack_[6][6][6][6];
  /** STACKING ENERGIES : TERMINAL MISMATCHES AND BASE-PAIRS */
  int tstkh_[6][6][6][6];
  /** STACKING ENERGIES : TERMINAL MISMATCHES AND BASE-PAIRS (unclear what distinction is to tstkh?) */
  int tstki_[6][6][6][6];
  /** tetraloops. Storing data from tloop.dat, e.g., the first is GGGGAC -3.0 (then tloop_[1][1]=-300).
   * the left column is a numerical value calculated from the sequence of the loop using the
   * tonumi function, e.g., for tloop_[1][1] =7343 for sequence GGGGAC. */
  int tloop_[s_maxtloop+1][2];
  /** Number of tetraloops */
  int numoftloops_;
  /** Data tables for symetric interior loops of size 4 (int22.dat). Key iloop22[a][b][c][d][j][l][k][m] =
   * a j l b
   * c k m d*/
  int iloop22_[6][6][6][6][6][6][6][6];
  /** Data tables for asymmetric interior loops of size 3  */
  int iloop21_[6][6][6][6][6][6][6];
  /** The lookup table for single mismatches. */
  int iloop11_[6][6][6][6][6][6];
  /** The first of the tables for coaxial stacking with an intervening mismatch. This is the stack with the open backbone.*/
  int coax_[6][6][6][6];
  /** The second of the tables for coaxial stacking with an intervening mismatch. This is the stack with the continuous backbone.*/
  int tstackcoax_[6][6][6][6];
  /** The first of the tables for coaxial stacking with an intervening mismatch. This is the stack with the open backbone.*/
  int coaxstack_[6][6][6][6];
  /** For terminal mismatch stacking in exterior loops, i.e. loops that contain the ends of the sequence. (tstack.dat). */
  int tstack_[6][6][6][6];
  /** Used for terminal stacking in a multibranch loop. */
  int tstkm_[6][6][6][6];
  /** terminal AU penalty  (miscloop.dat) */
  int auend_;
  /** bonus for GGG hairpin  */
  int gubonus_;
  /** c hairpin intercept  */
  int cint_;
  /** c hairpin slope (miscloop.dat)*/
  int cslope_;
  /** c hairpin of 3 */
  int c3_;
  /** constant multi-loop penalty for efn2 */
  int efn2a_;
  /** efn2 multibranched loops free base penalty. */
  int efn2b_;
  /** efn2 efn2 multibranched loops helix penalty */ 
  int efn2c_;
  /** from triloop.dat (empty in mfold but two entries in Matthews 2004) 
   * A lookup table for hairpin loops of three. It is applied analogously to the tetraloop table in Tloop.dat.*/
  int triloop_[s_maxtloop+1][2];
  /** Number of triloops entered from triloop.dat */
  int numoftriloops_;
  /** Intermolecular initiation free energy (miscloop.dat) */
  int init_;
  /** GAIL Rule (Grossly Asymmetric Interior Loop Rule) (on/off <-> 1/0)  */
  int gail_;
  /** Extrapolation for large loops based on polymer theory 
   * internal, bulge or hairpin loops > 30: dS(T)=dS(30)+param*ln(n/30) (from miscloop.dat) */
  float prelog_;
  /**
   * @param directory_path The path to the directory containing the thermodynamic datafiles.
   */
public:
  Datatable(const char* directory_path);
  ~Datatable();

  const char* get_data_dir() const;

  int get_destabilizing_energy_internal_loop(unsigned int loop_size) const;
  int get_destabilizing_energy_bulge_loop(unsigned int loop_size) const;
  int get_destabilizing_energy_hairpin_loop(unsigned int loop_size) const;
  int get_stack_energy(int w, int x, int y, int z) const;
  int get_tstackh_energy(int w, int x, int y, int z) const;
  int get_tstacki_energy(int w, int x, int y, int z) const;
  int get_tetraloop_energy(const char *seq) const;
  float get_prelog() const;
  int get_maxpen() const;
  int get_poppen(unsigned int i) const;
  int get_constant_multiloop_penalty() const;
  int get_constant_efn2_multiloop_penalty() const;
  int get_terminal_AU_penalty() const;
  int get_GU_bonus() const;
  int get_c_hairpin_intercept() const;
  int get_c_hairpin_slope() const;
  int get_c_hairpin_of_3() const;
  int get_intermolecular_initiation_free_energy() const;
  int get_GAIL() const;
  int get_dangle_energy(int i, int j, int k, int l) const;
  int get_iloop22(int a, int b, int c, int d, int j, int k, int l, int m) const;
  int get_coaxial_energy(int i, int j, int k, int l) const;
  int get_tstack_coaxial_energy(int i, int j, int k, int l) const;

  void efn2(RNAStructure *ct, int structnum);
  int erg1(int i, int j, int ip, int jp, RNAStructure *ct);
  int erg3(int i, int j, RNAStructure *ct, int dbl);
  
private:
  void input_data();
  void input_loop_dat(const std::string &path);
  void input_stack_dat(const std::string &path);
  void input_tstackh_dat(const std::string &path);
  void input_tstacki_dat(const std::string &path);
  void input_tloop_dat(const std::string &path);
  void input_miscloop_dat(const std::string &path);
  void input_dangle_dat(const std::string &path);
  void input_int22_dat(const std::string &path);
  void input_int21_dat(const std::string &path);
  void input_coaxial_dat(const std::string &path);
  void input_triloop_dat(const std::string &path);
  void input_tstackcoax_dat(const std::string &path);
  void input_coaxstack_dat(const std::string &path);
  void input_tstack_dat(const std::string &path);
  void input_tstackm_dat(const std::string &path);
  void input_int11_dat(const std::string &path);
  
//this function calculates whether a terminal pair i, j requires the end penalty
  int penalty2(int i, int j) const; 



};





#endif
/* eof */
