/*
 * rnx
 *
 * Copyright (C) 2015 Peter N Robinson
 *
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software, to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to permit
 * persons to whom the Software is furnished to do so, subject to the following
 * conditions:
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

/**
 * @file
 *
 * @brief This is the main driver file required for the rnx application.
 *    

 * @author Peter Robinson
 * @version 0.0.2 (Jan 18 2016)
 *
 */

#include <iostream>
#include "optionparser.h"
#include "Sequence.h"
#include "EnergyFunction2.h"
#include "RNAStructure.h"

enum  optionIndex { UNKNOWN, HELP, CT };

const option::Descriptor usage[] =
  {
    {UNKNOWN, 0,"" , ""    ,option::ArgType::NONE, "USAGE: rnx [options]\n\n"
     "Options:" },
    {HELP,    0,"" , "help",option::ArgType::NONE, "  --help  \tPrint usage and exit." },
    {UNKNOWN, 0,"" ,  ""   ,option::ArgType::NONE, "\nExamples:\n"
     "  rnx CTfile\n"
     "  rnx -unk --plus -ppp file1 file2\n" },
    //{0,0,0,0,0,0}
  };

int main(int argc, char* argv[]) {

  std::cout << "rnx " << std::endl;
  argc-=(argc>0); argv+=(argc>0); // skip program name argv[0] if present
  //option::Stats  stats(usage, argc, argv);
  //option::Option options[stats.options_max], buffer[stats.buffer_max];
  //option::Parser parse(usage, argc, argv, options, buffer);
  /*
  if (parse.error())
    return 1;
  
  if (options[HELP] || argc == 0) {
    option::printUsage(std::cout, usage);
    return 0;
  }

  if (options[CT]) {
    for (option::Option* opt = options[CT]; opt; opt = opt->next()) {
      const char *fname = opt->arg;
      std::cout << "FFF \"" << fname << "\"\n";
    }
    
  }
  */
  
  const char * fname = "./testdata/u3.ct";
  /*
  for (option::Option* opt = options[UNKNOWN]; opt; opt = opt->next())
     std::cout << "Unknown option: " << opt->name << "\n";
  
  for (int i = 0; i < parse.nonOptionsCount(); ++i)
    std::cout << "Non-option #" << i << ": " << parse.nonOption(i) << "\n";
  */

  const char *dir = "./dat";
  Datatable dattab(dir);
  RNAStructure rnastruct(fname);

  dattab.efn2(&rnastruct, 1);

  const char *newout = "testout.ct";
  rnastruct.ctout (newout);
  return 0;
}
