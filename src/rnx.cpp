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


const std::vector<option::Descriptor> usage = {
     { 'h', "help",  option::ArgType::NONE, "Print usage message and exit" },
     { 'c', "ct", option::ArgType::STRING, "CT file" },
     { 'd', "data", option::ArgType::STRING, "Data directory with folding parameters" },
   };



int main(int argc, char* argv[]) {

 
  option::Parser parser(usage, argc, argv);
  std::string cmd = parser.get_command_string();
  std::cout << cmd << std::endl;
  const char * fname = "./testdata/u3.ct";
  

  const char *dir = "./dat";
  Datatable dattab(dir);
  RNAStructure rnastruct(fname);

  dattab.efn2(&rnastruct, 1);

  const char *newout = "testout.ct";
  rnastruct.ctout (newout);
  return 0;
}
