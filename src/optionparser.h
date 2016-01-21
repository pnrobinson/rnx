/*
 * An option parser for the rnx application. 
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



#ifndef OPTIONPARSER_H_
#define OPTIONPARSER_H_

#include <vector>
#include <sstream>
#include <string>


namespace option
{


  class Arg;

  class CheckArg{};
  /**
   * The four types of allowed argument: None, string, integer, and float value
   */
  enum class ArgType { NONE, STRING, INTEGER, FLOAT };
  
  struct Descriptor {
    /**
     * @brief  a short option character (without the leading @c - ).
     */
    const char shortopt;
    
    /**
     * @brief The long option name (without the leading @c -- ).
     *
     * If this Descriptor should not have a long option name, use the empty
     * string "". NULL is not permitted here!
     */
    const char* const longopt;
    
    /**
     * @brief One of NONE, STRING, INTEGER, or FLOAT.
     */
    const ArgType argtype;
    
    /**
     * @brief The usage text associated with the options in this Descriptor.
     *
     * You can use option::printUsage() to format your usage message based on
     * the @c help texts. You can use dummy Descriptors where
     * @ref shortopt and @ref longopt are both the empty string to add text to
     * the usage that is not related to a specific option.
     *
     * See option::printUsage() for special formatting characters you can use in
     * @c help to get a column layout.
     */
    const char* help;
  };



 

  class Arg {
    
  public :
    enum ArgType type;

  };


  /**
   * A single, completely parsed option.
   */
  class Option {

  private:

  public:
    Option(const char * args);
    Option (const Option &orig);
    void operator= (const Option &orig);
    /*ArgType type ()	const {
      return desc == 0 ? 0 : desc->argtype;
      }*/
 


  };




  class Parser {
    int op_count_; 
    int nonop_count_;
    std::string command_string_;
    std::string prog_name_;
    const char** nonop_args_; 
    bool err_;
    std::vector<Option> optionlist_;
    //Descriptor *desc_;
  public:
    /*
      Parser():  op_count_(0), nonop_count_(0), nonop_args_(0), err_(false) {
      }
    */
    Parser(std::vector<Descriptor> usage, int argc, const char** argv){
      input_command_string(argc,argv);
      parse_options(usage,argc,argv);
    }
    int optionsCount() { return op_count_; };
    int nonOptionsCount() {return nonop_count_; }
    //Returns the number of non-option arguments that remained at the end of the most recent parse() that actually encountered non-option arguments.
    //const char **nonOptions() {return nonop_args_; }
    //Returns a pointer to an array of non-option arguments (only valid if nonOptionsCount() >0 ).
    //const char * nonOption (int i);
    //Returns nonOptions()[i] (without checking if i is in range!).
    bool error() { return err_; }
    std::string get_command_string() const { return command_string_; }
    bool hasOption(char c) const;
    bool hasOption(const char *p) const;
    std::string get_program_name() const { return prog_name_; }

  private:
    void input_command_string(int argc, const char** argv) {
      std::stringstream os;
      if (argc==0) { return; }
      os << argv[0];
      for (unsigned int i=1;i<argc;++i) {
	os << " " << argv[i];
      }
      command_string_=os.str();
    }
    void parse_options(std::vector<Descriptor> usage,int argc, const char** argv) {
      if (argv==NULL || argv[0]==NULL)
	return; // Nothing to parse, should not happen.
      prog_name_ = argv[0];
      for (unsigned int i=1;i<argc;++i) {
	const char *p = argv[i];
	unsigned int len = strlen(p);
	if (*p=='-') {
	  if (len>2 && *(p+1)=='-') {
	    // long option
	    const char *q = p+2;
	    std::cout << p << ":Len=" << len << "  q=" << q << "\n";
	  }
	}
      }
    }
  };

  bool Parser::hasOption(char c) const {
    return false;
  }
  
}
// namespace option

#endif /* OPTIONPARSER_H_ */
