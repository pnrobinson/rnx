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
 * 
 * @author Peter Robinson
 * @version 0.1 (22 Jan 2016)
 */



#ifndef OPTIONPARSER_H_
#define OPTIONPARSER_H_

#include <vector>
#include <sstream>
#include <string>
#include <cstring>
//#include <unistd.h>

namespace option
{
  /**
   * The four types of allowed argument: None, string, integer, and float value
   */
  enum class ArgType { NONE, STRING, INTEGER, FLOAT };
  /**
   * A simple struct that is used to define the 
   * expected command line arguments.
   */
  struct Descriptor {
    /**
     * A short option character (without the leading @c - ).
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

  /**
   * A single, completely parsed option.
   */
  class Option {

  private:
    /** The short option for this argument */
    char shortopt_;
    /**
     * @brief The long option name (without the "--");
     */
    std::string longopt_;
    /**
     * @brief One of NONE, STRING, INTEGER, or FLOAT.
     */
    ArgType argtype_;
    /**
     * The value of the option (or NULL if there is no value )
     */
    std::string value_;
  public:
    /**
     * Initialize the option based on the Description that matches the 
     * flag and its argument (if any). Note that we strip the "--" from
     * the flag if needed.
     * @param desc The Description object that matched the flag
     * @param arg The corresponding argument (except for ArgType::NONE).
     */
    Option(const Descriptor & desc,
	   const char *  arg 
	   ):shortopt_(desc.shortopt), value_(arg), argtype_(desc.argtype) {
      size_t len  = strlen(desc.longopt);
      size_t i=0;
      if (0<len && desc.longopt[0]=='-') i=1;
      if (1<len && desc.longopt[1]=='-') i=2;
      longopt_=std::string(desc.longopt+i,desc.longopt+len);
    }
    Option (const Option &orig) {
      shortopt_ = orig.shortopt_;
      longopt_=orig.longopt_;
      argtype_ = orig.argtype_;
      value_=orig.value_;
    }
    ~Option() {
      /* no-op */
    }
    inline char get_shortopt() const { return shortopt_; }
    inline std::string get_longopt() const { return longopt_; }
    inline std::string get_value() const { return value_; }
    void operator= (const Option &orig);
  };



  /**
   * This class coordinates the parsing of the command line arguments.
   */
  class Parser {
    /** A copy of the entire command used to run the program. */
    std::string command_string_;
    /** A copy of the name of the program. */
    std::string prog_name_;
    /** A flag to indicate that some error occured. Usually, client code should check this and die with
	an error message if there has ben an error. */
    bool err_;
    /** List of all of the successfully parsed Options */
    std::vector<Option> option_list_;
    /** List of all arguments that do not have an option flag (i.e. that are not proceeded by -c or --foo)*/
    std::vector<std::string> nonoption_list_;
   
  public:
    /**
     * Note that the signature const char*const* is needed to allow us to 
     * pass a non-const char** as an argument.
     */
    Parser(std::vector<Descriptor> usage, int argc, const char* const * argv){
      input_command_string(argc,argv);
      parse_options(usage,argc,argv);
    }
    int options_count() { return option_list_.size(); };
    int nonoption_count() {return nonoption_list_.size(); }
    bool error() { return err_; }
    std::string get_command_string() const { return command_string_; }
    bool has_option(char c) const;
    bool has_option(const char *p) const;
    std::string get_value(char c) const;
    std::string get_value(const char *p) const;
    std::string get_nonoption(int n) const;
    std::string get_program_name() const { return prog_name_; }

  private:
    void input_command_string(int argc, const char* const* argv) {
      std::stringstream os;
      if (argc==0) { return; }
      os << argv[0];
      for (unsigned int i=1;i<argc;++i) {
	os << " " << argv[i];
      }
      command_string_=os.str();
    }
 
    void parse_options(std::vector<Descriptor> usage,int argc, const char* const* argv) {
      if (argv==NULL || argv[0]==NULL)
	return; // Nothing to parse, should not happen.
      prog_name_ = argv[0];
      for (unsigned int i=1;i<argc;++i) {
	const char *p = argv[i];
	unsigned int len = strlen(p);
	if (len>2 && i+1<argc && 0==strncmp(p,"--",2)) {
	  const char *q=p+2;
	  for(std::vector<struct Descriptor>::iterator it = usage.begin(); it != usage.end(); ++it) {
	    if (!strcmp(q,(it->longopt))) {
	      option_list_.push_back(Option(*it, argv[i+1]));
	    }
	  }
	  i++; // advance i since we read argv[i+1] already
	} else if (len>1 && i+1<argc && *p=='-') {
	  const char q=p[1];
	  for(std::vector<struct Descriptor>::iterator it = usage.begin(); it != usage.end(); ++it) {
	    if (q ==it->shortopt) {
	      option_list_.push_back(Option(*it, argv[i+1]));
	    }
	  }
	  i++;  // advance i since we read argv[i+1] already
	} else {
	  // An argument that was not preceded by a short or long option.
	  nonoption_list_.push_back(std::string(argv[i]));
	}
      }
    }
  };

  /**
   * @param c A short option
   * @return true if the Option is valid (Because it was initialized so using a Descriptor).
   */
  bool Parser::has_option(char c) const {
    for (std::vector<Option>::const_iterator it = option_list_.begin(); it != option_list_.end(); ++it) {
      if (c == it->get_shortopt())
	return true;
    }
    return false;
  }

  /**
   * @param p A long option
   * @return true if the Option is valid (Because it was initialized so using a Descriptor).
   */
  bool Parser::has_option(const char *p) const {
    std::string longopt(p);
    for (std::vector<Option>::const_iterator it = option_list_.begin(); it != option_list_.end(); ++it) {
      if (longopt == it->get_longopt())
	return true;
    }
    return false;
  }

  /**
   * Return the value of the argument corresponding to the short option c
   * @param c The short option.
   */
  std::string Parser::get_value(char c) const {
    for (std::vector<Option>::const_iterator it = option_list_.begin(); it != option_list_.end(); ++it) {
      if (c == it->get_shortopt())
	return it->get_value();
    }
    return NULL;
  }

   /**
   * Return the value of the argument corresponding to the short option c
   * @param l The long option.
   */
  std::string Parser::get_value(const char* l) const {
    std::string longopt(l);
    for (std::vector<Option>::const_iterator it = option_list_.begin(); it != option_list_.end(); ++it) {
      if (longopt == it->get_longopt())
	return it->get_value();
    }
    return NULL;
  }

  /**
   * @param i index of the nonoption
   * Return the corresponding nonoption (argument not preceeded by - or -- option)
   */
  std::string Parser::get_nonoption(int i) const {
    if (i>= nonoption_list_.size())
      return "";
    return nonoption_list_[i];
  }
  
}
// namespace option

#endif /* OPTIONPARSER_H_ */
