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


namespace option
{


  class Arg;

  class CheckArg{};

  enum class ArgType { NONE, STRING, INTEGER, FLOAT };
  
  //typedef ArgStatus(* CheckArg)(const Option &option, bool msg)

  struct Descriptor {
    /**
     * @brief Index of this option's linked list in the array filled in by the parser.
     *
     * @par Tip:
     * Use an enum rather than plain ints for better readability, as shown in the example
     * at Descriptor.
     */
    const unsigned index;
    
    /**
     * @brief Used to distinguish between options with the same @ref index.
     * See @ref index for details.
     *
     * It is recommended that you use an enum rather than a plain int to make your
     * code more readable.
     */
    const int type;
    
    /**
     * @brief Each char in this string will be accepted as a short option character.
     *
     * The string must not include the minus character @c '-' or you'll get undefined
     * behaviour.
     *
     * If this Descriptor should not have short option characters, use the empty
     * string "". NULL is not permitted here!
     *
     * See @ref longopt for more information.
     */
    const char* const shortopt;
    
    /**
     * @brief The long option name (without the leading @c -- ).
     *
     * If this Descriptor should not have a long option name, use the empty
     * string "". NULL is not permitted here!
     *
     * While @ref shortopt allows multiple short option characters, each
     * Descriptor can have only a single long option name. If you have multiple
     * long option names referring to the same option use separate Descriptors
     * that have the same @ref index and @ref type. You may repeat
     * short option characters in such an alias Descriptor but there's no need to.
     *
     * @par Dummy Descriptors:
     * You can use dummy Descriptors with an
     * empty string for both @ref shortopt and @ref longopt to add text to
     * the usage that is not related to a specific option. See @ref help.
     * The first dummy Descriptor will be used for unknown options (see below).
     *
     * @par Unknown Option Descriptor:
     * The first dummy Descriptor in the list of Descriptors,
     * whose @ref shortopt and @ref longopt are both the empty string, will be used
     * as the Descriptor for unknown options. An unknown option is a string in
     * the argument vector that is not a lone minus @c '-' but starts with a minus
     * character and does not match any Descriptor's @ref shortopt or @ref longopt. @n
     * Note that the dummy descriptor's @ref check_arg function @e will be called and
     * its return value will be evaluated as usual. I.e. if it returns @ref ARG_ILLEGAL
     * the parsing will be aborted with <code>Parser::error()==true</code>. @n
     * if @c check_arg does not return @ref ARG_ILLEGAL the descriptor's
     * @ref index @e will be used to pick the linked list into which
     * to put the unknown option. @n
     * If there is no dummy descriptor, unknown options will be dropped silently.
     *
     */
    const char* const longopt;
    
    /**
     * @brief For each option that matches @ref shortopt or @ref longopt this function
     * will be called to check a potential argument to the option.
     *
     * This function will be called even if there is no potential argument. In that case
     * it will be passed @c NULL as @c arg parameter. Do not confuse this with the empty
     * string.
     *
     * See @ref CheckArg for more information.
     */
    const  ArgType check_arg;
    
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
     *
     * @attention
     * Must be UTF-8-encoded. If your compiler supports C++11 you can use the "u8"
     * prefix to make sure string literals are properly encoded.
     */
    const char* help;
  };



 

  class Arg {
    
  public :
    enum ArgType type;

  };



  
}
// namespace option

#endif /* OPTIONPARSER_H_ */
