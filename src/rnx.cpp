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
 * @brief This is the only file required to use The Lean Mean C++ Option Parser.
 *        Just \#include it and you're set.
 *
 * The Lean Mean C++ Option Parser handles the program's command line arguments 
 * (argc, argv).
 * It supports the short and long option formats of getopt(), getopt_long() 
 * and getopt_long_only() but has a more convenient interface.
 * The following features set it apart from other option parsers:
 *
 * @par Highlights:
 * <ul style="padding-left:1em;margin-left:0">
 * <li> It is a header-only library. Just <code>\#include "optionparser.h"</code> and you're set.
 * <li> It is freestanding. There are no dependencies whatsoever, not even the
 *      C or C++ standard library.
 * <li> It has a usage message formatter that supports column alignment and
 *      line wrapping. This aids localization because it adapts to
 *      translated strings that are shorter or longer (even if they contain
 *      Asian wide characters).
 * <li> Unlike getopt() and derivatives it doesn't force you to loop through
 *     options sequentially. Instead you can access options directly like this:
 *     <ul style="margin-top:.5em">
 *     <li> Test for presence of a switch in the argument vector:
 *      @code if ( options[QUIET] ) ... @endcode
 *     <li> Evaluate --enable-foo/--disable-foo pair where the last one used wins:
 *     @code if ( options[FOO].last()->type() == DISABLE ) ... @endcode
 *     <li> Cumulative option (-v verbose, -vv more verbose, -vvv even more verbose):
 *     @code int verbosity = options[VERBOSE].count(); @endcode
 *     <li> Iterate over all --file=&lt;fname> arguments:
 *     @code for (Option* opt = options[FILE]; opt; opt = opt->next())
 *   fname = opt->arg; ... @endcode
 *     <li> If you really want to, you can still process all arguments in order:
 *     @code
 *   for (int i = 0; i < p.optionsCount(); ++i) {
 *     Option& opt = buffer[i];
 *     switch(opt.index()) {
 *       case HELP:    ...
 *       case VERBOSE: ...
 *       case FILE:    fname = opt.arg; ...
 *       case UNKNOWN: ...
 *     @endcode
 *     </ul>
 * </ul> @n
 * Despite these features the code size remains tiny. 
 * It is smaller than <a href="http://uclibc.org">uClibc</a>'s GNU getopt() and just a
 * couple 100 bytes larger than uClibc's SUSv3 getopt(). @n
 * (This does not include the usage formatter, of course. But you don't have to use that.)
 *
 * @par Download:
 * Tarball with examples and test programs:
 * <a style="font-size:larger;font-weight:bold" href="http://sourceforge.net/projects/optionparser/files/optionparser-1.3.tar.gz/download">optionparser-1.3.tar.gz</a> @n
 * Just the header (this is all you really need):
 * <a style="font-size:larger;font-weight:bold" href="http://optionparser.sourceforge.net/optionparser.h">optionparser.h</a>
 *
 * @par Changelog:
 * <b>Version 1.3:</b> Compatible with Microsoft Visual C++. @n
 * <b>Version 1.2:</b> Added @ref option::Option::namelen "Option::namelen" and removed the extraction
 *                     of short option characters into a special buffer. @n
 *                     Changed @ref option::Arg::Optional "Arg::Optional" to accept arguments if they are attached
 *                     rather than separate. This is what GNU getopt() does and how POSIX recommends
 *                     utilities should interpret their arguments.@n
 * <b>Version 1.1:</b> Optional mode with argument reordering as done by GNU getopt(), so that
 *                     options and non-options can be mixed. See
 *                     @ref option::Parser::parse() "Parser::parse()".
 *
 * @par Feedback:
 * Send questions, bug reports, feature requests etc. to: <tt><b>optionparser-feedback<span id="antispam">&nbsp;(a)&nbsp;</span>lists.sourceforge.net</b></tt>
 * @htmlonly <script type="text/javascript">document.getElementById("antispam").innerHTML="@"</script> @endhtmlonly
 *
 *
 * @par Example program:
 * (Note: @c option::* identifiers are links that take you to their documentation.)
 * @code
 * #error EXAMPLE SHORTENED FOR READABILITY. BETTER EXAMPLES ARE IN THE .TAR.GZ!
 * #include <iostream>
 * #include "optionparser.h"
 *
 * enum  optionIndex { UNKNOWN, HELP, PLUS };
 * const option::Descriptor usage[] =
 * {
 *  {UNKNOWN, 0,"" , ""    ,option::Arg::None, "USAGE: example [options]\n\n"
 *                                             "Options:" },
 *  {HELP,    0,"" , "help",option::Arg::None, "  --help  \tPrint usage and exit." },
 *  {PLUS,    0,"p", "plus",option::Arg::None, "  --plus, -p  \tIncrement count." },
 *  {UNKNOWN, 0,"" ,  ""   ,option::Arg::None, "\nExamples:\n"
 *                                             "  example --unknown -- --this_is_no_option\n"
 *                                             "  example -unk --plus -ppp file1 file2\n" },
 *  {0,0,0,0,0,0}
 * };
 *
 * int main(int argc, char* argv[])
 * {
 *   argc-=(argc>0); argv+=(argc>0); // skip program name argv[0] if present
 *   option::Stats  stats(usage, argc, argv);
 *   option::Option options[stats.options_max], buffer[stats.buffer_max];
 *   option::Parser parse(usage, argc, argv, options, buffer);
 *
 *   if (parse.error())
 *     return 1;
 *
 *   if (options[HELP] || argc == 0) {
 *     option::printUsage(std::cout, usage);
 *     return 0;
 *   }
 *
 *   std::cout << "--plus count: " <<
 *     options[PLUS].count() << "\n";
 *
 *   for (option::Option* opt = options[UNKNOWN]; opt; opt = opt->next())
 *     std::cout << "Unknown option: " << opt->name << "\n";
 *
 *   for (int i = 0; i < parse.nonOptionsCount(); ++i)
 *     std::cout << "Non-option #" << i << ": " << parse.nonOption(i) << "\n";
 * }
 * @endcode
 *
 * @par Option syntax:
 * @li The Lean Mean C++ Option Parser follows POSIX <code>getopt()</code> conventions and supports
 *     GNU-style <code>getopt_long()</code> long options as well as Perl-style single-minus
 *     long options (<code>getopt_long_only()</code>).
 * @li short options have the format @c -X where @c X is any character that fits in a char.
 * @li short options can be grouped, i.e. <code>-X -Y</code> is equivalent to @c -XY.
 * @li a short option may take an argument either separate (<code>-X foo</code>) or
 *     attached (@c -Xfoo). You can make the parser accept the additional format @c -X=foo by
 *     registering @c X as a long option (in addition to being a short option) and
 *     enabling single-minus long options.
 * @li an argument-taking short option may be grouped if it is the last in the group, e.g.
 *     @c -ABCXfoo or <code> -ABCX foo </code> (@c foo is the argument to the @c -X option).
 * @li a lone minus character @c '-' is not treated as an option. It is customarily used where
 *     a file name is expected to refer to stdin or stdout.
 * @li long options have the format @c --option-name.
 * @li the option-name of a long option can be anything and include any characters.
 *     Even @c = characters will work, but don't do that.
 * @li [optional] long options may be abbreviated as long as the abbreviation is unambiguous.
 *     You can set a minimum length for abbreviations.
 * @li [optional] long options may begin with a single minus. The double minus form is always
 *     accepted, too.
 * @li a long option may take an argument either separate (<code> --option arg </code>) or
 *     attached (<code> --option=arg </code>). In the attached form the equals sign is mandatory.
 * @li an empty string can be passed as an attached long option argument: <code> --option-name= </code>.
 *     Note the distinction between an empty string as argument and no argument at all.
 * @li an empty string is permitted as separate argument to both long and short options.
 * @li Arguments to both short and long options may start with a @c '-' character. E.g.
 *     <code> -X-X </code>, <code>-X -X</code> or <code> --long-X=-X </code>. If @c -X
 *     and @c --long-X take an argument, that argument will be @c "-X" in all 3 cases.
 * @li If using the built-in @ref option::Arg::Optional "Arg::Optional", optional arguments must
 *     be attached.
 * @li the special option @c -- (i.e. without a name) terminates the list of
 *     options. Everything that follows is a non-option argument, even if it starts with
 *     a @c '-' character. The @c -- itself will not appear in the parse results.
 * @li the first argument that doesn't start with @c '-' or @c '--' and does not belong to
 *     a preceding argument-taking option, will terminate the option list and is the
 *     first non-option argument. All following command line arguments are treated as
 *     non-option arguments, even if they start with @c '-' . @n
 *     NOTE: This behaviour is mandated by POSIX, but GNU getopt() only honours this if it is
 *     explicitly requested (e.g. by setting POSIXLY_CORRECT). @n
 *     You can enable the GNU behaviour by passing @c true as first argument to
 *     e.g. @ref option::Parser::parse() "Parser::parse()".
 * @li Arguments that look like options (i.e. @c '-' followed by at least 1 character) but
 *     aren't, are NOT treated as non-option arguments. They are treated as unknown options and
 *     are collected into a list of unknown options for error reporting. @n
 *     This means that in order to pass a first non-option
 *     argument beginning with the minus character it is required to use the
 *     @c -- special option, e.g.
 *     @code
 *     program -x -- --strange-filename
 *     @endcode
 *     In this example, @c --strange-filename is a non-option argument. If the @c --
 *     were omitted, it would be treated as an unknown option. @n
 *     See @ref option::Descriptor::longopt for information on how to collect unknown options.
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
    {UNKNOWN, 0,"" , ""    ,option::Arg::None, "USAGE: rnx [options]\n\n"
     "Options:" },
    {HELP,    0,"" , "help",option::Arg::None, "  --help  \tPrint usage and exit." },
    {UNKNOWN, 0,"" ,  ""   ,option::Arg::None, "\nExamples:\n"
     "  rnx CTfile\n"
     "  rnx -unk --plus -ppp file1 file2\n" },
    {0,0,0,0,0,0}
  };

int main(int argc, char* argv[]) {

  std::cout << "rnx " << std::endl;
  argc-=(argc>0); argv+=(argc>0); // skip program name argv[0] if present
  option::Stats  stats(usage, argc, argv);
  option::Option options[stats.options_max], buffer[stats.buffer_max];
  option::Parser parse(usage, argc, argv, options, buffer);
  
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

  
  const char * fname = "./testdata/u3.ct";

  for (option::Option* opt = options[UNKNOWN]; opt; opt = opt->next())
     std::cout << "Unknown option: " << opt->name << "\n";
  
  for (int i = 0; i < parse.nonOptionsCount(); ++i)
    std::cout << "Non-option #" << i << ": " << parse.nonOption(i) << "\n";


  const char *dir = "./dat";
  Datatable dattab(dir);
  RNAStructure rnastruct(fname);

  dattab.efn2(&rnastruct, 1);

  const char *newout = "testout.ct";
  rnastruct.ctout (newout);
  return 0;
}
