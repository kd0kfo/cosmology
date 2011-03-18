#ifndef ARGS_PARSER
#define ARGS_PARSER

#include <getopt.h>
static char program_doc[] = "interpreter for cosmological calculations";
static char usage_doc[] = "<options>";
struct good_option{
	struct option getopts_options;
	std::string help_string,arg_type;
};

static struct good_option options[] = {
  {{"todo",no_argument,0,0},"Prints a list of things that need to be done to improve this program.",""},
  {{"version",no_argument,0,1},"Printer the current version",""},
  {{"run",required_argument,0,'c'},"Compiles a specific file.","FILE"},
  {{"output",required_argument,0,'o'},"Name of the output file.","FILE"},
  {{"help",no_argument,0,'h'},"This help dialog.",""},
  {{"verbose",optional_argument,0,'v'},"Run verbosely.","LEVEL"},
  {{0,0,0,0},"",""}  
};

struct assembler_arguments
{
  std::string source_filename,output_filename; 
  int verbosity;
};

void args_usage()
{
	using namespace std;
	struct good_option curr_opt;
	size_t opt_index = 0;
	cout << PROGRAM_NAME << "-- " << program_doc << endl;
	cout << "Usage: " << PROGRAM_NAME << usage_doc << endl;
	curr_opt = options[opt_index++];
	
	cout << endl;
	while(curr_opt.getopts_options.name != 0)
	{
		struct option& gnu_opt = curr_opt.getopts_options;
		cout << "--" << gnu_opt.name;
		if(gnu_opt.val > 0x20)
			cout << ", -" << (char) gnu_opt.val;
		cout << " : " << curr_opt.help_string;
		curr_opt = options[opt_index++];
		if(gnu_opt.has_arg == required_argument)cout << "[";
		if(gnu_opt.has_arg == optional_argument)cout << "<";
		if(gnu_opt.has_arg != no_argument)cout << curr_opt.arg_type;
		if(gnu_opt.has_arg == required_argument)cout << "]";
		if(gnu_opt.has_arg == optional_argument)cout << ">";
		cout << endl;
	}
	exit(0);
}

enum ARG_ERROR{NO_ERROR = 0,ARGS_UNKNOWN = 1};

void load_verbosity(const char* verbose_level,int& verbose_holder)
{
  if(verbose_level == 0)
    {
      verbose_holder = 0;
      return;
    }
  std::istringstream buff(verbose_level);
  buff >> verbose_holder;
  if(buff.fail())
    verbose_holder = 0;
}

static ARG_ERROR parse_opt(int key, struct good_option* the_options, int index,struct assembler_arguments& args)
{
	using namespace std;
	struct good_option option = the_options[index];

	switch(key)
	{
	case 't':
	  cout << Build::todo << endl;
	  args_usage();
	case 1:
	  cout << PROGRAM_NAME << " version: " << Build::getVersion() << endl;
	  cout << "Last Build: " << Build::getBuild() << endl;
	  args_usage();
	case 'v':
	  load_verbosity(optarg,args.verbosity);
	  break;
	case 'c':
	  args.source_filename =  optarg;
	  break;
	case 'o':
	  args.output_filename = optarg;
	  break;
	default:
	  return ARGS_UNKNOWN;
	}
	return NO_ERROR;
}


void parse_args(int argc, char** argv, struct assembler_arguments& args)
{
	int c,option_index,opt_result;
	struct option* long_opts;
	struct good_option curr_opt;
	size_t num_opts = 0;

	curr_opt = options[num_opts];
	while(curr_opt.getopts_options.name != 0)
		curr_opt = options[++num_opts];
	num_opts++;//for the null opt at the end.
	long_opts = new struct option[num_opts];
	for(size_t i = 0;i<num_opts;i++)
	  long_opts[i] = options[i].getopts_options;

	while(true)
	{
	  c = getopt_long(argc,argv,"tvc:o:",long_opts,&option_index);
	  if(c == -1)
	    break;//end of args
	  parse_opt(c,options,optind,args);
	}
}

#endif
