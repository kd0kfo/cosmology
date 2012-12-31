/**
 * 
 * This file is part of physcalc, an interactive calculator utility
 * designed to carry out lensing calculations and to manipulate plane
 * data structures.
 *
 * This file contains parser code. It is meant to be used with
 * yacc (specificly GNU bison). For an excellent "A Compact Guide
 * to Lex & Yacc" by Tom Niemann" available at epaperpress.com. 
 * Portions of source code from that paper are used here with 
 * permission given in the Preface (page 3).
 * 
 * Copyright 2007, 2010 David Coss, PhD
 *
 * physcalc is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * physcalc is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with physcalc.  If not, see <http://www.gnu.org/licenses/>.
 */
#include "functions.h"
#include "symrec.h"
#include "physcalc.yacc.h"
#include "physcalc.h"
#include "version.h"

#include "libmygl/version.h"

#include "libdnstd/version.h"

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <cstdio>
#include <errno.h>
#include <dirent.h>
#include <getopt.h>
#include <limits.h>
#include <cmath>

const char* PROMPT_STRING = ">";
const char* program_name = "physcalc";
struct calcval PHYSCALC_ans;
extern int yyparse ();
extern char *get_radix();

int PHYSCALC_is_interactive;

struct init
{
  char const *fname;
  double (*fnct) (double);
  void (*v_c_fnct)(const char*);
};
     
struct init const arith_fncts[] =
  {
    "abs",  fabs,0,
    "acos", acos,0,
    "asin", asin,0,
    "atan", atan,0,
    "ceil", ceil,0,
    "cos",  cos,0,
    "cosh",  cosh,0,
    "exp",  exp,0,
    "floor", floor,0,
    "sin",  sin,0,
    "sinh",  sinh,0,
    "sqrt", sqrt,0,
    "tan",  tan,0,
    "tanh",  tanh,0,
    "log",  log10,0,
    "ln",   log,0,
    0, 0, 0
  };

struct init const other_fncts[] = 
  {
    "ls",0,do_ls,
    0,0,0
  };
     
/* The symbol table: a chain of `struct symrec'.  */
symrec *sym_table;
     
/**
 * Initialize function table
 */
void
init_table (void)
{
  int i;
  symrec *ptr;
  const struct init *curr = arith_fncts;

  i = 0;
  while (curr[i].fname != 0)
    {
      ptr = putsym (curr[i].fname, FNCT);
      if(curr[i].fnct != 0)
	ptr->value.fnctptr =	 curr[i].fnct;
      else if(curr[i].v_c_fnct != 0)
	ptr->value.v_c_fnctptr =	 curr[i].v_c_fnct;
      ptr->plane_fnctptr = NULL;

      i++;
      if(curr == arith_fncts && curr[i].fname == 0)
	{
	  curr = other_fncts;i = 0;
	}
    }

  ptr = putsym("print",FNCT);ptr->plane_fnctptr = print_plane;ptr->isPlaneFunction = true;
  ptr = putsym("clear",FNCT);ptr->plane_fnctptr = clear_plane;ptr->isPlaneFunction = true;
  ptr = putsym("add",FNCT);ptr->plane_fnctptr = add_planes;ptr->isPlaneFunction = true;
  ptr = putsym("subtract",FNCT);ptr->plane_fnctptr = subtract_planes;ptr->isPlaneFunction = true;
  ptr = putsym("multiply",FNCT);ptr->plane_fnctptr = multiply_planes;ptr->isPlaneFunction = true;
  ptr = putsym("open",FNCT);ptr->plane_fnctptr = open_plane;ptr->isPlaneFunction = true;
  ptr = putsym("save",FNCT);ptr->plane_fnctptr = save_plane;ptr->isPlaneFunction = true;
  ptr = putsym("copy",FNCT);ptr->plane_fnctptr = copy_plane;ptr->isPlaneFunction = true;
  ptr = putsym("fourier",FNCT);ptr->plane_fnctptr = fourier_plane;ptr->isPlaneFunction = true;
  ptr = putsym("ifourier",FNCT);ptr->plane_fnctptr = ifourier_plane;ptr->isPlaneFunction = true;
}
    

void do_funct(symrec *rec, struct calcval param, struct calcval *result)
{
  if(rec == NULL || rec->plane_fnctptr != NULL)
    return;

  result->re = (*rec->value.fnctptr)(param.re);
  result->im = 0;
}

void do_funct(symrec *rec, symrec *param, struct calcval *result)
{
  if(rec == NULL || param == NULL)
    return;

  if(rec->plane_fnctptr != 0)
    {
      symrec* newrec = (*rec->plane_fnctptr)(&param,1);
      if(newrec != NULL)
	{
	  newrec->next = (struct symrec *)sym_table;
	  sym_table = newrec;
	}
      return;
    }

  if(rec->value.v_c_fnctptr != 0 && rec->isPlaneFunction)
    {
      (*rec->value.v_c_fnctptr)(param->name);
      return;
    }
  result->re = (*rec->value.fnctptr)(param->value.var[0]);
  result->im = 0.0;
}


void do_funct2(symrec* rec, symrec* param, symrec* param2, struct calcval *result)
{
  if(rec == NULL || param == NULL || param2 == NULL)
    return;
  if(rec->plane_fnctptr == NULL ||  !rec->isPlaneFunction)
    return;

  symrec* params[] = {param,param2};
  symrec* newrec = (*rec->plane_fnctptr)(params,2);
  if(newrec != NULL)
    {
      newrec->next = (struct symrec *)sym_table;
      sym_table = newrec;
    }

}

void handle_plane(symrec *rec,double& i, double& j, struct calcval *result)
{
  rec->isPlane = true;
  if(rec->value.planeptr == NULL)
    {
      rec->value.planeptr = create_plane(i,j); 
      return;
    }
  try{
    const Double& complex = rec->value.planeptr->getValue((int)i,(int)j);
    result->re = complex.getValue(0);
    result->im = complex.getValue(1);
  }
  catch(DavidException de)
    {
      std::cout << de.what() << std::endl;
    }
}

void print_complex(struct calcval cnumber)
{
  double intpart = 1, absval = fabs(cnumber.re);
  
  if(PHYSCALC_is_interactive)
    printf("\tans: ");
  if(cnumber.im != 0)
    {
      math::Complex buff(cnumber.re,cnumber.im);
      printf("%s",buff.str().c_str());
    }
  else if(modf(cnumber.re,&intpart) == 0 && absval < 1.0e+9)
    printf(get_radix(),(long int)cnumber.re);
  else if(absval > 1.0e+9 || absval < 1.0e-5)
    printf("%e",cnumber.re);
  else
    printf("%f",cnumber.re);
  
  if(PHYSCALC_is_interactive)
    printf("\n%s",PROMPT_STRING);
  fflush(stdout);
  
}

void do_funct(symrec* fnct,  struct calcval *result)
{
  if(fnct == NULL)
    return;
  if(fnct->type != FNCT)
    return;
  if(fnct->value.v_c_fnctptr == 0)
    return;
  (*fnct->value.v_c_fnctptr)(NULL);
}

void print_help(const char *keyword)
{
  if(keyword == NULL || strcmp(keyword,"print") == 0)
    printf("print -- summarizes data contined in the specified variable.\n");
  if(keyword == NULL || strcmp(keyword,"clear") == 0)
    printf("clear -- clears plane variable contents. \n");
  if(keyword == NULL || strcmp(keyword,"add") == 0)
    printf("add -- adds two planes. \n\tUsage: add(A,B)\n\t\t A <- A+B\n");
  if(keyword == NULL || strcmp(keyword,"substract") == 0)
    printf("substracts -- subtracts two planes. \n\tUsage: subtract(A,B)\n\t\t A <- A-B\n");
  if(keyword == NULL || strcmp(keyword,"multiply") == 0)
    printf("add -- multiply two planes. \n\tUsage: multiply(A,B)\n\t\t A <- A*B\n");
  if(keyword == NULL || strcmp(keyword,"open") == 0)
    printf("open -- opens a plane and stores it in the specified variable.\n\tUsage: open(A,filename);\n");
  if(keyword == NULL || strcmp(keyword,"save") == 0)
    printf("save -- saves a plane save the specified file.\n\tUsage: save(A,filename);\n");
  if(keyword == NULL || strcmp(keyword,"copy") == 0)
    printf("copy -- copies a plane.\n\tUsage: copy(A,B)\n\t\t B <- A\n");
  if(keyword == NULL || strcmp(keyword,"fourier") == 0)
    printf("fourier -- fourier transforms a plane. Performed in place.\nUsage: fourier(A)\n\t\tA <- F(A)\n");
  if(keyword == NULL || strcmp(keyword,"ifourier") == 0)
    printf("ifourier -- inverse fourier transforms a plane. Performed in place.\nUsage: ifourier(A)\n\t\tA <- F^-1 (A)\n");
  if(keyword == NULL || strcmp(keyword,"radix") == 0)
    printf("radix -- set or print the current radix being used for display. Does not affect input.\nOptions: dec, oct, hex\n");
  if(keyword == NULL)
    {
      printf("\nThe following standard functions are also included:\n"); 
      const struct init *curr = arith_fncts;
      while(curr != NULL && curr->fname != 0)
	{
	  printf("%s",curr->fname);
	  curr = curr++;
	  if(curr != NULL && curr->fname != 0)
	    printf(",");
	}
      printf("\n");
    }
}

void do_ls(const char *path)
{
  DIR *dp;
  struct dirent *ep;
	
  if(path == NULL)
    dp = opendir ("./");
  else
    dp = opendir(path);
  if (dp != NULL)
    {
      while (ep = readdir (dp))
	if(ep->d_name && ep->d_name[0] != '.')
	  puts (ep->d_name);
      closedir (dp);
    }
  else
    printf ("Couldn't open the directory\n");

	
}

void print_copyright()
{
  printf("%s v%s  Copyright (C) 2012 David Coss, PhD\n",program_name,VERSION);
  printf("This program comes with ABSOLUTELY NO WARRANTY. Use and redistribution rights are granted under the terms of the GNU General Public License; for details visit http://www.gnu.org/copyleft/gpl.html or read the COPYING file distributed with this program.\n");
}

static const char short_opts[] = "c:f:hv";
enum{PHYSCALC_ARG_TAB_LONG,PHYSCALC_ARG_TAB_SHORT,PHYSCALC_ARG_TAB_ALL,PHYSCALC_BUILD_INFO};
struct option long_opts[] =
  {
    {"calc",1,NULL,'c'},
    {"file",1,NULL,'f'},
    {"help",0,NULL,'h'},
    {"tab-long",0,NULL,PHYSCALC_ARG_TAB_LONG},
    {"tab-short",0,NULL,PHYSCALC_ARG_TAB_SHORT},
    {"tab-all",0,NULL,PHYSCALC_ARG_TAB_ALL},
    {"version",0,NULL,'v'},
    {"build",0,NULL,PHYSCALC_BUILD_INFO},
    {NULL,0,NULL,0}
  };
 
void print_help()
{
  size_t curr = 0;
  printf("\nUsage: ./%s [options]\n",program_name);
  printf("\nOptions:\n\n");
  while(long_opts[curr].name != NULL)
    {
      printf("--%s",long_opts[curr].name);
      if(long_opts[curr].val >= 'a' && long_opts[curr].val <= 'z')
	printf(", -%c",long_opts[curr].val);
      printf("\t");
      switch(long_opts[curr].val)
	{
	case 'c':
	  printf("Calculate the specified expression");
	  break;
	case 'f':
	  printf("Specify a file to be used as input");
	  break;
	case 'h':
	  printf("Display this message");
	  break;
	case 'v':
	  printf("Display the version of %s",program_name);
	  break;
	case PHYSCALC_BUILD_INFO:
	  printf("Display build information");
	  break;
	default:
	  break;
	}
      printf("\n");
      curr++;
    }

  std::cout << std::endl;
}

void print_tab_complete(int short_or_long)
{
  struct option *the_opts = long_opts;
  
  while(the_opts != NULL && the_opts->name != NULL)
    {
      if(short_or_long == PHYSCALC_ARG_TAB_LONG)
	printf("--%s ",the_opts->name);
      else if(the_opts->val >= 'a' && the_opts->val <= 'z')
	printf("-%c ",the_opts->val);

      if(short_or_long == PHYSCALC_ARG_TAB_ALL)
	printf("--%s ",the_opts->name);
      
      the_opts++;
    }
}

void print_build_info()
{
  printf("Git Commit: %s\n",build_git_sha);
  printf("Build Time: %s\n\n",build_git_time);
  
  printf("mygl commit: %s\n",libmygl::build_git_sha);
  printf("mygl Build Time: %s\n\n",libmygl::build_git_time);

  printf("dnstd commit: %s\n",libdnstd::build_git_sha);
  printf("dnstd Build Time: %s\n\n",libdnstd::build_git_time);

}

int main (int argc, char **argv)
{
  symrec *funct;
  FILE *input = NULL;
  int optflag, weblength;
  char *strweblength;
  PHYSCALC_is_interactive = 0;
  init_table();

#ifdef WEB_USAGE
  printf("Content-Type: text/plain;charset=us-ascii\n\n");
  if((strweblength = getenv("CONTENT_LENGTH")) != NULL)
    {
      if(sscanf(strweblength,"%d",&weblength) == 1)
	{
	  char webbuff[(weblength + 1)*sizeof(char)];
	  memset(webbuff,0,sizeof(char)*(weblength+1));
	  input = tmpfile();
	  if(input == NULL)
	    {
	      printf("Could not calculate.");
	      exit(0);
	    }
	  fread(webbuff,sizeof(char),weblength,stdin);
	  fwrite(webbuff,sizeof(char),weblength,input);
	  rewind(input);
	}
    }
#endif

  while((optflag = getopt_long(argc,argv,short_opts,long_opts,NULL)) != -1)
    {
      switch(optflag)
	{
	case 'c':
	  input = tmpfile();
	  fputs(optarg,input);
	  rewind(input);
	  break;
	case 'f':
	  input = fopen(optarg,"r");
	  if(input == NULL)
	    {
	      fprintf(stderr,"Could not open %s for reading.\n",optarg);
	      fprintf(stderr,"Reason: %s\n",strerror(errno));
	      exit(errno);
	    }
	  break;
	case PHYSCALC_ARG_TAB_LONG: case PHYSCALC_ARG_TAB_SHORT:case PHYSCALC_ARG_TAB_ALL:
	  print_tab_complete(optflag);
	  exit(0);
	  break;
	case PHYSCALC_BUILD_INFO: case 'v':
	  print_copyright();
	  if(optflag == PHYSCALC_BUILD_INFO)
	    {
	      printf("\n");
	      print_build_info();
	    }
	  exit(0);
	case 'h':default:
	  print_copyright();
	  print_help();
	  exit(0);	  
	}
    }

  if(input == NULL)
    {
      print_copyright();
      printf("\n");
      printf("Welcome! To exit, type quit.\nALL statements must end with a semicolon.\n%s",PROMPT_STRING);
      PHYSCALC_is_interactive = 1;
      yyparse();
      printf("Good bye.\n");
    }
  else
    {
      yyrestart(input);
      yyparse();
      return PHYSCALC_ans.re;
    }
  return 0;
}



