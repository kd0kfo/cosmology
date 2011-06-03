#include <cstdio>
#include <dirent.h>
#include "functions.h"
#include "symrec.h"
#include "physcalc.yacc.h"
#include "physcalc.h"

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

const char* PROMPT_STRING = ">";

extern int yyparse ();

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
  using namespace std;
  cout << "\tans: ";
  if(cnumber.im != 0)
    {
      math::Complex buff(cnumber.re,cnumber.im);
      cout << buff;
    }
  else
    cout << cnumber.re;
  cout <<   endl << ">";
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
  printf("physcalc v%s  Copyright (C) 2010 David Coss, PhD\n",PACKAGE_VERSION);
  printf("This program comes with ABSOLUTELY NO WARRANTY; for details visit http://www.gnu.org/licenses/gpl-2.0.html or read the COPYING file distributed with this program.\n");
}

 
int
main (int argc, char **argv)
{
  symrec *funct;
  FILE *input = NULL;
  print_copyright();
  init_table ();
  if(argc > 1)
    {
      input = fopen(argv[1],"r");
      if(input != NULL)
	yyrestart(input);
    }

  if(input == NULL)
    {
      printf("\n");
      printf("Welcome! To exit, type quit.\nALL statements must end with a semicolon.\n%s",PROMPT_STRING);
      yyparse();
      printf("Good bye.\n");
    }
  else
    {
      yyparse();      
    }
  return 0;
}



