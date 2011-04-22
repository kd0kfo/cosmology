#ifndef SYMREC_H
#define SYMREC_H 1

#include <vector>
#include <string>
#include "libmygl/plane.h"

struct calcval{double re,im;};
/* Function type.  */
     typedef double (*func_t) (double);
     
     /* Data type for links in the chain of symbols.  */
     struct symrec
     {
       typedef struct symrec* (*plane_func_t) (struct symrec **vars,size_t size);

       char *name;  /* name of symbol */
       int type;    /* type of symbol: either VAR or FNCT */
       union
       {
         double var[2];      /* value of a VAR */
         func_t fnctptr;  /* value of a FNCT */
	 void (*v_c_fnctptr)(const char*);
	 Plane<math::Complex>* planeptr;
       } value;
       plane_func_t plane_fnctptr;
       bool isPlane;
       struct symrec *next;  /* link field */
     };
     
     typedef struct symrec symrec;
     
     /* The symbol table: a chain of `struct symrec'.  */
     extern symrec *sym_table;
     
     symrec *putsym (char const *, int);
     symrec *getsym (char const *);

#endif

